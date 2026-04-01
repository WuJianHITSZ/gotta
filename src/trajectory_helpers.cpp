// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace {

bool has_name(const CharacterVector& names, const std::string& target) {
  for (int i = 0; i < names.size(); ++i) {
    if (names[i] == target) {
      return true;
    }
  }
  return false;
}

bool is_spatial_mesh(SEXP mesh) {
  CharacterVector classes = Rf_getAttrib(mesh, R_ClassSymbol);
  return has_name(classes, "SpatialMesh");
}

NumericMatrix as_matrix(SEXP x) {
  if (Rf_isMatrix(x)) {
    return as<NumericMatrix>(x);
  }
  if (Rf_isFrame(x)) {
    return as<NumericMatrix>(x);
  }
  if (Rf_isVector(x)) {
    NumericVector v = as<NumericVector>(x);
    NumericMatrix m(v.size(), 1);
    for (int i = 0; i < v.size(); ++i) {
      m(i, 0) = v[i];
    }
    return m;
  }
  stop("Unsupported feature/layout type");
}

NumericMatrix get_layout_matrix(SEXP mesh, const std::string& name) {
  S4 mesh_s4(mesh);
  List cartogram = as<List>(mesh_s4.slot("cartogram"));
  CharacterVector cart_names = cartogram.names();
  if (has_name(cart_names, name)) {
    List cart_entry = cartogram[name];
    if (!cart_entry.containsElementNamed("uv")) {
      stop("Cartogram '%s' missing uv entry", name);
    }
    return as_matrix(cart_entry["uv"]);
  }

  List layout = as<List>(mesh_s4.slot("layout"));
  CharacterVector layout_names = layout.names();
  if (has_name(layout_names, name)) {
    return as_matrix(layout[name]);
  }

  stop("Layout '%s' not found in mesh@cartogram or mesh@layout", name);
}

NumericMatrix correlate3x3(const NumericMatrix& img, const double kernel[3][3]) {
  int rows = img.nrow();
  int cols = img.ncol();
  NumericMatrix out(rows, cols);
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      double sum = 0.0;
      for (int kr = -1; kr <= 1; ++kr) {
        for (int kc = -1; kc <= 1; ++kc) {
          int rr = r + kr;
          int cc = c + kc;
          if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) {
            continue;
          }
          sum += kernel[kr + 1][kc + 1] * img(rr, cc);
        }
      }
      out(r, c) = sum;
    }
  }
  return out;
}

NumericMatrix gaussian_blur(const NumericMatrix& img, double sigma) {
  if (sigma <= 0.0) {
    return clone(img);
  }

  int radius = static_cast<int>(std::ceil(2.0 * sigma));
  if (radius < 1) {
    radius = 1;
  }
  int size = 2 * radius + 1;
  std::vector<double> kernel(size * size);
  double sum = 0.0;
  double denom = 2.0 * sigma * sigma;

  for (int r = -radius; r <= radius; ++r) {
    for (int c = -radius; c <= radius; ++c) {
      double val = std::exp(-(r * r + c * c) / denom);
      kernel[(r + radius) * size + (c + radius)] = val;
      sum += val;
    }
  }

  for (double &val : kernel) {
    val /= sum;
  }

  int rows = img.nrow();
  int cols = img.ncol();
  NumericMatrix out(rows, cols);

  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      double acc = 0.0;
      for (int kr = -radius; kr <= radius; ++kr) {
        for (int kc = -radius; kc <= radius; ++kc) {
          int rr = r + kr;
          int cc = c + kc;
          if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) {
            continue;
          }
          double kval = kernel[(kr + radius) * size + (kc + radius)];
          acc += kval * img(rr, cc);
        }
      }
      out(r, c) = acc;
    }
  }

  return out;
}

List gradient_spatial_impl(SEXP mesh, const std::string& feature_name, Nullable<double> sigma) {
  if (!is_spatial_mesh(mesh)) {
    stop("mesh must be a SpatialMesh object");
  }
  if (feature_name.empty()) {
    stop("feature.name must be a single character string");
  }

  S4 mesh_s4(mesh);
  DataFrame vertex_metadata = as<DataFrame>(mesh_s4.slot("vertex_metadata"));
  CharacterVector vm_names = vertex_metadata.names();

  SEXP feature_sexp = R_NilValue;
  if (has_name(vm_names, feature_name)) {
    feature_sexp = vertex_metadata[feature_name];
  } else {
    feature_sexp = get_layout_matrix(mesh, feature_name);
  }

  if (feature_sexp == R_NilValue) {
    stop("Feature or layout '%s' not found in mesh@vertex_metadata or mesh@layout", feature_name);
  }

  NumericMatrix feature = as_matrix(feature_sexp);
  NumericMatrix array_coords = get_layout_matrix(mesh, "array_coords");
  NumericMatrix layout = get_layout_matrix(mesh, "spatial_coords");
  IntegerVector boundary = mesh_s4.slot("boundary");

  if (feature.nrow() != array_coords.nrow()) {
    stop("Feature has %d rows but array_coords has %d rows. They must match.",
         feature.nrow(), array_coords.nrow());
  }

  double Dy[3][3] = {
    {0.0, 0.0, 0.0},
    {-0.5, 0.0, 0.5},
    {0.0, 0.0, 0.0}
  };
  double Dx[3][3] = {
    {-1.0 / (2.0 * std::sqrt(3.0)), -1.0 / (2.0 * std::sqrt(3.0)), 0.0},
    {0.0, 0.0, 0.0},
    {0.0, 1.0 / (2.0 * std::sqrt(3.0)), 1.0 / (2.0 * std::sqrt(3.0))}
  };

  int imsize1 = 0;
  int imsize2 = 0;
  for (int i = 0; i < array_coords.nrow(); ++i) {
    int v1 = static_cast<int>(std::round(array_coords(i, 0)));
    int v2 = static_cast<int>(std::round(array_coords(i, 1)));
    if (v1 > imsize1) imsize1 = v1;
    if (v2 > imsize2) imsize2 = v2;
  }

  int n_points = array_coords.nrow();
  int n_features = feature.ncol();
  if (n_features == 0) {
    n_features = 1;
  }

  NumericMatrix image_x_final(n_points, n_features);
  NumericMatrix image_y_final(n_points, n_features);

  for (int k = 0; k < n_features; ++k) {
    NumericMatrix img_mat(imsize1, imsize2);
    for (int i = 0; i < n_points; ++i) {
      int row = static_cast<int>(std::round(array_coords(i, 0))) - 1;
      int col = static_cast<int>(std::round(array_coords(i, 1))) - 1;
      if (row < 0 || row >= imsize1 || col < 0 || col >= imsize2) {
        continue;
      }
      img_mat(row, col) = feature(i, k);
    }

    NumericMatrix grad_x = correlate3x3(img_mat, Dx);
    NumericMatrix grad_y = correlate3x3(img_mat, Dy);

    if (sigma.isNotNull()) {
      double s = as<double>(sigma);
      if (s > 0) {
        grad_x = gaussian_blur(grad_x, s);
        grad_y = gaussian_blur(grad_y, s);
      }
    }

    for (int i = 0; i < n_points; ++i) {
      int row = static_cast<int>(std::round(array_coords(i, 0))) - 1;
      int col = static_cast<int>(std::round(array_coords(i, 1))) - 1;
      if (row < 0 || row >= imsize1 || col < 0 || col >= imsize2) {
        continue;
      }
      image_x_final(i, k) = grad_x(row, col);
      image_y_final(i, k) = grad_y(row, col);
    }
  }

  if (boundary.size() > 0) {
    for (int b = 0; b < boundary.size(); ++b) {
      int row = boundary[b] - 1;
      if (row < 0 || row >= n_points) {
        continue;
      }
      for (int k = 0; k < n_features; ++k) {
        image_x_final(row, k) = 0.0;
        image_y_final(row, k) = 0.0;
      }
    }
  }

  NumericMatrix trajectory(n_points, 2 * n_features);
  for (int i = 0; i < n_points; ++i) {
    for (int k = 0; k < n_features; ++k) {
      trajectory(i, k) = image_x_final(i, k);
      trajectory(i, k + n_features) = image_y_final(i, k);
    }
  }

  return List::create(
    _["trajectory"] = trajectory,
    _["layout"] = layout
  );
}

}  // namespace

// [[Rcpp::export]]
List gradient_spatial(SEXP mesh, std::string feature_name, Nullable<double> sigma = R_NilValue) {
  return gradient_spatial_impl(mesh, feature_name, sigma);
}

// [[Rcpp::export]]
List gradient_layout(SEXP mesh, std::string feature_name, std::string layout_name,
                     Nullable<double> sigma = R_NilValue) {
  if (!is_spatial_mesh(mesh)) {
    stop("mesh must be a SpatialMesh object");
  }
  if (feature_name.empty()) {
    stop("feature.name must be a single character string");
  }
  if (layout_name.empty()) {
    stop("layout.name must be a single character string");
  }

  S4 mesh_s4(mesh);
  DataFrame vertex_metadata = as<DataFrame>(mesh_s4.slot("vertex_metadata"));
  CharacterVector vm_names = vertex_metadata.names();
  if (!has_name(vm_names, feature_name)) {
    stop("Feature '%s' not found in mesh@vertex_metadata", feature_name);
  }

  SEXP feature_value = vertex_metadata[feature_name];
  List result_feature = gradient_spatial_impl(mesh, feature_name, sigma);
  NumericMatrix trajectory_feature = result_feature["trajectory"];

  List result_layout = gradient_spatial_impl(mesh, layout_name, sigma);
  NumericMatrix trajectory_layout = result_layout["trajectory"];

  int n_points = trajectory_layout.nrow();
  if (trajectory_feature.nrow() != n_points) {
    stop("Dimension mismatch: trajectory_feature has %d rows but trajectory_layout has %d rows",
         trajectory_feature.nrow(), n_points);
  }
  if (trajectory_layout.ncol() != 4) {
    stop("trajectory_layout must have 4 columns (got %d). Layout must be 2D.", trajectory_layout.ncol());
  }
  if (trajectory_feature.ncol() != 2) {
    stop("trajectory_feature must have 2 columns (got %d). Feature must be scalar.", trajectory_feature.ncol());
  }

  NumericMatrix trajectory(n_points, 2);
  NumericVector distortion(n_points);
  NumericVector divergence(n_points);
  NumericVector curl(n_points);

  const double eps = 1e-12;

  for (int i = 0; i < n_points; ++i) {
    double a11 = trajectory_layout(i, 0);
    double a12 = trajectory_layout(i, 1);
    double a21 = trajectory_layout(i, 2);
    double a22 = trajectory_layout(i, 3);
    double det = a11 * a22 - a12 * a21;

    if (std::abs(det) < eps) {
      trajectory(i, 0) = 0.0;
      trajectory(i, 1) = 0.0;
      distortion[i] = 0.0;
      divergence[i] = 0.0;
      curl[i] = 0.0;
      continue;
    }

    double b1 = trajectory_feature(i, 0);
    double b2 = trajectory_feature(i, 1);

    trajectory(i, 0) = (b1 * a22 - b2 * a12) / det;
    trajectory(i, 1) = (a11 * b2 - a21 * b1) / det;

    distortion[i] = det;
    divergence[i] = a11 + a22;
    curl[i] = a12 - a21;
  }

  IntegerVector boundary = mesh_s4.slot("boundary");
  if (boundary.size() > 0) {
    for (int b = 0; b < boundary.size(); ++b) {
      int row = boundary[b] - 1;
      if (row < 0 || row >= n_points) {
        continue;
      }
      trajectory(row, 0) = 0.0;
      trajectory(row, 1) = 0.0;
      distortion[row] = 0.0;
      divergence[row] = 0.0;
      curl[row] = 0.0;
    }
  }

  return List::create(
    _["trajectory"] = trajectory,
    _["distortion"] = distortion,
    _["layout"] = get_layout_matrix(mesh, layout_name),
    _["feature"] = feature_value,
    _["divergence"] = divergence,
    _["curl"] = curl
  );
}

// [[Rcpp::export]]
List jacobian_layout(SEXP mesh, std::string xy_name, std::string uv_name,
                     Nullable<double> sigma = R_NilValue) {
  if (!is_spatial_mesh(mesh)) {
    stop("mesh must be a SpatialMesh object");
  }

  get_layout_matrix(mesh, xy_name);
  get_layout_matrix(mesh, uv_name);

  List result_xy = gradient_spatial_impl(mesh, xy_name, sigma);
  NumericMatrix xy = result_xy["trajectory"];

  List result_uv = gradient_spatial_impl(mesh, uv_name, sigma);
  NumericMatrix uv = result_uv["trajectory"];

  int n_points = xy.nrow();
  NumericVector distortion(n_points);
  NumericVector divergence(n_points);
  NumericVector curl(n_points);

  const double eps = 1e-12;

  for (int i = 0; i < n_points; ++i) {
    double a11 = xy(i, 0);
    double a12 = xy(i, 1);
    double a21 = xy(i, 2);
    double a22 = xy(i, 3);
    double detA = a11 * a22 - a12 * a21;

    if (std::abs(detA) < eps) {
      distortion[i] = 0.0;
      divergence[i] = 0.0;
      curl[i] = 0.0;
      continue;
    }

    double b11 = uv(i, 0);
    double b12 = uv(i, 1);
    double b21 = uv(i, 2);
    double b22 = uv(i, 3);

    double invA11 = a22 / detA;
    double invA12 = -a12 / detA;
    double invA21 = -a21 / detA;
    double invA22 = a11 / detA;

    double j11 = invA11 * b11 + invA12 * b21;
    double j12 = invA11 * b12 + invA12 * b22;
    double j21 = invA21 * b11 + invA22 * b21;
    double j22 = invA21 * b12 + invA22 * b22;

    distortion[i] = j11 * j22 - j12 * j21;
    divergence[i] = j11 + j22;
    curl[i] = j12 - j21;
  }

  S4 mesh_s4(mesh);
  IntegerVector boundary = mesh_s4.slot("boundary");
  if (boundary.size() > 0) {
    for (int b = 0; b < boundary.size(); ++b) {
      int row = boundary[b] - 1;
      if (row < 0 || row >= n_points) {
        continue;
      }
      distortion[row] = 0.0;
      divergence[row] = 0.0;
      curl[row] = 0.0;
    }
  }

  return List::create(
    _["distortion"] = distortion,
    _["divergence"] = divergence,
    _["curl"] = curl
  );
}
