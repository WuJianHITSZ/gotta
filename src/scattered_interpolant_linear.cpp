// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// ===== Utilities =============================================================

inline bool any_na_row(const NumericMatrix& M, int i){
  for (int j = 0; j < M.ncol(); ++j)
    if (NumericMatrix::is_na(M(i,j))) return true;
  return false;
}

// Compute squared distance
inline double sqdist2(double x1, double y1, double x2, double y2){
  double dx = x1 - x2, dy = y1 - y2;
  return dx*dx + dy*dy;
}

// ===== Model builder: si_create_linear ======================================

/*
  Build a model for linear (barycentric) interpolation over a 2D triangulation.

  xy  : (n x 2) numeric matrix of sample locations
  z   : (n)     numeric vector of values
  tri : (m x 3) integer matrix of triangles (1-based indices)
  eps : determinant tolerance for degenerate triangles

  Returns a list with:
   - xy, z, tri (as given)
   - inv : (m x 4) row-major inverse of edge matrices per triangle
   - p3  : (m x 2) anchor vertex per triangle
   - valid : (m) logical; triangle usable?
   - centroid : (m x 2) centroids of triangles (for nearest-triangle fallback)
   - eps : stored tolerance
*/
// [[Rcpp::export]]
Rcpp::List si_create_linear(const Rcpp::NumericMatrix& xy,
                            const Rcpp::NumericVector& z,
                            const Rcpp::IntegerMatrix& tri,
                            const double eps = 1e-12)
{
  if (xy.ncol() != 2)
    stop("si_create_linear: 'xy' must be an n x 2 matrix.");
  if (xy.nrow() != z.size())
    stop("si_create_linear: 'xy' and 'z' must have the same number of rows/length.");
  if (tri.ncol() != 3)
    stop("si_create_linear: 'tri' must be an m x 3 matrix (1-based indices).");
  if (!R_finite(eps) || eps <= 0.0)
    stop("si_create_linear: 'eps' must be a finite positive number.");

  const int n = xy.nrow();
  const int m = tri.nrow();

  NumericMatrix inv(m, 4);     // [a11 a12 a21 a22]
  NumericMatrix p3 (m, 2);     // (x3, y3)
  LogicalVector valid(m, false);
  NumericMatrix centroid(m, 2);

  for (int t = 0; t < m; ++t) {
    int i1 = tri(t,0) - 1;
    int i2 = tri(t,1) - 1;
    int i3 = tri(t,2) - 1;

    if (i1 < 0 || i1 >= n || i2 < 0 || i2 >= n || i3 < 0 || i3 >= n) {
      continue; // invalid indices
    }

    double x1 = xy(i1,0), y1 = xy(i1,1);
    double x2 = xy(i2,0), y2 = xy(i2,1);
    double x3 = xy(i3,0), y3 = xy(i3,1);

    if (NumericMatrix::is_na(x1) || NumericMatrix::is_na(y1) ||
        NumericMatrix::is_na(x2) || NumericMatrix::is_na(y2) ||
        NumericMatrix::is_na(x3) || NumericMatrix::is_na(y3)) {
      continue; // NA coordinates -> invalid triangle
    }

    // store centroid (useful for fallback triangle)
    centroid(t,0) = (x1 + x2 + x3) / 3.0;
    centroid(t,1) = (y1 + y2 + y3) / 3.0;

    // anchor p3
    p3(t,0) = x3;
    p3(t,1) = y3;

    // Edge matrix columns: v1 = p1-p3, v2 = p2-p3
    double a = x1 - x3; // A = [[a b],[c d]]
    double c = y1 - y3;
    double b = x2 - x3;
    double d = y2 - y3;

    double det = a*d - b*c;
    if (!R_finite(det) || std::fabs(det) <= eps) {
      valid[t] = false;
      // inv left as zeros
      continue;
    }

    double inv_a11 =  d / det;
    double inv_a12 = -b / det;
    double inv_a21 = -c / det;
    double inv_a22 =  a / det;

    inv(t,0) = inv_a11;
    inv(t,1) = inv_a12;
    inv(t,2) = inv_a21;
    inv(t,3) = inv_a22;

    valid[t] = true;
  }

  return Rcpp::List::create(
    _["xy"]       = xy,
    _["z"]        = z,
    _["tri"]      = tri,
    _["inv"]      = inv,
    _["p3"]       = p3,
    _["valid"]    = valid,
    _["centroid"] = centroid,
    _["eps"]      = eps
  );
}

// ===== Predictor: si_predict_linear =========================================

/*
  Predict values at query points using the linear model:

  model : list from si_create_linear()
  query : either numeric(2) or (m x 2) numeric matrix
  inside_tol : tolerance for "inside triangle" test (allows small negatives)

  Behavior:
   - If a valid triangle contains the query (all barycentric >= -inside_tol),
     interpolate linearly by barycentric weights.
   - Otherwise, choose the valid triangle with nearest centroid and evaluate
     the same affine combination — linear extrapolation.

  Returns numeric vector of length m (or 1 for length-2 input).
*/
// [[Rcpp::export]]
Rcpp::NumericVector si_predict_linear(const Rcpp::List& model,
                                      SEXP query,
                                      double inside_tol = 1e-12)
{
  // Extract model parts
  if (!model.containsElementNamed("xy") ||
      !model.containsElementNamed("z")  ||
      !model.containsElementNamed("tri")||
      !model.containsElementNamed("inv")||
      !model.containsElementNamed("p3") ||
      !model.containsElementNamed("valid"))
    stop("si_predict_linear: 'model' must come from si_create_linear().");

  NumericMatrix xy       = model["xy"];
  NumericVector z        = model["z"];
  IntegerMatrix tri      = model["tri"];
  NumericMatrix inv      = model["inv"];
  NumericMatrix p3       = model["p3"];
  LogicalVector valid    = model["valid"];

  NumericMatrix centroid;
  if (model.containsElementNamed("centroid")) {
    // *** FIX: explicitly coerce to avoid ambiguous operator= ***
    centroid = Rcpp::as<Rcpp::NumericMatrix>(model["centroid"]);
  } else {
    // build centroids on the fly (shouldn't happen if created by si_create_linear)
    centroid = NumericMatrix(tri.nrow(), 2);
    for (int t=0;t<tri.nrow();++t){
      int i1 = tri(t,0)-1, i2 = tri(t,1)-1, i3 = tri(t,2)-1;
      centroid(t,0) = (xy(i1,0)+xy(i2,0)+xy(i3,0))/3.0;
      centroid(t,1) = (xy(i1,1)+xy(i2,1)+xy(i3,1))/3.0;
    }
  }

  // Coerce query into an m x 2 matrix
  NumericMatrix Q;
  bool single_point = false;
  if (Rf_isNumeric(query) && Rf_isMatrix(query)) {
    Q = Rcpp::as<Rcpp::NumericMatrix>(query);
    if (Q.ncol() != 2) stop("si_predict_linear: query matrix must have 2 columns.");
  } else if (Rf_isNumeric(query) && Rf_isVector(query)) {
    NumericVector v = Rcpp::as<Rcpp::NumericVector>(query);
    if (v.size() != 2) stop("si_predict_linear: numeric query must be length-2.");
    Q = NumericMatrix(1,2);
    Q(0,0) = v[0]; Q(0,1) = v[1];
    single_point = true;
  } else {
    stop("si_predict_linear: 'query' must be numeric(2) or an m x 2 numeric matrix.");
  }

  const int m = Q.nrow();
  NumericVector out(m, NA_REAL);

  for (int qi = 0; qi < m; ++qi) {
    if (any_na_row(Q, qi)) {
      out[qi] = NA_REAL;
      continue;
    }

    double xq = Q(qi,0), yq = Q(qi,1);

    // 1) Try to find a containing triangle (linear interpolation)
    bool done = false;
    for (int t = 0; t < tri.nrow(); ++t) {
      if (!valid[t]) continue;

      // Anchor and inverse for this triangle
      double x3 = p3(t,0), y3 = p3(t,1);
      double a11 = inv(t,0), a12 = inv(t,1);
      double a21 = inv(t,2), a22 = inv(t,3);

      // Compute l1,l2; l3 = 1 - l1 - l2
      double vx = xq - x3, vy = yq - y3;
      double l1 = a11*vx + a12*vy;
      double l2 = a21*vx + a22*vy;
      double l3 = 1.0 - l1 - l2;

      if (l1 >= -inside_tol && l2 >= -inside_tol && l3 >= -inside_tol) {
        int i1 = tri(t,0) - 1;
        int i2 = tri(t,1) - 1;
        int i3 = tri(t,2) - 1;
        out[qi] = l1*z[i1] + l2*z[i2] + l3*z[i3];
        done = true;
        break;
      }
    }
    if (done) continue;

    // 2) Outside all triangles: linear extrapolation via nearest triangle plane
    double best_d2 = R_PosInf;
    int best_t = -1;
    for (int t = 0; t < tri.nrow(); ++t) {
      if (!valid[t]) continue;
      double d2 = sqdist2(xq, yq, centroid(t,0), centroid(t,1));
      if (d2 < best_d2) { best_d2 = d2; best_t = t; }
    }

    if (best_t >= 0) {
      double x3 = p3(best_t,0), y3 = p3(best_t,1);
      double a11 = inv(best_t,0), a12 = inv(best_t,1);
      double a21 = inv(best_t,2), a22 = inv(best_t,3);

      double vx = xq - x3, vy = yq - y3;
      double l1 = a11*vx + a12*vy;
      double l2 = a21*vx + a22*vy;
      double l3 = 1.0 - l1 - l2;

      int i1 = tri(best_t,0) - 1;
      int i2 = tri(best_t,1) - 1;
      int i3 = tri(best_t,2) - 1;
      out[qi] = l1*z[i1] + l2*z[i2] + l3*z[i3];
    } else {
      // No valid triangles at all
      out[qi] = NA_REAL;
    }
  }

  return out;
}
