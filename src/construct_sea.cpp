//
//  construct_sea.cpp
//  GOTTA_v1.6.1
//
//  Created by Jian Wu on 2025/11/13.
//


// construct_sea.cpp
#include <Rcpp.h>
#include <cmath>
#include <vector>

// User-provided helpers (already implemented elsewhere)
#include "isinpolygon.h"
#include "compute_bd.h"

using namespace Rcpp;

// [[Rcpp::export]]
List construct_sea(const NumericMatrix& vertex, const IntegerMatrix& face) {
  // vertex: n x 2
  // face  : m x 3 (1-based)
  
  // Basic checks
  if (vertex.ncol() != 2) {
    stop("vertex must be an n x 2 matrix.");
  }
  if (face.ncol() != 3) {
    stop("face must be an m x 3 matrix.");
  }
  
  int n_vert = vertex.nrow();
  
  // ------------------------------------------------------------
  // 1. Boundary + normalize to unit disk
  // ------------------------------------------------------------
  // boundary vertex indices (1-based)
  IntegerVector bd = compute_bd(face);
  int n_bd = bd.size();
  if (n_bd < 3) {
    stop("Boundary has fewer than 3 vertices; cannot proceed.");
  }
  
  // z = vertex[,1] + i * vertex[,2]
  NumericVector zx(n_vert), zy(n_vert);
  for (int i = 0; i < n_vert; ++i) {
    zx[i] = vertex(i, 0);
    zy[i] = vertex(i, 1);
  }
  
  // z <- z - mean(z)
  double mean_x = 0.0, mean_y = 0.0;
  for (int i = 0; i < n_vert; ++i) {
    mean_x += zx[i];
    mean_y += zy[i];
  }
  mean_x /= n_vert;
  mean_y /= n_vert;
  for (int i = 0; i < n_vert; ++i) {
    zx[i] -= mean_x;
    zy[i] -= mean_y;
  }
  
  // z <- z / max(Mod(z))
  double max_mod = 0.0;
  for (int i = 0; i < n_vert; ++i) {
    double r = std::sqrt(zx[i] * zx[i] + zy[i] * zy[i]);
    if (r > max_mod) max_mod = r;
  }
  if (max_mod <= 0.0) {
    stop("All vertices are at the origin; cannot normalize.");
  }
  for (int i = 0; i < n_vert; ++i) {
    zx[i] /= max_mod;
    zy[i] /= max_mod;
  }
  
  // z_bdy <- z[bd]
  NumericVector zbx(n_bd), zby(n_bd);
  for (int i = 0; i < n_bd; ++i) {
    int idx = bd[i] - 1;  // convert 1-based to 0-based
    if (idx < 0 || idx >= n_vert) {
      stop("Boundary index out of range.");
    }
    zbx[i] = zx[idx];
    zby[i] = zy[idx];
  }
  
  // mean boundary edge length: mean(Mod(z_bdy - z_bdy[c(2:length,1)]))
  double edge_sum = 0.0;
  for (int i = 0; i < n_bd; ++i) {
    int j = (i + 1) % n_bd;
    double dx = zbx[i] - zbx[j];
    double dy = zby[i] - zby[j];
    edge_sum += std::sqrt(dx * dx + dy * dy);
  }
  double edge_size = edge_sum / n_bd;
  
  // shrink slightly: z <- z * (1 - 3*edge_size)
  double shrink = 1.0 - 3.0 * edge_size;
  for (int i = 0; i < n_vert; ++i) {
    zx[i] *= shrink;
    zy[i] *= shrink;
  }
  // update boundary values from shrunken z
  for (int i = 0; i < n_bd; ++i) {
    int idx = bd[i] - 1;
    zbx[i] = zx[idx];
    zby[i] = zy[idx];
  }
  
  // unit circle sampling: circle = exp(1i * seq(0.01, 2*pi, by=edge_size))
  std::vector<double> theta_vec;
  for (double t = 0.01; t <= 2.0 * M_PI + 1e-12; t += edge_size) {
    theta_vec.push_back(t);
  }
  int n_circle = static_cast<int>(theta_vec.size());
  NumericVector circx(n_circle), circy(n_circle);
  for (int i = 0; i < n_circle; ++i) {
    circx[i] = std::cos(theta_vec[i]);
    circy[i] = std::sin(theta_vec[i]);
  }
  
  // ------------------------------------------------------------
  // 2. Generate interior grid points (meshgrid style)
  // ------------------------------------------------------------
  // xs <- seq(-1.1, 1.1, by=edge_size)
  std::vector<double> xs_vec, ys_vec;
  for (double v = -1.1; v <= 1.1 + 1e-12; v += edge_size) {
    xs_vec.push_back(v);
  }
  for (double v = -1.1; v <= 1.1 + 1e-12; v += edge_size) {
    ys_vec.push_back(v);
  }
  int nx = static_cast<int>(xs_vec.size());
  int ny = static_cast<int>(ys_vec.size());
  
  // X, Y: ny x nx
  NumericMatrix X(ny, nx), Y(ny, nx);
  for (int j = 0; j < nx; ++j) {
    for (int i = 0; i < ny; ++i) {
      X(i, j) = xs_vec[j];
      Y(i, j) = ys_vec[i];
    }
  }
  
  // X(2:2:end, ) += edge_size/2
  for (int i = 1; i < ny; i += 2) { // R's 2,4,6,... => indices 1,3,5,... in 0-based
    for (int j = 0; j < nx; ++j) {
      X(i, j) += edge_size / 2.0;
    }
  }
  
  // Xv <- as.vector(X); Yv <- as.vector(Y)  (column-major)
  int n_grid = nx * ny;
  NumericVector Xv(n_grid), Yv(n_grid);
  int idx = 0;
  for (int j = 0; j < nx; ++j) {
    for (int i = 0; i < ny; ++i) {
      Xv[idx] = X(i, j);
      Yv[idx] = Y(i, j);
      ++idx;
    }
  }
  
  // ------------------------------------------------------------
  // 3. Choose new vertices between shrunken boundary and circle
  // ------------------------------------------------------------
  // cp for inner boundary: (1 + edge_size/2) * (Re, Im) of z_bdy
  NumericMatrix cp_in(n_bd, 2);
  double scale_in = 1.0 + edge_size / 2.0;
  for (int i = 0; i < n_bd; ++i) {
    cp_in(i, 0) = scale_in * zbx[i];
    cp_in(i, 1) = scale_in * zby[i];
  }
  
  // cp for outer circle: (1 - edge_size/2) * (Re, Im) of circle
  NumericMatrix cp_out(n_circle, 2);
  double scale_out = 1.0 - edge_size / 2.0;
  for (int i = 0; i < n_circle; ++i) {
    cp_out(i, 0) = scale_out * circx[i];
    cp_out(i, 1) = scale_out * circy[i];
  }
  
  // xy = cbind(Xv, Yv)
  NumericMatrix xy(n_grid, 2);
  for (int i = 0; i < n_grid; ++i) {
    xy(i, 0) = Xv[i];
    xy(i, 1) = Yv[i];
  }
  
  // insidez <- isinpolygon(cp_in, xy)
  // insidec <- isinpolygon(cp_out, xy)
  LogicalVector insidez = isinpolygon(cp_in, xy);
  LogicalVector insidec = isinpolygon(cp_out, xy);
  
  // newv <- which(insidec & !insidez)
  std::vector<int> newv_idx;
  newv_idx.reserve(n_grid);
  for (int i = 0; i < n_grid; ++i) {
    if (insidec[i] && !insidez[i]) {
      // store R-style 1-based index of this grid point
      newv_idx.push_back(i + 1);
    }
  }
  
  int n_new = static_cast<int>(newv_idx.size());
  NumericMatrix new_vertices(n_new, 2);
  for (int k = 0; k < n_new; ++k) {
    int gi = newv_idx[k] - 1; // back to 0-based
    new_vertices(k, 0) = Xv[gi];
    new_vertices(k, 1) = Yv[gi];
  }
  
  // ------------------------------------------------------------
  // 4. Triangulate annulus with RTriangle
  // ------------------------------------------------------------
  int n_z = n_vert;
  int n_disk = n_z + n_new + n_circle;
  
  // vertex_disk = rbind( cbind(Re(z), Im(z)),
  //                      new_vertices,
  //                      cbind(Re(circle), Im(circle)) )
  NumericMatrix vertex_disk(n_disk, 2);
  // original vertices
  for (int i = 0; i < n_z; ++i) {
    vertex_disk(i, 0) = zx[i];
    vertex_disk(i, 1) = zy[i];
  }
  // new_vertices
  for (int i = 0; i < n_new; ++i) {
    vertex_disk(n_z + i, 0) = new_vertices(i, 0);
    vertex_disk(n_z + i, 1) = new_vertices(i, 1);
  }
  // circle vertices
  for (int i = 0; i < n_circle; ++i) {
    vertex_disk(n_z + n_new + i, 0) = circx[i];
    vertex_disk(n_z + n_new + i, 1) = circy[i];
  }
  
  // circ_idx <- seq(n_z + n_new + 1, n_disk)
  IntegerVector circ_idx(n_circle);
  for (int i = 0; i < n_circle; ++i) {
    circ_idx[i] = n_z + n_new + 1 + i; // R-style 1-based
  }
  
  // S_outer: segments for outer circle
  IntegerMatrix S_outer(n_circle, 2);
  for (int i = 0; i < n_circle; ++i) {
    int a = circ_idx[i];
    int b = circ_idx[(i + 1) % n_circle];
    S_outer(i, 0) = a;
    S_outer(i, 1) = b;
  }
  
  // S_inner: segments for inner boundary (hole)
  IntegerMatrix S_inner(n_bd, 2);
  for (int i = 0; i < n_bd; ++i) {
    int a = bd[i];
    int b = bd[(i + 1) % n_bd];
    S_inner(i, 0) = a;
    S_inner(i, 1) = b;
  }
  
  // S = rbind(S_outer, S_inner)
  IntegerMatrix S(n_circle + n_bd, 2);
  for (int i = 0; i < n_circle; ++i) {
    S(i, 0) = S_outer(i, 0);
    S(i, 1) = S_outer(i, 1);
  }
  for (int i = 0; i < n_bd; ++i) {
    S(n_circle + i, 0) = S_inner(i, 0);
    S(n_circle + i, 1) = S_inner(i, 1);
  }
  
  // hole_pt <- cbind(mean(Re(z_bdy)), mean(Im(z_bdy)))
  double mean_bx = 0.0, mean_by = 0.0;
  for (int i = 0; i < n_bd; ++i) {
    mean_bx += zbx[i];
    mean_by += zby[i];
  }
  mean_bx /= n_bd;
  mean_by /= n_bd;
  NumericMatrix hole_pt(1, 2);
  hole_pt(0, 0) = mean_bx;
  hole_pt(0, 1) = mean_by;
  
  // Call RTriangle::pslg and RTriangle::triangulate from C++
  Environment RTriangle = Environment::namespace_env("RTriangle");
  Function pslg = RTriangle["pslg"];
  Function triangulate = RTriangle["triangulate"];
  
  SEXP pslg_obj = pslg(Named("P") = vertex_disk,
                       Named("S") = S,
                       Named("H") = hole_pt);
  
  List tri = triangulate(pslg_obj);
  if (tri.containsElementNamed("T") == FALSE) {
    stop("RTriangle::triangulate() did not return 'T'.");
  }
  SEXP Tsexp = tri["T"];
  if (Rf_isNull(Tsexp)) {
    stop("RTriangle::triangulate() returned NULL T.");
  }
  IntegerMatrix face_new(Tsexp);  // m_new x 3
  
  // ------------------------------------------------------------
  // 5. Remove triangles whose centroids lie inside original bd
  // ------------------------------------------------------------
  int m_new = face_new.nrow();
  NumericMatrix centroids(m_new, 2);
  for (int i = 0; i < m_new; ++i) {
    int a = face_new(i, 0) - 1; // 0-based
    int b = face_new(i, 1) - 1;
    int c = face_new(i, 2) - 1;
    double cx = (vertex_disk(a, 0) + vertex_disk(b, 0) + vertex_disk(c, 0)) / 3.0;
    double cy = (vertex_disk(a, 1) + vertex_disk(b, 1) + vertex_disk(c, 1)) / 3.0;
    centroids(i, 0) = cx;
    centroids(i, 1) = cy;
  }
  
  // cp for boundary: cbind(Re(z_bdy), Im(z_bdy))
  NumericMatrix cp_bd(n_bd, 2);
  for (int i = 0; i < n_bd; ++i) {
    cp_bd(i, 0) = zbx[i];
    cp_bd(i, 1) = zby[i];
  }
  
  LogicalVector IN2 = isinpolygon(cp_bd, centroids);
  
  // Filter faces where !IN2
  std::vector<int> keep_faces;
  keep_faces.reserve(m_new);
  for (int i = 0; i < m_new; ++i) {
    if (!IN2[i]) {
      keep_faces.push_back(i);
    }
  }
  int m_keep = static_cast<int>(keep_faces.size());
  IntegerMatrix face_new_filtered(m_keep, 3);
  for (int k = 0; k < m_keep; ++k) {
    int i = keep_faces[k];
    face_new_filtered(k, 0) = face_new(i, 0);
    face_new_filtered(k, 1) = face_new(i, 1);
    face_new_filtered(k, 2) = face_new(i, 2);
  }
  
  // face_disk <- rbind(face, face_new_filtered)
  int m_old = face.nrow();
  IntegerMatrix face_disk(m_old + m_keep, 3);
  for (int i = 0; i < m_old; ++i) {
    face_disk(i, 0) = face(i, 0);
    face_disk(i, 1) = face(i, 1);
    face_disk(i, 2) = face(i, 2);
  }
  for (int i = 0; i < m_keep; ++i) {
    face_disk(m_old + i, 0) = face_new_filtered(i, 0);
    face_disk(m_old + i, 1) = face_new_filtered(i, 1);
    face_disk(m_old + i, 2) = face_new_filtered(i, 2);
  }
  
  // ------------------------------------------------------------
  // 6–7. Reflection "sea" block is commented out in original R,
  // so we stop here and return vertex_disk and face_disk.
  // ------------------------------------------------------------
  
  return List::create(
    _["vertex"] = vertex_disk,
    _["face"]   = face_disk
  );
}
