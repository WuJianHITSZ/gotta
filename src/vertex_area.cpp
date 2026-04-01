// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
using namespace Rcpp;

// ----- low-level helpers (column-major safe) -----

// squared distance between two vertex rows (1-based)
inline double sqdist_row(const NumericMatrix& V, int ai, int bi) {
  const int a = ai - 1, b = bi - 1;
  const int d = V.ncol();
  double s = 0.0;
  for (int k = 0; k < d; ++k) {
    const double diff = V(a, k) - V(b, k);
    s += diff * diff;
  }
  return s;
}

// dot of two edge vectors (vA - vB) · (vC - vD), rows are 1-based
inline double dot_edges(const NumericMatrix& V, int A, int B, int C, int D) {
  const int a = A - 1, b = B - 1, c = C - 1, d = D - 1;
  const int dim = V.ncol();
  double s = 0.0;
  for (int k = 0; k < dim; ++k) {
    const double u = V(a, k) - V(b, k);
    const double v = V(c, k) - V(d, k);
    s += u * v;
  }
  return s;
}

// clamp tiny negative to zero
inline double nonneg(double x) { return x < 0.0 ? 0.0 : x; }

// ----- face areas (Heron; works for any d >= 2) -----
// [[Rcpp::export]]
NumericVector face_area(const IntegerMatrix& face, const NumericMatrix& vertex) {
  const int nf = face.nrow();
  if (face.ncol() != 3) stop("face must be n x 3 (triangles).");
  if (vertex.ncol() < 2) stop("vertex must have dimension d >= 2.");

  NumericVector fa(nf);
  for (int f = 0; f < nf; ++f) {
    const int i = face(f, 0), j = face(f, 1), k = face(f, 2);
    const double a = std::sqrt(sqdist_row(vertex, j, i)); // |vj-vi|
    const double b = std::sqrt(sqdist_row(vertex, k, j)); // |vk-vj|
    const double c = std::sqrt(sqdist_row(vertex, i, k)); // |vi-vk|
    const double s = 0.5 * (a + b + c);
    fa[f] = std::sqrt(nonneg(s * (s - a) * (s - b) * (s - c)));
  }
  return fa;
}

// ----- main: per-vertex area (one_ring / mixed) -----
// [[Rcpp::export]]
NumericVector vertex_area(const IntegerMatrix& face,
                              const NumericMatrix& vertex,
                              std::string type = "one_ring") {
  const int nf = face.nrow();
  const int nv = vertex.nrow();

  if (face.ncol() != 3) stop("face must be n x 3 (triangles).");
  if (vertex.ncol() < 2) stop("vertex must have dimension d >= 2.");

  // normalize type
  std::transform(type.begin(), type.end(), type.begin(),
                 [](unsigned char c){ return std::tolower(c); });

  // per-face area
  NumericVector fa = face_area(face, vertex);

  NumericVector va(nv, 0.0);

  // -------- one_ring: sum incident face areas (not divided by 3) --------
  if (type == "one_ring") {
    for (int f = 0; f < nf; ++f) {
      const int i1 = face(f, 0), i2 = face(f, 1), i3 = face(f, 2);
      const double A = fa[f];
      va[i1 - 1] += A;
      va[i2 - 1] += A;
      va[i3 - 1] += A;
    }
    return va;
  }

  if (type != "mixed") {
    stop("vertex_area: type must be 'one_ring' or 'mixed'.");
  }

  // -------- mixed (Voronoi/mixed area with obtuse handling) -------------
  // Precompute edge squared lengths per face
  std::vector<double> l2_12(nf), l2_23(nf), l2_31(nf);
  for (int f = 0; f < nf; ++f) {
    const int i1 = face(f, 0), i2 = face(f, 1), i3 = face(f, 2);
    l2_12[f] = sqdist_row(vertex, i2, i1);
    l2_23[f] = sqdist_row(vertex, i3, i2);
    l2_31[f] = sqdist_row(vertex, i1, i3);
  }

  // Per-face cotangents at each corner:
  // c1 = cot(angle at i1) = cot( (i2 - i1), (i3 - i1) )
  // c2 = cot(angle at i2) = cot( (i3 - i2), (i1 - i2) )
  // c3 = cot(angle at i3) = cot( (i1 - i3), (i2 - i3) )
  const double eps = std::numeric_limits<double>::epsilon();
  std::vector<double> c1(nf), c2(nf), c3(nf);

  for (int f = 0; f < nf; ++f) {
    const int i1 = face(f, 0), i2 = face(f, 1), i3 = face(f, 2);

    // for c1
    {
      const double cs = dot_edges(vertex, i2, i1, i3, i1);
      const double d1 = l2_12[f];
      const double d2 = l2_31[f]; // |i1 - i3|^2 == |i3 - i1|^2
      const double rad = d1 * d2 - cs * cs;
      const double denom = std::sqrt(rad > 0.0 ? rad : 0.0);
      c1[f] = cs / std::max(denom, eps);
    }
    // for c2
    {
      const double cs = dot_edges(vertex, i3, i2, i1, i2);
      const double d1 = l2_23[f];
      const double d2 = l2_12[f];
      const double rad = d1 * d2 - cs * cs;
      const double denom = std::sqrt(rad > 0.0 ? rad : 0.0);
      c2[f] = cs / std::max(denom, eps);
    }
    // for c3
    {
      const double cs = dot_edges(vertex, i1, i3, i2, i3);
      const double d1 = l2_31[f];
      const double d2 = l2_23[f];
      const double rad = d1 * d2 - cs * cs;
      const double denom = std::sqrt(rad > 0.0 ? rad : 0.0);
      c3[f] = cs / std::max(denom, eps);
    }
  }

  // preliminary per-face-per-vertex areas
  std::vector<double> vaf1(nf), vaf2(nf), vaf3(nf);
  for (int f = 0; f < nf; ++f) {
    vaf1[f] = (l2_12[f] * c3[f] + l2_31[f] * c2[f]) / 8.0;
    vaf2[f] = (l2_23[f] * c1[f] + l2_12[f] * c3[f]) / 8.0;
    vaf3[f] = (l2_31[f] * c2[f] + l2_23[f] * c1[f]) / 8.0;
  }

  // obtuse handling (match your R logic)
  for (int f = 0; f < nf; ++f) {
    if (c1[f] < 0.0) { vaf1[f] = fa[f] / 2.0; vaf2[f] = fa[f] / 4.0; vaf3[f] = fa[f] / 4.0; }
    if (c2[f] < 0.0) { vaf1[f] = fa[f] / 4.0; vaf2[f] = fa[f] / 2.0; vaf3[f] = fa[f] / 4.0; }
    if (c3[f] < 0.0) { vaf1[f] = fa[f] / 4.0; vaf2[f] = fa[f] / 4.0; vaf3[f] = fa[f] / 2.0; }
  }

  // accumulate to vertices
  for (int f = 0; f < nf; ++f) {
    const int i1 = face(f, 0), i2 = face(f, 1), i3 = face(f, 2);
    va[i1 - 1] += vaf1[f];
    va[i2 - 1] += vaf2[f];
    va[i3 - 1] += vaf3[f];
  }

  return va;
}
