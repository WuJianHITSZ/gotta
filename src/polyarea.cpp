// polyarea_rcpp.cpp
#include <Rcpp.h>
using namespace Rcpp;

// ---- helper: polygon area (shoelace) on first m rows of matrix ----
static inline double polygon_area_matrix(const NumericMatrix& poly, int m) {
  // Assumes poly has 2 columns and m >= 3
  double A = 0.0;
  // Use standard shoelace with wrap-around
  for (int i = 0; i < m - 1; ++i) {
    A += poly(i, 0) * poly(i + 1, 1) - poly(i + 1, 0) * poly(i, 1);
  }
  // close the polygon from last to first
  A += poly(m - 1, 0) * poly(0, 1) - poly(0, 0) * poly(m - 1, 1);
  A *= 0.5;
  return A >= 0 ? A : -A;  // absolute value
}

// ---- helper: check for any NA in the first m rows of a 2-col matrix ----
static inline bool has_na_first_m(const NumericMatrix& poly, int m) {
  for (int i = 0; i < m; ++i) {
    if (NumericMatrix::is_na(poly(i, 0)) || NumericMatrix::is_na(poly(i, 1))) {
      return true;
    }
  }
  return false;
}

// [[Rcpp::export]]
NumericVector polyarea(List pa) {
  const R_xlen_t n = pa.size();
  NumericVector out(n);

  // Preserve shape if input is a list-array (like MATLAB cell array dims)
  SEXP dims = pa.attr("dim");
  if (!Rf_isNull(dims)) {
    out.attr("dim") = dims;            // same shape as input list
    out.attr("dimnames") = pa.attr("dimnames");
  }

  for (R_xlen_t k = 0; k < n; ++k) {
    SEXP el = pa[k];
    if (!Rf_isMatrix(el))
      stop("Each element of 'pa' must be a numeric matrix (n x 2).");

    NumericMatrix poly(el);
    if (poly.ncol() != 2)
      stop("Each polygon must be an n x 2 matrix.");

    int M = poly.nrow();
    if (M == 0)
      stop("Empty polygon matrix at index %lld.", static_cast<long long>(k + 1));

    // Duplicate-last-vertex check (exact equality, matches the MEX)
    int m = M;
    if (M >= 2 &&
        poly(0, 0) == poly(M - 1, 0) &&
        poly(0, 1) == poly(M - 1, 1)) {
      m = M - 1;
    }

    if (m <= 2)
      stop("Polygons must have at least 3 vertices (index %lld).",
           static_cast<long long>(k + 1));

    if (has_na_first_m(poly, m)) {
      out[k] = NA_REAL;                // propagate NA if any vertex is NA
      continue;
    }

    out[k] = polygon_area_matrix(poly, m);
  }

  return out;
}

// Optional convenience: allow a single matrix input too
// [[Rcpp::export]]
double polyarea_matrix(const NumericMatrix& poly) {
  if (poly.ncol() != 2)
    stop("Polygon must be an n x 2 matrix.");
  int M = poly.nrow();
  int m = M;
  if (M >= 2 &&
      poly(0, 0) == poly(M - 1, 0) &&
      poly(0, 1) == poly(M - 1, 1)) {
    m = M - 1;
  }
  if (m <= 2) stop("Polygon must have at least 3 vertices.");
  if (has_na_first_m(poly, m)) return NA_REAL;
  return polygon_area_matrix(poly, m);
}
