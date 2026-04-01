#include <Rcpp.h>
using namespace Rcpp;

// Ray-casting test: return true if (X,Y) inside polygon with N vertices
static inline bool in_polygon(int N, const double* xp, const double* yp, double X, double Y) {
  bool c = false;
  for (int i = 0, j = N - 1; i < N; j = i++) {
    const bool cond =
      (((yp[i] <= Y) && (Y < yp[j])) || ((yp[j] <= Y) && (Y < yp[i]))) &&
      (X < (xp[j] - xp[i]) * (Y - yp[i]) / (yp[j] - yp[i]) + xp[i]);
    if (cond) c = !c;
  }
  return c;
}

// [[Rcpp::export]]
LogicalVector isinpolygon(const NumericMatrix& polygon, const NumericMatrix& xy) {
  // polygon: m x 2 ; xy: n x 2
  if (polygon.ncol() != 2 || polygon.nrow() < 3)
    stop("polygon must be an m x 2 matrix with m >= 3");
  if (xy.ncol() != 2)
    stop("xy must be an n x 2 matrix");

  const int N = polygon.nrow();
  // Copy columns to contiguous buffers (faster & simple)
  std::vector<double> xp(N), yp(N);
  for (int i = 0; i < N; ++i) {
    xp[i] = polygon(i, 0);
    yp[i] = polygon(i, 1);
  }

  const int M = xy.nrow();
  LogicalVector out(M);
  for (int k = 0; k < M; ++k) {
    double X = xy(k, 0);
    double Y = xy(k, 1);
    if (Rcpp::NumericVector::is_na(X) || Rcpp::NumericVector::is_na(Y)) {
      out[k] = NA_LOGICAL;
    } else {
      out[k] = in_polygon(N, xp.data(), yp.data(), X, Y);
    }
  }
  return out;
}
