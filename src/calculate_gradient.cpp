//
//  calculate_gradient_rcpp.cpp
//  mex2r
//
//  Created by Jian Wu on 2025/11/2.
//


// calculate_gradient_rcpp.cpp
// One-file build: merges isinpolygon, polyarea, polybool, and adds calculate_gradient()
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "polybool.h"
#include "polyarea.h"
#include "isinpolygon.h"

using namespace Rcpp;

/* =========================
   isinpolygon  (from upload)
   ========================= */
//static inline bool in_polygon_rc(int N, const double* xp, const double* yp, double X, double Y) {
//  bool c = false;
//  for (int i = 0, j = N - 1; i < N; j = i++) {
//    const bool cond =
//      (((yp[i] <= Y) && (Y < yp[j])) || ((yp[j] <= Y) && (Y < yp[i]))) &&
//      (X < (xp[j] - xp[i]) * (Y - yp[i]) / (yp[j] - yp[i]) + xp[i]);
//    if (cond) c = !c;
//  }
//  return c;
//}
//
//// [[Rcpp::export]]
//LogicalVector isinpolygon(const NumericMatrix& polygon, const NumericMatrix& xy) {
//  if (polygon.ncol() != 2 || polygon.nrow() < 3)
//    stop("polygon must be an m x 2 matrix with m >= 3");
//  if (xy.ncol() != 2)
//    stop("xy must be an n x 2 matrix");
//
//  const int N = polygon.nrow();
//  std::vector<double> xp(N), yp(N);
//  for (int i = 0; i < N; ++i) { xp[i] = polygon(i,0); yp[i] = polygon(i,1); }
//
//  const int M = xy.nrow();
//  LogicalVector out(M);
//  for (int k = 0; k < M; ++k) {
//    double X = xy(k,0), Y = xy(k,1);
//    if (Rcpp::NumericVector::is_na(X) || Rcpp::NumericVector::is_na(Y)) out[k] = NA_LOGICAL;
//    else out[k] = in_polygon_rc(N, xp.data(), yp.data(), X, Y);
//  }
//  return out;
//}

/* =========================
   polyarea  (from upload)
   ========================= */
//static inline double polygon_area_matrix(const NumericMatrix& poly, int m) {
//  double A = 0.0;
//  for (int i = 0; i < m - 1; ++i)
//    A += poly(i,0) * poly(i+1,1) - poly(i+1,0) * poly(i,1);
//  A += poly(m-1,0) * poly(0,1) - poly(0,0) * poly(m-1,1);
//  A *= 0.5;
//  return A >= 0 ? A : -A;
//}
//static inline bool has_na_first_m(const NumericMatrix& poly, int m) {
//  for (int i = 0; i < m; ++i)
//    if (NumericMatrix::is_na(poly(i,0)) || NumericMatrix::is_na(poly(i,1))) return true;
//  return false;
//}
//
//// hej [[Rcpp::export]]
//NumericVector polyarea(List pa) {
//  const R_xlen_t n = pa.size();
//  NumericVector out(n);
//
//  SEXP dims = pa.attr("dim");
//  if (!Rf_isNull(dims)) { out.attr("dim") = dims; out.attr("dimnames") = pa.attr("dimnames"); }
//
//  for (R_xlen_t k = 0; k < n; ++k) {
//    if (!Rf_isMatrix(pa[k])) stop("Each element of 'pa' must be an n x 2 numeric matrix.");
//    NumericMatrix poly(pa[k]);
//    if (poly.ncol()!=2) stop("Each polygon must be an n x 2 matrix.");
//    int M = poly.nrow();
//    if (M == 0) stop("Empty polygon matrix at index %lld.", (long long)k+1);
//    int m = M;
//    if (M >= 2 && poly(0,0)==poly(M-1,0) && poly(0,1)==poly(M-1,1)) m = M-1;
//    if (m <= 2) stop("Polygon must have at least 3 vertices (index %lld).", (long long)k+1);
//    if (has_na_first_m(poly, m)) { out[k] = NA_REAL; continue; }
//    out[k] = polygon_area_matrix(poly, m);
//  }
//  return out;
//}
//// hej [[Rcpp::export]]
//double polyarea_matrix(const NumericMatrix& poly) {
//  if (poly.ncol()!=2) stop("Polygon must be an n x 2 matrix.");
//  int M = poly.nrow(), m = M;
//  if (M >= 2 && poly(0,0)==poly(M-1,0) && poly(0,1)==poly(M-1,1)) m = M-1;
//  if (m <= 2) stop("Polygon must have at least 3 vertices.");
//  if (has_na_first_m(poly, m)) return NA_REAL;
//  return polygon_area_matrix(poly, m);
//}

/* =========================
   polybool  (from upload, header-only shim for single file)
   ========================= */
// Clipper (bundled) — single-header include. If you keep Clipper as separate files,
// just place "clipper.hpp/clipper.cpp" in src/ and replace this small embed with #include "clipper.hpp".
// For portability here, we include the header only and rely on the header's implementation.
// If your local Clipper needs the .cpp, put it alongside and add #include "clipper.hpp" instead.
//#define CLIPPER_LIB_NAMESPACE clipperlib_embed
//#include "clipper.hpp"    // make sure this header is available in your src/ (or vendor it)
//
//using namespace ClipperLib;
//
//static ClipType parseOp(const std::string& op) {
//  if (op.rfind("or",   0)==0) return ctUnion;
//  if (op.rfind("and",  0)==0) return ctIntersection;
//  if (op.rfind("notb", 0)==0) return ctDifference;
//  if (op.rfind("diff", 0)==0) return ctDifference;
//  if (op.rfind("xor",  0)==0) return ctXor;
//  stop("polybool: unknown op. Use 'or','and','notb'/'diff','xor'.");
//}
//static void enforceOrientation(Path& path, bool isHole) {
//  bool orient = Orientation(path);
//  if ((!isHole && !orient) || (isHole && orient)) ReversePath(path);
//}
//static Paths toPaths(SEXP rPoly, LogicalVector holeFlags, double ug) {
//  if (!(ug > 0.0)) stop("polybool: 'ug' must be > 0.");
//  Paths out;
//  if (Rf_isMatrix(rPoly)) {
//    NumericMatrix M(rPoly);
//    if (M.ncol()!=2) stop("polybool: polygon matrices must be n x 2.");
//    out.resize(1); out[0].resize(M.nrow());
//    for (int i=0;i<M.nrow();++i) out[0][i] = IntPoint((cInt)llround(ug*M(i,0)), (cInt)llround(ug*M(i,1)));
//    bool isHole = (holeFlags.size()==1) ? (holeFlags[0]==TRUE) : false;
//    enforceOrientation(out[0], isHole);
//  } else if (Rf_isNewList(rPoly)) {
//    List L(rPoly); int n = L.size(); out.resize(n);
//    LogicalVector hf = holeFlags.size()==0 ? LogicalVector(n, false) : holeFlags;
//    if (hf.size()!=n) stop("polybool: holeFlags length mismatch.");
//    for (int k=0;k<n;++k) {
//      if (L[k]==R_NilValue) stop("polybool: empty polygon at index %d.", k+1);
//      NumericMatrix M(L[k]);
//      if (M.ncol()!=2) stop("polybool: polygon matrices must be n x 2.");
//      out[k].resize(M.nrow());
//      for (int i=0;i<M.nrow();++i) out[k][i] = IntPoint((cInt)llround(ug*M(i,0)), (cInt)llround(ug*M(i,1)));
//      enforceOrientation(out[k], hf[k]==TRUE);
//    }
//  } else {
//    stop("polybool: pa/pb must be a matrix (n x 2) or a list of matrices.");
//  }
//  return out;
//}
//static List fromPaths(const Paths& pc, double ug) {
//  const double iug = 1.0/ug;
//  const size_t K = pc.size();
//  List pc_out(K); LogicalVector hc(K);
//  for (size_t k=0;k<K;++k) {
//    const Path& P = pc[k]; const int vnu = (int)P.size();
//    NumericMatrix M(vnu,2);
//    for (int i=0;i<vnu;++i) { M(i,0)=iug*(double)P[i].X; M(i,1)=iug*(double)P[i].Y; }
//    pc_out[k]=M; hc[k] = !Orientation(pc[k]);
//  }
//  return List::create(_["pc"]=pc_out, _["hc"]=hc);
//}
//// hej [[Rcpp::export]]
//List polybool(SEXP pa, SEXP pb, std::string op,
//              Rcpp::Nullable<LogicalVector> ha = R_NilValue,
//              Rcpp::Nullable<LogicalVector> hb = R_NilValue,
//              Rcpp::Nullable<NumericVector> ug_in = R_NilValue) {
//  if (TYPEOF(pa)==NILSXP || TYPEOF(pb)==NILSXP) stop("polybool: need pa, pb, op.");
//  double ug = 1e6;
//  if (ug_in.isNotNull()) { NumericVector tmp(ug_in.get()); if (tmp.size()!=1) stop("'ug' must be length-1."); ug = tmp[0]; }
//  LogicalVector ha_vec, hb_vec; if (ha.isNotNull()) ha_vec = LogicalVector(ha.get()); if (hb.isNotNull()) hb_vec = LogicalVector(hb.get());
//  ClipType clip_op = parseOp(op);
//  Paths A = toPaths(pa, ha_vec, ug), B = toPaths(pb, hb_vec, ug), PC;
//  Clipper C; C.AddPaths(A, ptSubject, true); C.AddPaths(B, ptClip, true);
//  if (!C.Execute(clip_op, PC, pftEvenOdd, pftEvenOdd)) stop("polybool: Clipper failed.");
//  return fromPaths(PC, ug);
//}

/* =========================
   small helpers for calculate_gradient
   ========================= */
// area of triangle given 3x2 (or 3x3 with z=0) matrix of vertices
static inline double tri_area2(const NumericMatrix& tri) {
  // Shoelace in 2D
  double x1 = tri(0,0), y1 = tri(0,1);
  double x2 = tri(1,0), y2 = tri(1,1);
  double x3 = tri(2,0), y3 = tri(2,1);
  return 0.5 * std::fabs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
}

// Return mean of a NumericVector (naive)
static inline double mean_numeric(const NumericVector& v) {
  if (v.size()==0) return NA_REAL;
  double s=0.0; for (int i=0;i<v.size();++i) s += v[i];
  return s / (double)v.size();
}

// Call R function sigma(.) and return average of its numeric result
static inline double sigma_mean(Function sigma, const SEXP arg) {
  SEXP resS = sigma(arg);
  NumericVector res(resS);
  return mean_numeric(res);
}

// Call sigma on a 1x2 point (as numeric(2)) and return a scalar (mean if vector)
static inline double sigma_on_point(Function sigma, double x, double y) {
  NumericVector pt = NumericVector::create(x, y);
  return sigma_mean(sigma, pt);
}

/* =========================
   calculate_gradient (main)
   ========================= */
// [[Rcpp::export]]
NumericVector calculate_gradient(NumericMatrix cp,
                                     List pd,
                                     Function sigma,
                                     bool more_accurate = false) {
  // validate pd
  if (!pd.containsElementNamed("dpe") || !pd.containsElementNamed("cell"))
    stop("pd must contain fields 'dpe' (matrix) and 'cell' (list).");
  NumericMatrix dpe = pd["dpe"];
  List cell = pd["cell"];
  const int nc = cell.size();

  if (cp.ncol()!=2) stop("cp must be an n x 2 matrix.");
  if (dpe.ncol()!=2) stop("pd$dpe must be an n x 2 matrix.");

  // in2: which dpe points are inside cp
  LogicalVector in2 = isinpolygon(cp, dpe);
  std::vector<bool> in1(nc, true);

  // mark cells fully inside
  for (int i=0;i<nc;++i) {
    IntegerVector idx = cell[i];
    for (int k=0;k<idx.size();++k) {
      int r = idx[k] - 1; // 1-based -> 0-based
      if (r < 0 || r >= dpe.nrow()) continue;
      if (in2[r] != TRUE) { in1[i] = false; break; }
    }
  }

  // if cp is explicitly closed with last==first, drop the last
  if (cp.nrow() >= 2) {
    if (cp(0,0)==cp(cp.nrow()-1,0) && cp(0,1)==cp(cp.nrow()-1,1)) {
      cp = cp(Range(0, cp.nrow()-2), _);
    }
  }

  NumericVector D(nc, 0.0);

  for (int i=0;i<nc;++i) {
    IntegerVector idx = cell[i];
    const int m = idx.size();
    if (m == 0) { D[i]=0.0; continue; }

    // Build ci (cell polygon in dpe order)
    NumericMatrix ci(m, 2);
    for (int k=0;k<m;++k) {
      int r = idx[k]-1; if (r<0 || r>=dpe.nrow()) stop("cell index out of range.");
      ci(k,0) = dpe(r,0); ci(k,1) = dpe(r,1);
    }

    if (more_accurate) {
      // If not fully inside, intersect first
      if (!in1[i]) {
        List res = polybool(cp, ci, "and");     // returns list with elements $pc, $hc
        List pc = res["pc"];
        if (pc.size()==0) { D[i]=0.0; continue; }
        NumericMatrix ci2 = pc[0];
        // ensure closed
        if (!(ci2(0,0)==ci2(ci2.nrow()-1,0) && ci2(0,1)==ci2(ci2.nrow()-1,1))) {
          NumericMatrix tmp(ci2.nrow()+1, 2);
          for (int r=0;r<ci2.nrow();++r){ tmp(r,0)=ci2(r,0); tmp(r,1)=ci2(r,1); }
          tmp(ci2.nrow(),0)=ci2(0,0); tmp(ci2.nrow(),1)=ci2(0,1);
          ci = tmp;
        } else {
          ci = ci2;
        }
      }

      if (ci.nrow() < 4) { D[i]=0.0; continue; } // need at least 3 distinct + repeat of first

      // centroid of rim (excluding closing vertex)
      NumericMatrix rim = ci(Range(0, ci.nrow()-2), _);
      double cx = 0.0, cy = 0.0;
      for (int r=0;r<rim.nrow();++r){ cx += rim(r,0); cy += rim(r,1); }
      cx /= (double)rim.nrow(); cy /= (double)rim.nrow();

      double acc = 0.0;
      for (int j=0;j<ci.nrow()-1;++j) {
        NumericMatrix tri(3,2);
        tri(0,0)=cx;        tri(0,1)=cy;
        tri(1,0)=ci(j,0);   tri(1,1)=ci(j,1);
        tri(2,0)=ci(j+1,0); tri(2,1)=ci(j+1,1);

        const double a_tri = tri_area2(tri);

        // sigma at triangle mean
        double mx = (tri(0,0)+tri(1,0)+tri(2,0))/3.0;
        double my = (tri(0,1)+tri(1,1)+tri(2,1))/3.0;
        const double mu_tri = sigma_on_point(sigma, mx, my);

        acc += a_tri * mu_tri;
      }
      D[i] = acc;
      continue;
    }

    // Approximate path
    if (in1[i]) {
      // remove closing point
      NumericMatrix ci_core(ci.nrow()-1, 2);
      for (int r=0;r<ci_core.nrow();++r){ ci_core(r,0)=ci(r,0); ci_core(r,1)=ci(r,1); }

      // (2*mean(sigma(ci_core)) + sigma(colMeans(ci_core)))/3
      double sig_mean_boundary = sigma_mean(sigma, ci_core);
      double cx = 0.0, cy = 0.0;
      for (int r=0;r<ci_core.nrow();++r){ cx+=ci_core(r,0); cy+=ci_core(r,1); }
      cx /= (double)ci_core.nrow(); cy /= (double)ci_core.nrow();
      double sig_at_centroid = sigma_on_point(sigma, cx, cy);
      double mui = (2.0*sig_mean_boundary + sig_at_centroid)/3.0;

      double area = polyarea(List::create(ci_core))[0];
      D[i] = area * mui;
    } else {
      // intersect then integrate on intersection
      List res = polybool(cp, ci, "and");
      List pc = res["pc"];
      if (pc.size()==0) { D[i]=0.0; continue; }
      NumericMatrix xy = pc[0];
      if (xy.nrow() < 3) { D[i]=0.0; continue; }

      double sig_mean_boundary = sigma_mean(sigma, xy);
      double cx = 0.0, cy = 0.0;
      for (int r=0;r<xy.nrow();++r){ cx+=xy(r,0); cy+=xy(r,1); }
      cx /= (double)xy.nrow(); cy /= (double)xy.nrow();
      double sig_at_centroid = sigma_on_point(sigma, cx, cy);
      double mui = (2.0*sig_mean_boundary + sig_at_centroid)/3.0;

      double area = polyarea(List::create(xy))[0];
      D[i] = area * mui;
    }
  }

  return D;
}
