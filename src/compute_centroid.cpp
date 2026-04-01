// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#define CLIPPER_LIB_NAMESPACE clipperlib_embed
#include "clipper.hpp"

#include "polybool.h"

using namespace Rcpp;
using namespace ClipperLib;

// =====================================================
// === polybool(): polygon Boolean operation (embedded)
// =====================================================
//static ClipType parseOp(const std::string& op) {
//  if (op.rfind("or",   0)==0) return ctUnion;
//  if (op.rfind("and",  0)==0) return ctIntersection;
//  if (op.rfind("notb", 0)==0) return ctDifference;
//  if (op.rfind("diff", 0)==0) return ctDifference;
//  if (op.rfind("xor",  0)==0) return ctXor;
//  stop("polybool: unknown op. Use 'or','and','notb'/'diff','xor'.");
//}
//
//static void enforceOrientation(Path& path, bool isHole) {
//  bool orient = Orientation(path);
//  if ((!isHole && !orient) || (isHole && orient)) ReversePath(path);
//}
//
//static Paths toPaths(SEXP rPoly, LogicalVector holeFlags, double ug) {
//  if (!(ug > 0.0)) stop("polybool: 'ug' must be > 0.");
//  Paths out;
//  if (Rf_isMatrix(rPoly)) {
//    NumericMatrix M = as<NumericMatrix>(rPoly);
//    if (M.ncol()!=2) stop("polybool: polygon matrices must be n x 2.");
//    out.resize(1); out[0].resize(M.nrow());
//    for (int i=0;i<M.nrow();++i)
//      out[0][i] = IntPoint((cInt)llround(ug*M(i,0)), (cInt)llround(ug*M(i,1)));
//    bool isHole = (holeFlags.size()==1) ? (holeFlags[0]==TRUE) : false;
//    enforceOrientation(out[0], isHole);
//  } else if (Rf_isNewList(rPoly)) {
//    List L = as<List>(rPoly);
//    int n = L.size(); out.resize(n);
//    LogicalVector hf = holeFlags.size()==0 ? LogicalVector(n, false) : holeFlags;
//    if (hf.size()!=n) stop("polybool: holeFlags length mismatch.");
//    for (int k=0;k<n;++k) {
//      if (L[k]==R_NilValue) stop("polybool: empty polygon at index %d.", k+1);
//      NumericMatrix M = as<NumericMatrix>(L[k]);
//      if (M.ncol()!=2) stop("polybool: polygon matrices must be n x 2.");
//      out[k].resize(M.nrow());
//      for (int i=0;i<M.nrow();++i)
//        out[k][i] = IntPoint((cInt)llround(ug*M(i,0)), (cInt)llround(ug*M(i,1)));
//      enforceOrientation(out[k], hf[k]==TRUE);
//    }
//  } else {
//    stop("polybool: pa/pb must be a matrix (n x 2) or a list of matrices.");
//  }
//  return out;
//}
//
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
//
//// [[Rcpp::export]]
//List polybool(SEXP pa, SEXP pb, std::string op,
//              Rcpp::Nullable<LogicalVector> ha = R_NilValue,
//              Rcpp::Nullable<LogicalVector> hb = R_NilValue,
//              Rcpp::Nullable<NumericVector> ug_in = R_NilValue) {
//  if (TYPEOF(pa)==NILSXP || TYPEOF(pb)==NILSXP)
//    stop("polybool: need pa, pb, op.");
//  double ug = 1e6;
//  if (ug_in.isNotNull()) {
//    NumericVector tmp(ug_in.get());
//    if (tmp.size()!=1) stop("'ug' must be length-1.");
//    ug = tmp[0];
//  }
//  LogicalVector ha_vec, hb_vec;
//  if (ha.isNotNull()) ha_vec = LogicalVector(ha.get());
//  if (hb.isNotNull()) hb_vec = LogicalVector(hb.get());
//  ClipType clip_op = parseOp(op);
//  Paths A = toPaths(pa, ha_vec, ug), B = toPaths(pb, hb_vec, ug), PC;
//  Clipper C; C.AddPaths(A, ptSubject, true); C.AddPaths(B, ptClip, true);
//  if (!C.Execute(clip_op, PC, pftEvenOdd, pftEvenOdd))
//    stop("polybool: Clipper failed.");
//  return fromPaths(PC, ug);
//}

// =====================================================
// === compute_centroid(): main exported function
// =====================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_centroid(const NumericMatrix& cp, const List& pd) {
  if (!pd.containsElementNamed("cell") || !pd.containsElementNamed("dpe"))
    stop("pd must contain fields 'cell' (list) and 'dpe' (matrix).");

  List cell = pd["cell"];
  NumericMatrix dpe = pd["dpe"];
  const int nc = cell.size();

  if (dpe.ncol() != 2) stop("pd$dpe must be an n x 2 matrix.");
  if (cp.ncol()  != 2) stop("cp must be an m x 2 matrix.");

  NumericMatrix cc(nc, 2);

  for (int i = 0; i < nc; ++i) {
    IntegerVector ci = cell[i];
    const int m = ci.size();
    if (m == 0) { cc(i,0)=NA_REAL; cc(i,1)=NA_REAL; continue; }

    NumericMatrix di(m, 2);
    for (int k = 0; k < m; ++k) {
      int r = ci[k] - 1;
      if (r < 0 || r >= dpe.nrow())
        stop("cell index out of range (cell %d)", i+1);
      di(k,0) = dpe(r,0);
      di(k,1) = dpe(r,1);
    }

    if (m >= 2 && ci[0] == ci[m-1]) {
      NumericMatrix di2(m-1, 2);
      for (int r=0;r<m-1;++r){ di2(r,0)=di(r,0); di2(r,1)=di(r,1); }
      di = di2;
    }

    List pb = polybool(cp, di, std::string("and"));

    double cx = NA_REAL, cy = NA_REAL;
    if (pb.containsElementNamed("pc")) {
      List pc = pb["pc"];
      if (pc.size()>0 && Rf_isMatrix(pc[0])) {
        NumericMatrix M = as<NumericMatrix>(pc[0]);  // <-- fixed here
        double sx=0.0, sy=0.0;
        for (int r=0;r<M.nrow();++r){ sx+=M(r,0); sy+=M(r,1); }
        cx = sx / (double)M.nrow();
        cy = sy / (double)M.nrow();
        cc(i,0)=cx; cc(i,1)=cy;
        continue;
      }
    }

    double sx=0.0, sy=0.0;
    for (int r=0;r<di.nrow();++r){ sx+=di(r,0); sy+=di(r,1); }
    cc(i,0)=sx/(double)di.nrow();
    cc(i,1)=sy/(double)di.nrow();
  }

  return cc;
}
