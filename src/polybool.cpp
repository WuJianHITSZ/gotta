//
//  polybool2_rcpp.cpp
//  mex2r
//
//  Created by Jian Wu on 2025/11/1.
//


// polybool_rcpp.cpp
// One-file Rcpp entry point that mirrors MATLAB's polybool.m and inlines polyboolmex logic.
// Exposes: polybool(pa, pb, op, ha = NULL, hb = NULL, ug = 1e6)
//
// Dependencies: Rcpp, Clipper (clipper.hpp + clipper.cpp in the same directory or your package src/)
//
// Usage in R (after sourceCpp):
//   pc_hc <- polybool(pa, pb, "and")
//   pc_hc$pc  # list of polygons (matrices n x 2)
//   pc_hc$hc  # logical vector: TRUE means "is a hole boundary"
//
// License of Clipper: Boost Software License (per Clipper project).

#include <Rcpp.h>
#include "clipper.hpp"

using namespace Rcpp;
using namespace ClipperLib;

// Helper: convert op string to Clipper ClipType
static ClipType parseOp(const std::string& op) {
  if (op.rfind("or", 0) == 0)        return ctUnion;
  if (op.rfind("and", 0) == 0)       return ctIntersection;
  if (op.rfind("notb", 0) == 0)      return ctDifference;
  if (op.rfind("diff", 0) == 0)      return ctDifference;
  if (op.rfind("xor", 0) == 0)       return ctXor;
  stop("polybool: unknown boolean set operation '%s'. Expected one of: 'or','and','notb'/'diff','xor'.", op);
}

// Ensure polygon orientation matches hole flag convention (like MATLAB MEX)
// Exterior: positive orientation; Hole: negative orientation.
static void enforceOrientation(Path& path, bool isHole) {
  bool orient = Orientation(path); // Clipper's Orientation
  if ((!isHole && !orient) || (isHole && orient)) {
    ReversePath(path);
  }
}

// Convert an R object (matrix or list-of-matrices) into Clipper Paths,
// applying scaling (ug) and enforcing orientation according to hole flags.
// - rPoly: either NumericMatrix (n x 2) or List of NumericMatrix
// - holeFlags: LogicalVector of same length as number of polygons (all FALSE if empty)
// - ug: user->grid scale (>0)
static Paths toPaths(SEXP rPoly, LogicalVector holeFlags, double ug) {
  if (!(ug > 0.0)) stop("polybool: 'ug' must be > 0.");
//  const double iug = 1.0 / ug;

  Paths out;
  if (Rf_isMatrix(rPoly)) {
    NumericMatrix M(rPoly);
    if (M.ncol() != 2) stop("polybool: polygon matrices must have 2 columns (x,y).");
    out.resize(1);
    out[0].resize(M.nrow());
    // scale to integer grid (round-half-up)
    for (int i = 0; i < M.nrow(); ++i) {
      cInt X = (cInt)std::llround(ug * M(i, 0));
      cInt Y = (cInt)std::llround(ug * M(i, 1));
      out[0][i] = IntPoint(X, Y);
    }
    bool isHole = (holeFlags.size() == 1) ? (holeFlags[0] == TRUE) : false;
    enforceOrientation(out[0], isHole);
  } else if (Rf_isNewList(rPoly)) {
    List L(rPoly);
    const int n = L.size();
    out.resize(n);
    // If no hole flags provided, default FALSE for all (exterior)
    LogicalVector hf = holeFlags;
    if (hf.size() == 0) {
      hf = LogicalVector(n, false);
    } else if (hf.size() != n) {
      stop("polybool: hole flag vector length must match number of polygons.");
    }

    for (int k = 0; k < n; ++k) {
      if (L[k] == R_NilValue) stop("polybool: empty polygon at index %d.", k + 1);
      NumericMatrix M = as<NumericMatrix>(L[k]);
      if (M.ncol() != 2) stop("polybool: polygon matrices must have 2 columns (x,y).");
      const int vnu = M.nrow();
      out[k].resize(vnu);
      for (int i = 0; i < vnu; ++i) {
        cInt X = (cInt)std::llround(ug * M(i, 0));
        cInt Y = (cInt)std::llround(ug * M(i, 1));
        out[k][i] = IntPoint(X, Y);
      }
      enforceOrientation(out[k], hf[k] == TRUE);
    }
  } else {
    stop("polybool: pa/pb must be a matrix (n x 2) or a list of such matrices.");
  }
  return out;
}

// Convert Clipper Paths back to R: list of (n x 2) numeric matrices, de-scaling by iug.
// Also compute hole flags: hc[k] = !Orientation(pc[k]) (same convention used in mex code).
static List fromPaths(const Paths& pc, double ug) {
  const double iug = 1.0 / ug;
  const size_t K = pc.size();
  List pc_out(K);
  LogicalVector hc(K);

  for (size_t k = 0; k < K; ++k) {
    const Path& P = pc[k];
    const int vnu = static_cast<int>(P.size());
    NumericMatrix M(vnu, 2);
    for (int i = 0; i < vnu; ++i) {
      M(i, 0) = iug * static_cast<double>(P[i].X);
      M(i, 1) = iug * static_cast<double>(P[i].Y);
    }
    pc_out[k] = M;
    // same hole flag convention as in the MEX: TRUE means "is a hole boundary"
    hc[k] = !Orientation(pc[k]);
  }
  return List::create(_["pc"] = pc_out, _["hc"] = hc);
}

// [[Rcpp::export]]
List polybool(SEXP pa,
              SEXP pb,
              std::string op,
              Rcpp::Nullable<LogicalVector> ha = R_NilValue,
              Rcpp::Nullable<LogicalVector> hb = R_NilValue,
              Rcpp::Nullable<NumericVector> ug_in = R_NilValue) {

  // Validate required args
  if (TYPEOF(pa) == NILSXP || TYPEOF(pb) == NILSXP) {
    stop("polybool: expecting at least pa, pb, and op.");
  }

  // Defaults matching MATLAB polybool.m
  double ug = 1e6;
  if (ug_in.isNotNull()) {
    NumericVector tmp(ug_in.get());
    if (tmp.size() != 1) stop("polybool: 'ug' must be a single numeric value.");
    ug = tmp[0];
  }

  // Prepare hole flags. If missing, we fill after we know counts inside toPaths.
  LogicalVector ha_vec, hb_vec;
  if (ha.isNotNull()) ha_vec = LogicalVector(ha.get());
  if (hb.isNotNull()) hb_vec = LogicalVector(hb.get());

  // Parse operation
  ClipType clip_op = parseOp(op);

  // Convert inputs to Clipper Paths with scaling and orientation
  Paths A = toPaths(pa, ha_vec, ug);
  Paths B = toPaths(pb, hb_vec, ug);

  // Clip
  Paths PC;
  Clipper C;
  C.AddPaths(A, ptSubject, true);
  C.AddPaths(B, ptClip,    true);

  // Even-Odd fill rule to mirror the MEX function
  bool ok = C.Execute(clip_op, PC, pftEvenOdd, pftEvenOdd);
  if (!ok) stop("polybool: Clipper library failed to execute operation.");

  // Return list(pc = list of polygons, hc = hole flags)
  return fromPaths(PC, ug);
}
