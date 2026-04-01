// calculate_hessian_rcpp.cpp
// One-file build: isinpolygon + intersectEdgePolygon + calculate_hessian
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <unordered_map>

#include "isinpolygon.h"

using namespace Rcpp;

/* ============================================================
 * isinpolygon (merged)
 * ============================================================ */
//
//// Ray-casting test: return true if (X,Y) inside polygon with N vertices
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
//    if (NumericVector::is_na(X) || NumericVector::is_na(Y)) out[k] = NA_LOGICAL;
//    else out[k] = in_polygon_rc(N, xp.data(), yp.data(), X, Y);
//  }
//  return out;
//}

/* ============================================================
 * intersectEdgePolygon (fast C++ port)
 * ============================================================ */

// 2D helpers
static inline double cross2(const double ax, const double ay, const double bx, const double by) {
  return ax*by - ay*bx;
}
static inline double dot2(const double ax, const double ay, const double bx, const double by) {
  return ax*bx + ay*by;
}
static inline bool on01(double x, double tol) { return x >= -tol && x <= 1.0 + tol; }

// Solve segment-segment intersections for a single polygon edge q..q+s
static inline void seg_seg_intersections(const double* p, const double* r,
                                         const double* q, const double* s,
                                         double tol,
                                         std::vector< std::array<double,2> >& outPts) {
  const double rxs   = cross2(r[0], r[1], s[0], s[1]);
  const double qmpx  = q[0] - p[0];
  const double qmpy  = q[1] - p[1];
  const double qmpxr = cross2(qmpx, qmpy, r[0], r[1]);

  if (std::fabs(rxs) > tol) {
    const double t = cross2(qmpx, qmpy, s[0], s[1]) / rxs;
    const double u = qmpxr / rxs;
    if (on01(t, tol) && on01(u, tol)) {
      outPts.push_back({ p[0] + t*r[0], p[1] + t*r[1] });
    }
    return;
  }

  // Parallel case
  if (std::fabs(qmpxr) > tol) return; // parallel & non-colinear

  // Colinear: project q and q+s onto r and clip overlap with [0,1]
  const double r2 = dot2(r[0], r[1], r[0], r[1]);
  if (r2 < tol) {
    // r is degenerate (edge is a point)
    const double s2 = dot2(s[0], s[1], s[0], s[1]);
    if (s2 < tol) {
      const double dx = p[0]-q[0], dy = p[1]-q[1];
      if (dx*dx + dy*dy <= tol*tol) outPts.push_back({p[0], p[1]});
      return;
    } else {
      const double u = dot2(p[0]-q[0], p[1]-q[1], s[0], s[1]) / s2;
      if (on01(u, tol)) outPts.push_back({p[0], p[1]});
      return;
    }
  }

  const double t0 = dot2(q[0]-p[0], q[1]-p[1], r[0], r[1]) / r2;
  const double t1 = dot2(q[0]+s[0]-p[0], q[1]+s[1]-p[1], r[0], r[1]) / r2;
  const double tmin = std::min(t0, t1);
  const double tmax = std::max(t0, t1);

  const double a = std::max(0.0, tmin);
  const double b = std::min(1.0, tmax);
  if (b + tol < a) return; // no overlap

  const double Ax = p[0] + a*r[0], Ay = p[1] + a*r[1];
  const double Bx = p[0] + b*r[0], By = p[1] + b*r[1];
  const double dx = Ax - Bx, dy = Ay - By;
  if (dx*dx + dy*dy <= tol*tol) outPts.push_back({Ax, Ay});
  else { outPts.push_back({Ax, Ay}); outPts.push_back({Bx, By}); }
}

// [[Rcpp::export]]
NumericMatrix intersectEdgePolygon(const NumericVector& edge, const NumericMatrix& poly, double tol = 1e-9) {
  if (edge.size() != 4) stop("edge must be numeric length-4 c(x1,y1,x2,y2)");
  if (poly.ncol() != 2) stop("poly must be an N x 2 matrix");

  double p[2] = { edge[0], edge[1] };
  double p2[2] = { edge[2], edge[3] };
  double r[2] = { p2[0]-p[0], p2[1]-p[1] };

  const int n = poly.nrow();
  if (n < 2) {
    NumericMatrix out(0,2); colnames(out) = CharacterVector::create("x","y"); return out;
  }

  std::vector< std::array<double,2> > pts;
  pts.reserve(n);

  for (int i=0; i<n; ++i) {
    const int j = (i+1 < n) ? (i+1) : 0;
    double q[2]  = { poly(i,0), poly(i,1) };
    double q2[2] = { poly(j,0), poly(j,1) };
    double s[2]  = { q2[0]-q[0], q2[1]-q[1] };
    seg_seg_intersections(p, r, q, s, tol, pts);
  }

  if (pts.empty()) {
    NumericMatrix out(0,2); colnames(out) = CharacterVector::create("x","y"); return out;
  }

  // compute t along edge to sort and deduplicate by tolerance grid
  std::vector<double> tvals; tvals.reserve(pts.size());
  const double r2 = r[0]*r[0] + r[1]*r[1];
  for (auto &pt : pts) {
    double t = (r2 > tol) ? ((pt[0]-p[0])*r[0] + (pt[1]-p[1])*r[1]) / r2 : 0.0;
    tvals.push_back(t);
  }

  // dedup (grid)
  std::vector<size_t> keep; keep.reserve(pts.size());
  std::unordered_map<long long, bool> seen;
  const double inv = 1.0 / tol;
  for (size_t k=0;k<pts.size();++k) {
    long long kx = llround(pts[k][0] * inv);
    long long ky = llround(pts[k][1] * inv);
    long long key = (kx<<32) ^ (ky & 0xffffffffLL);
    if (!seen.count(key)) { seen[key]=true; keep.push_back(k); }
  }

  // sort by t
  std::vector<size_t> ord(keep.begin(), keep.end());
  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b){ return tvals[a] < tvals[b]; });

  NumericMatrix M(ord.size(), 2);
  for (size_t k=0;k<ord.size();++k){ M(k,0)=pts[ord[k]][0]; M(k,1)=pts[ord[k]][1]; }
  colnames(M) = CharacterVector::create("x","y");
  return M;
}

/* ============================================================
 * calculate_hessian (MAIN)
 * ============================================================ */

static inline double vnorm2_2d(double x, double y){ return std::sqrt(x*x + y*y); }

static inline double sigma_at_point(Function sigma, double x, double y) {
  NumericMatrix pt(1,2); pt(0,0)=x; pt(0,1)=y;
  NumericVector res = sigma(pt);
  if (res.size()==0) return NA_REAL;
  double s=0.0; for (int i=0;i<res.size();++i) s+=res[i];
  return s / (double)res.size();
}

// [[Rcpp::export]]
SEXP calculate_hessian(NumericMatrix cp, List pd, Function sigma) {
  // Validate pd
  if (!pd.containsElementNamed("cell") || !pd.containsElementNamed("dpe") || !pd.containsElementNamed("uv"))
    stop("pd must contain fields 'cell', 'dpe', and 'uv'.");
  NumericMatrix dpe = pd["dpe"];
  List          cell = pd["cell"];
  NumericMatrix uv   = pd["uv"];

  if (cp.ncol()!=2) stop("cp must be an n x 2 polygon.");
  if (dpe.ncol()!=2) stop("pd$dpe must be an m x 2 matrix.");
  if (uv.ncol()!=2)  stop("pd$uv must be an nc x 2 matrix.");

  const int nc = cell.size();
  if (nc == 0) {
    Function sparseMatrix("sparseMatrix");
    return sparseMatrix(_["i"]=IntegerVector(0), _["j"]=IntegerVector(0),
                        _["x"]=NumericVector(0), _["dims"]=IntegerVector::create(0,0));
  }

  // Count oriented edges: ne = sum(length(ci)) - nc
  long long ne = 0;
  for (int i=0;i<nc;++i) { IntegerVector ci = cell[i]; ne += std::max(0, (int)ci.size()-1); }
  if (ne <= 0) {
    Function sparseMatrix("sparseMatrix");
    return sparseMatrix(_["i"]=IntegerVector(nc), _["j"]=IntegerVector(nc),
                        _["x"]=NumericVector(nc,0.0), _["dims"]=IntegerVector::create(nc,nc));
  }

  // Build oriented edge triplets (I, J, owner)
  std::vector<int> I(ne), J(ne), OWN(ne);
  long long k = 0;
  for (int i=0;i<nc;++i) {
    IntegerVector ci = cell[i];
    if (ci.size() < 2) continue;
    for (int j=0;j<ci.size()-1;++j) {
      I[k]   = ci[j];
      J[k]   = ci[j+1];
      OWN[k] = i+1; // 1-based owner cell id
      ++k;
    }
  }
  const long long E = k; // actual # edges collected

  // Map each directed edge (I->J) to its owner cell
  auto key2 = [](int a, int b)->long long { return ( (long long)a << 32 ) ^ (long long)b; };
  std::unordered_map<long long, int> owner_of; owner_of.reserve((size_t)E*2);
  for (long long e = 0; e < E; ++e) {
    owner_of.emplace(key2(I[e], J[e]), OWN[e]);
  }

  // in_poly for each dpe point
  LogicalVector in_poly = isinpolygon(cp, dpe);

  // p = sigma(dpe); accept scalar or vector
  NumericVector p_sigma = sigma(dpe);
  NumericVector p(dpe.nrow());
  if (p_sigma.size()==1) {
    std::fill(p.begin(), p.end(), p_sigma[0]);
  } else {
    if (p_sigma.size() != dpe.nrow()) stop("sigma(dpe) must return scalar or vector of length nrow(dpe).");
    for (int r=0;r<p.size();++r) p[r] = p_sigma[r];
  }

  // Prepare output triplets for H (off-diagonals first)
  std::vector<int>    I2; I2.reserve((size_t)E + (size_t)nc);
  std::vector<int>    J2; J2.reserve((size_t)E + (size_t)nc);
  std::vector<double> X2; X2.reserve((size_t)E + (size_t)nc);
  std::vector<double> rowSum(nc, 0.0);

  for (long long e=0; e<E; ++e) {
    const int i_dpe = I[e];
    const int j_dpe = J[e];
    const int owner = OWN[e];

    // neighbor is the owner of the REVERSE edge (j->i), if present
    int neigh = 0;
    auto it2 = owner_of.find( key2(j_dpe, i_dpe) );
    if (it2 != owner_of.end()) neigh = it2->second;

    // Optional: if reverse maps to same owner, treat as boundary
    if (neigh == owner) neigh = 0;

    // Record off-diagonal slot (owner, neigh)
    I2.push_back(owner);
    J2.push_back(neigh);

    // Compute lij depending on inside/outside
    double p1x = dpe(i_dpe-1,0), p1y = dpe(i_dpe-1,1);
    double p2x = dpe(j_dpe-1,0), p2y = dpe(j_dpe-1,1);
    const bool in_i = (in_poly[i_dpe-1] == TRUE);
    const bool in_j = (in_poly[j_dpe-1] == TRUE);

    double lij = 0.0;
    const int s_in = (int)in_i + (int)in_j;
    if (s_in == 2) {
      const double len = vnorm2_2d(p1x-p2x, p1y-p2y);
      lij = len * (p[i_dpe-1] + p[j_dpe-1]) / 2.0;
    } else if (s_in == 1) {
      NumericVector edge = NumericVector::create(p1x, p1y, p2x, p2y);
      NumericMatrix pi   = intersectEdgePolygon(edge, cp);
      if (pi.nrow() != 1 || pi.ncol() != 2) {
        // For robustness, if multiple points happen due to colinearity, take the closest to the inside endpoint
        if (pi.nrow() >= 1) {
          // keep the first as fallback
        } else {
          stop("calculate_hessian: edge expected to intersect boundary at exactly one point.");
        }
      }
      const double qx = pi(0,0), qy = pi(0,1);
      if (in_i) {
        const double len = vnorm2_2d(qx-p1x, qy-p1y);
        const double sig = sigma_at_point(sigma, qx, qy);
        lij = len * (sig + p[i_dpe-1]) / 2.0;
      } else {
        const double len = vnorm2_2d(qx-p2x, qy-p2y);
        const double sig = sigma_at_point(sigma, qx, qy);
        lij = len * (sig + p[j_dpe-1]) / 2.0;
      }
    } else {
      lij = 0.0;
    }

    double denom = 0.0;
    if (owner>0 && neigh>0) {
      const double dx = uv(owner-1,0) - uv(neigh-1,0);
      const double dy = uv(owner-1,1) - uv(neigh-1,1);
      denom = vnorm2_2d(dx, dy);
    }
    const double val = (denom==0.0) ? 0.0 : (-lij / denom);

    X2.push_back(val);
    rowSum[owner-1] += val;
  }

  // Add diagonal entries: H[ii,ii] += -rowSums(H)
  for (int i=0;i<nc;++i) {
    I2.push_back(i+1);
    J2.push_back(i+1);
    X2.push_back(-rowSum[i]);
  }

  // Build sparse H and symmetrize: (H + t(H)) / 2
  Function sparseMatrix("sparseMatrix");
  SEXP H = sparseMatrix(_["i"]=wrap(I2), _["j"]=wrap(J2), _["x"]=wrap(X2),
                        _["dims"]=IntegerVector::create(nc, nc));
  Function transpose("t");
  SEXP Ht = transpose(H);
  Function base_plus("+");
  SEXP Hsum = base_plus(H, Ht);
  Function base_mul("*");
  NumericVector half(1); half[0]=0.5;
  SEXP Hsym = base_mul(Hsum, half);

  return Hsym;
}
