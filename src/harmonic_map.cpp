// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <complex>

#include "compute_bd.h"
#include "vertex_area.h"

using namespace Rcpp;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

// ---- helper declarations ----

// squared distance between two vertex rows (1-based)
double sqdist_row(const Rcpp::NumericMatrix& V, int ai, int bi);

double nonneg(double x) { return x < 0.0 ? 0.0 : x; };

// =====================
// Laplace–Beltrami (cotangent)
// =====================

// Compute cotangent at vertex p with neighbors q,r using law of cosines; robust
inline double cot_from_pts(const NumericMatrix& V, int p, int q, int r) {
  // p,q,r are 1-based vertex indices
  const double a = std::sqrt(sqdist_row(V, q, r)); // |q-r|
  const double b = std::sqrt(sqdist_row(V, r, p)); // |r-p|
  const double c = std::sqrt(sqdist_row(V, p, q)); // |p-q|
  const double denom_bc = std::max(b*c, std::numeric_limits<double>::epsilon());
  const double cs = (b*b + c*c - a*a) / (2.0 * denom_bc);
  double cssq = cs*cs;
  if (cssq > 1.0) cssq = 1.0;
  const double ss = std::sqrt(std::max(0.0, 1.0 - cssq));
  if (ss <= std::numeric_limits<double>::epsilon()) return 0.0;
  return cs / ss;
}

struct EdgeKey { int a,b; };
struct EdgeKeyHash {
  std::size_t operator()(EdgeKey const& e) const noexcept {
    return (static_cast<std::size_t>(e.a) << 32) ^ static_cast<std::size_t>(e.b);
  }
};
struct EdgeKeyEq {
  bool operator()(EdgeKey const& x, EdgeKey const& y) const noexcept {
    return x.a==y.a && x.b==y.b;
  }
};

// [[Rcpp::export]]
SEXP laplace_beltrami(const IntegerMatrix& face,
                          const NumericMatrix& vertex,
                          std::string method = "Polthier") {
  const int nf = face.nrow();
  const int nv = vertex.nrow();
  if (face.ncol() != 3) stop("face must be n x 3.");
  if (vertex.ncol() < 2) stop("vertex must have dimension d >= 2.");

  std::transform(method.begin(), method.end(), method.begin(),
                 [](unsigned char c){ return std::tolower(c); });

  // --- accumulate edge weights ew for each undirected edge (i<j)
  std::unordered_map<EdgeKey, double, EdgeKeyHash, EdgeKeyEq> ew;
  ew.reserve(nf * 3);

  for (int f = 0; f < nf; ++f) {
    int i = face(f,0), j = face(f,1), k = face(f,2);

    // Edge (i,j), opposite vertex k
    {
      int a = std::min(i,j), b = std::max(i,j);
      EdgeKey key{a,b};
      ew[key] += cot_from_pts(vertex, k, i, j);
    }
    // Edge (j,k), opposite vertex i
    {
      int a = std::min(j,k), b = std::max(j,k);
      EdgeKey key{a,b};
      ew[key] += cot_from_pts(vertex, i, j, k);
    }
    // Edge (k,i), opposite vertex j
    {
      int a = std::min(k,i), b = std::max(k,i);
      EdgeKey key{a,b};
      ew[key] += cot_from_pts(vertex, j, k, i);
    }
  }

  // Build sparse L according to method
  std::vector<Triplet<double>> trips;
  trips.reserve(ew.size()*4 + nv); // symmetric off-diags + diagonals

  // For Meyer/Desbrun we need per-vertex area
  NumericVector va;
  if (method == "meyer") {
    va = vertex_area(face, vertex, "mixed");
  } else if (method == "desbrun") {
    va = vertex_area(face, vertex, "one_ring"); // then *3 later
  }

  VectorXd diag = VectorXd::Zero(nv);

  for (auto &kv : ew) {
    int i = kv.first.a;
    int j = kv.first.b;
    double w = kv.second;

    if (method == "polthier") {
      double wij = w / 2.0;
      trips.emplace_back(i-1, j-1, -wij);
      trips.emplace_back(j-1, i-1, -wij);
      diag[i-1] += wij;
      diag[j-1] += wij;

    } else if (method == "meyer") {
      double wij = (w / std::max(va[i-1],  std::numeric_limits<double>::epsilon())
                  + w / std::max(va[j-1],  std::numeric_limits<double>::epsilon())) / 2.0;
      trips.emplace_back(i-1, j-1, -wij);
      trips.emplace_back(j-1, i-1, -wij);
      diag[i-1] += wij;
      diag[j-1] += wij;

    } else if (method == "desbrun") {
      double wij = (w / std::max(va[i-1],  std::numeric_limits<double>::epsilon())
                  + w / std::max(va[j-1],  std::numeric_limits<double>::epsilon())) / 2.0 * 3.0;
      trips.emplace_back(i-1, j-1, -wij);
      trips.emplace_back(j-1, i-1, -wij);
      diag[i-1] += wij;
      diag[j-1] += wij;

    } else {
      stop("laplace_beltrami: unknown method. Use 'Polthier','Meyer','Desbrun'.");
    }
  }

  for (int i = 0; i < nv; ++i) {
    if (diag[i] != 0.0) trips.emplace_back(i, i, diag[i]);
  }

  SparseMatrix<double> L(nv, nv);
  L.setFromTriplets(trips.begin(), trips.end());
  L.makeCompressed();

  return Rcpp::wrap(L); // dgCMatrix
}

// =====================
// Rectangular harmonic map (corners -> unit square OR explicit UVs)
// =====================

// helper to rotate a closed boundary so it starts at 'firstCorner'
static inline void rotate_closed_boundary(IntegerVector& bd, int firstCorner) {
  int nbd = bd.size();
  if (nbd == 0) return;
  int pos = -1;
  for (int i = 0; i < nbd; ++i) if (bd[i] == firstCorner) { pos = i; break; }
  if (pos >= 0 && pos != 0) {
    IntegerVector bd2(nbd);
    for (int t = 0; t < nbd; ++t) bd2[t] = bd[(pos + t) % nbd];
    bd = bd2;
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rect_harmonic_map(const Rcpp::IntegerMatrix& face,
                                      const Rcpp::NumericMatrix& vertex,
                                      SEXP corner,
                                      std::string method = "Polthier") {

  const int nv = vertex.nrow();
  if (face.ncol() != 3) stop("face must be n x 3 (triangles).");
  if (vertex.ncol() < 2) stop("vertex must have dimension d >= 2.");

  // 1) boundary (closed)
  IntegerVector bd = compute_bd(face, false);
  if (bd.size() < 4) stop("Boundary has fewer than 4 vertices.");
  // make bd closed (last==first)
  if (bd[bd.size()-1] != bd[0]) bd.push_back(bd[0]);

  // 2) parse 'corner'
  IntegerVector corner_idx;
  NumericMatrix corner_uv;
  bool explicit_uv = false;

  if (Rf_isMatrix(corner)) {
    NumericMatrix C(corner);
    if (C.nrow() != 4 || C.ncol() != 3)
      stop("When 'corner' is a matrix, it must be 4x3: [index, u, v].");
    corner_idx = IntegerVector(4);
    for (int i = 0; i < 4; ++i) corner_idx[i] = (int)std::round(C(i,0));
    corner_uv = NumericMatrix(4,2);
    for (int i = 0; i < 4; ++i){ corner_uv(i,0) = C(i,1); corner_uv(i,1) = C(i,2); }
    explicit_uv = true;
  } else {
    IntegerVector idx = as<IntegerVector>(corner);
    if (idx.size() != 4) stop("When 'corner' is a vector, it must contain 4 indices.");
    corner_idx = clone(idx);
    explicit_uv = false;
  }

  // 3) rotate boundary to start at first corner
  rotate_closed_boundary(bd, corner_idx[0]);

  // ensure bd is still closed
  if (bd[bd.size()-1] != bd[0]) bd.push_back(bd[0]);

  // positions of corners along bd
  IntegerVector ck(4);
  for (int k = 0; k < 4; ++k) {
    int pos = -1;
    for (int i = 0; i < bd.size(); ++i) if (bd[i] == corner_idx[k]) { pos = i+1; break; }
    if (pos < 0) stop("A corner index was not found on boundary.");
    ck[k] = pos;
  }
  // verify increasing order along rotated boundary
  if (!(ck[0]<ck[1] && ck[1]<ck[2] && ck[2]<ck[3]))
    stop("Corner indices must be provided in boundary order.");

  const int nbd = bd.size();
  const int ck5 = nbd;

  // 4) assemble UV boundary values
  NumericMatrix uv(nv, 2);
  NumericMatrix uvbd(nbd, 2);

  if (explicit_uv) {
    std::vector<std::complex<double>> zc(4);
    for (int i = 0; i < 4; ++i) zc[i] = std::complex<double>(corner_uv(i,0), corner_uv(i,1));
    std::vector<std::complex<double>> zbd(nbd);
    auto fillSeg = [&](int s, int e, std::complex<double> A, std::complex<double> B){
      int count = e - s + 1;
      for (int t = 0; t < count; ++t) {
        double a = count==1?0.0:(double)t/(count-1);
        zbd[s-1 + t] = (1.0 - a)*A + a*B;
      }
    };
    fillSeg(ck[0], ck[1], zc[0], zc[1]);
    fillSeg(ck[1], ck[2], zc[1], zc[2]);
    fillSeg(ck[2], ck[3], zc[2], zc[3]);
    fillSeg(ck[3], ck5,  zc[3], zc[0]);
    for (int i = 0; i < nbd; ++i) { uvbd(i,0)=zbd[i].real(); uvbd(i,1)=zbd[i].imag(); }
  } else {
    // unit square boundary
    auto fillSeg = [&](int s, int e, double x0,double y0,double x1,double y1){
      int count = e - s + 1;
      for (int t = 0; t < count; ++t) {
        double a = count==1?0.0:(double)t/(count-1);
        uvbd(s-1+t,0) = (1.0-a)*x0 + a*x1;
        uvbd(s-1+t,1) = (1.0-a)*y0 + a*y1;
      }
    };
    // (0,0)->(1,0)->(1,1)->(0,1)->(0,0)
    fillSeg(ck[0], ck[1], 0,0, 1,0);
    fillSeg(ck[1], ck[2], 1,0, 1,1);
    fillSeg(ck[2], ck[3], 1,1, 0,1);
    fillSeg(ck[3], ck5,   0,1,  0,0);
  }

  // drop duplicated last boundary vertex
  NumericMatrix uvbd_use(nbd-1,2);
  IntegerVector bd_use(nbd-1);
  for (int i = 0; i < nbd-1; ++i) {
    uvbd_use(i,0) = uvbd(i,0);
    uvbd_use(i,1) = uvbd(i,1);
    bd_use[i]     = bd[i];
  }
  for (int i = 0; i < nbd-1; ++i) uv(bd_use[i]-1,0)=uvbd_use(i,0), uv(bd_use[i]-1,1)=uvbd_use(i,1);

  // 5) interior Laplace solve: L(in,in) * X = - L(in,bd) * UVbd
  SEXP L_sexp = laplace_beltrami(face, vertex, method);
  SparseMatrix<double> L = Rcpp::as<SparseMatrix<double>>(L_sexp);
  L.makeCompressed();

  std::vector<char> is_inside(nv, 1);
  for (int r = 0; r < bd_use.size(); ++r) is_inside[bd_use[r]-1] = 0;

  int n_in = 0;
  std::vector<int> map_row(nv, -1);
  for (int i = 0; i < nv; ++i) if (is_inside[i]) map_row[i] = n_in++;

  std::vector<Triplet<double>> tin;
  tin.reserve(L.nonZeros());
  MatrixXd rhs(n_in, 2);
  rhs.setZero();

  for (int k = 0; k < L.outerSize(); ++k) {
    for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      const int i = it.row();
      const int j = it.col();
      const double v = it.value();
      if (is_inside[i]) {
        if (is_inside[j]) {
          tin.emplace_back(map_row[i], map_row[j], v);
        } else {
          rhs(map_row[i], 0) -= v * uv(j,0);
          rhs(map_row[i], 1) -= v * uv(j,1);
        }
      }
    }
  }

  SparseMatrix<double> Ain(n_in, n_in);
  Ain.setFromTriplets(tin.begin(), tin.end());
  Ain.makeCompressed();

  Eigen::SimplicialLLT<SparseMatrix<double>> solver;
  solver.compute(Ain);
  if (solver.info() != Eigen::Success) stop("LLT factorization failed for interior Laplacian.");
  MatrixXd uvin = solver.solve(rhs);
  if (solver.info() != Eigen::Success) stop("Linear solve failed for interior Laplacian.");

  for (int i = 0; i < nv; ++i) if (is_inside[i]) { int r = map_row[i]; uv(i,0)=uvin(r,0); uv(i,1)=uvin(r,1); }

  return uv;
}

// =====================
// Disk harmonic map with a fixed boundary reference point
// =====================

// [[Rcpp::export]]
NumericMatrix disk_harmonic_map(const IntegerMatrix& face,
                                            const NumericMatrix& vertex,
                                            int fixed_point,
                                            std::string method = "Polthier") {
  const int nv = vertex.nrow();
  if (nv == 0) stop("Empty vertex.");
  if (fixed_point < 1 || fixed_point > nv) stop("fixed_point must be a valid vertex index (1..nV).");

  // 1) Boundary
  IntegerVector bd = compute_bd(face, false);
  if (bd.size() == 0) stop("No boundary detected.");

  // Rotate bd to start at fixed_point if present
  {
    int pos = -1;
    for (int i = 0; i < bd.size(); ++i) if (bd[i] == fixed_point) { pos = i; break; }
    if (pos >= 0 && pos != 0) {
      IntegerVector bd2(bd.size());
      int m = bd.size();
      for (int t = 0; t < m; ++t) bd2[t] = bd[(pos + t) % m];
      bd = bd2;
    }
  }

  // 2) centroid: nearest vertex to mean of boundary vertices
  std::vector<double> mean( vertex.ncol(), 0.0 );
  for (int r = 0; r < bd.size(); ++r) {
    int idx = bd[r] - 1;
    for (int k = 0; k < vertex.ncol(); ++k) mean[k] += vertex(idx, k);
  }
  for (int k = 0; k < (int)mean.size(); ++k) mean[k] /= (double)bd.size();

  int centroid = 1;
  double best = std::numeric_limits<double>::infinity();
  for (int i = 0; i < nv; ++i) {
    double s = 0.0;
    for (int k = 0; k < vertex.ncol(); ++k) {
      double d = vertex(i,k) - mean[k];
      s += d*d;
    }
    if (s < best) { best = s; centroid = i+1; }
  }

  // 3) map boundary to unit circle proportional to boundary edge length
  IntegerVector next_bd(bd.size());
  for (int i = 0; i < bd.size()-1; ++i) next_bd[i] = bd[i+1];
  next_bd[bd.size()-1] = bd[0];

  std::vector<double> bl(bd.size());
  double bl_sum = 0.0;
  for (int i = 0; i < bd.size(); ++i) {
    bl[i] = std::sqrt(sqdist_row(vertex, bd[i], next_bd[i]));
    bl_sum += bl[i];
  }

  std::vector<double> t(bd.size());
  double acc = 0.0;
  for (int i = 0; i < bd.size(); ++i) {
    acc += bl[i];
    t[i] = acc / bl_sum * 2.0 * M_PI;
  }
  // shift: t([end,1:end-1])
  if (!t.empty()) {
    double last = t.back();
    for (int i = (int)t.size()-1; i >= 1; --i) t[i] = t[i-1];
    t[0] = last;
  }

  NumericMatrix uv(nv, 2);
  for (int r = 0; r < bd.size(); ++r) {
    uv(bd[r]-1, 0) = std::cos(t[r]);
    uv(bd[r]-1, 1) = std::sin(t[r]);
  }

  // 4) interior solve: L_in,in * uv_in = - L_in,bd * uv_bd
  SEXP L_sexp = laplace_beltrami(face, vertex, method);
  SparseMatrix<double> L = Rcpp::as< SparseMatrix<double> >(L_sexp);
  L.makeCompressed();

  std::vector<char> is_inside(nv, 1);
  for (int r = 0; r < bd.size(); ++r) is_inside[bd[r]-1] = 0;

  int n_in = 0;
  std::vector<int> map_row(nv, -1);
  for (int i = 0; i < nv; ++i) if (is_inside[i]) map_row[i] = n_in++;

  std::vector<Triplet<double>> tin;
  tin.reserve(L.nonZeros());
  MatrixXd rhs(n_in, 2);
  rhs.setZero();

  for (int k = 0; k < L.outerSize(); ++k) {
    for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      const int i = it.row();
      const int j = it.col();
      const double v = it.value();
      if (is_inside[i]) {
        if (is_inside[j]) {
          tin.emplace_back(map_row[i], map_row[j], v);
        } else {
          rhs(map_row[i], 0) -= v * uv(j,0);
          rhs(map_row[i], 1) -= v * uv(j,1);
        }
      }
    }
  }

  SparseMatrix<double> Ain(n_in, n_in);
  Ain.setFromTriplets(tin.begin(), tin.end());
  Ain.makeCompressed();

  Eigen::SimplicialLLT<SparseMatrix<double>> solver;
  solver.compute(Ain);
  if(solver.info() != Eigen::Success) stop("LLT factorization failed for interior Laplacian.");
  MatrixXd uvin = solver.solve(rhs);
  if(solver.info() != Eigen::Success) stop("Linear solve failed for interior Laplacian.");

  for (int i = 0; i < nv; ++i) {
    if (is_inside[i]) {
      int r = map_row[i];
      uv(i,0) = uvin(r,0);
      uv(i,1) = uvin(r,1);
    }
  }

  // 5) Möbius transform: center using centroid, then rotate to orient fixed_point
  std::vector< std::complex<double> > z(nv);
  for (int i = 0; i < nv; ++i) z[i] = std::complex<double>(uv(i,0), uv(i,1));
  std::complex<double> z0 = z[centroid-1];
  std::complex<double> z0c = std::conj(z0);

  for (int i = 0; i < nv; ++i) {
    std::complex<double> num = z[i] - z0;
    std::complex<double> den = (std::complex<double>(1.0,0.0) - z0c * z[i]);
    if (std::abs(den) < 1e-15) den = std::complex<double>(1e-15, 0.0);
    z[i] = num / den;
  }

  // rotate so that fixed_point has angle -pi/2
  {
    double ang = std::arg(z[fixed_point-1]);
    std::complex<double> rot = std::polar(1.0, -(ang + M_PI/2.0));
    for (int i = 0; i < nv; ++i) z[i] *= rot;
  }

  for (int i = 0; i < nv; ++i) {
    uv(i,0) = z[i].real();
    uv(i,1) = z[i].imag();
  }

  return uv;
}
