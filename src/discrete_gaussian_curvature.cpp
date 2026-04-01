// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <unordered_map>
#include <cmath>
using namespace Rcpp;

inline double clamp(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

// Compute per-face internal angles (radians) and cosines from edge lengths.
// l: (#F x 3) where columns are [s23, s31, s12] (edges opposite v1, v2, v3)
// Returns a list: A (angles #F x 3), cA (cosines #F x 3)
// [[Rcpp::export]]
List internalangles_intrinsic(const NumericMatrix& l) {
  int F = l.nrow();
  NumericMatrix A(F, 3), cA(F, 3);

  for (int f = 0; f < F; ++f) {
    // Match MATLAB naming:
    // s23 = l(:,1); s31 = l(:,2); s12 = l(:,3);
    double s23 = l(f, 0);
    double s31 = l(f, 1);
    double s12 = l(f, 2);

    // Law of cosines (cosines at v1, v2, v3 respectively)
    // ca23 is angle at vertex 1 (opposite s23), etc.
    double ca23 = (s12*s12 + s31*s31 - s23*s23) / (2.0 * s12 * s31);
    double ca31 = (s23*s23 + s12*s12 - s31*s31) / (2.0 * s23 * s12);
    double ca12 = (s31*s31 + s23*s23 - s12*s12) / (2.0 * s31 * s23);

    ca23 = clamp(ca23, -1.0, 1.0);
    ca31 = clamp(ca31, -1.0, 1.0);
    ca12 = clamp(ca12, -1.0, 1.0);

    cA(f, 0) = ca23; cA(f, 1) = ca31; cA(f, 2) = ca12;

    A(f, 0) = std::acos(ca23);
    A(f, 1) = std::acos(ca31);
    A(f, 2) = std::acos(ca12);
  }
  return List::create(_["A"] = A, _["cA"] = cA);
}

// Helper: row-wise Euclidean norm of (v2 - v1) for each face row
inline double edge_len(const NumericMatrix& V, int i, int j) {
  // i, j are 0-based vertex indices
  int d = V.ncol();
  double s = 0.0;
  for (int k = 0; k < d; ++k) {
    double diff = V(i, k) - V(j, k);
    s += diff * diff;
  }
  return std::sqrt(s);
}

// Compute per-face internal angles (radians) from vertex coords & faces
// vertex: #V x d (d>=2; typically >=3)
// face:   #F x 3 (1-based indices)
// Returns a list: A (angles #F x 3), cA (cosines #F x 3), and l (edge lengths #F x 3)
// [[Rcpp::export]]
List internalangles(const NumericMatrix& vertex, const IntegerMatrix& face) {
  int F = face.nrow();
  NumericMatrix l(F, 3);

  for (int f = 0; f < F; ++f) {
    // 1-based -> 0-based
    int i1 = face(f, 0) - 1;
    int i2 = face(f, 1) - 1;
    int i3 = face(f, 2) - 1;

    // s12 = ||v2 - v1||, s13 = ||v3 - v1||, s23 = ||v3 - v2||
    double s12 = edge_len(vertex, i2, i1);
    double s13 = edge_len(vertex, i3, i1);
    double s23 = edge_len(vertex, i3, i2);

    // MATLAB: l = [s23 s13 s12];
    l(f, 0) = s23;
    l(f, 1) = s13;
    l(f, 2) = s12;
  }

  return internalangles_intrinsic(l);
}

// Find boundary vertices (unique of boundary edges) from triangle mesh
// Using undirected edge counts: those with count==1 are boundary edges.
// Returns a sorted unique integer vector of 1-based vertex ids.
// [[Rcpp::export]]
IntegerVector outline_vertices(const IntegerMatrix& face) {
  int F = face.nrow();

  struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
      // simple hash combine
      return std::hash<int>()(p.first * 73856093) ^ std::hash<int>()(p.second * 19349663);
    }
  };

  std::unordered_map<std::pair<int,int>, int, PairHash> counts;
  counts.reserve(static_cast<size_t>(F) * 3);

  auto add_edge = [&](int a, int b){
    // store as (min,max); inputs a,b are 1-based
    if (a > b) std::swap(a, b);
    std::pair<int,int> e(a,b);
    auto it = counts.find(e);
    if (it == counts.end()) counts.emplace(e, 1);
    else it->second += 1;
  };

  for (int f = 0; f < F; ++f) {
    int i1 = face(f, 0);
    int i2 = face(f, 1);
    int i3 = face(f, 2);
    add_edge(i1, i2);
    add_edge(i2, i3);
    add_edge(i3, i1);
  }

  // Collect unique vertices that appear in edges with count==1
  std::unordered_set<int> bv;
  bv.reserve(F); // rough
  for (const auto& kv : counts) {
    if (kv.second == 1) {
      bv.insert(kv.first.first);
      bv.insert(kv.first.second);
    }
  }

  // Convert to sorted IntegerVector
  IntegerVector out(bv.begin(), bv.end());
  std::sort(out.begin(), out.end());
  return out;
}

// Main function: discrete Gaussian curvature per vertex (radians).
// Implements K_G(i) = 2*pi - sum_{faces incident to i} (angle at i),
// and subtracts an additional pi for boundary vertices (Meyer et al. 2002, eq. 9 w/o inverse area).
// vertex: #V x d (d>=2; typically >=3)
// face:   #F x 3 (1-based indices)
// Returns NumericVector k (#V)
// [[Rcpp::export]]
NumericVector discrete_gaussian_curvature(const NumericMatrix& vertex,
                                          const IntegerMatrix& face) {
  int Vn = vertex.nrow();
  int F = face.nrow();

  // 1) Per-face angles
  List ia = internalangles(vertex, face);
  NumericMatrix A = ia["A"]; // #F x 3, radians

  // 2) Accumulate sum of angles per vertex
  NumericVector angsum(Vn, 0.0);
  for (int f = 0; f < F; ++f) {
    int i1 = face(f, 0) - 1; // 0-based
    int i2 = face(f, 1) - 1;
    int i3 = face(f, 2) - 1;

    // A(f,0),A(f,1),A(f,2) correspond to angles at vertices (i1,i2,i3)
    angsum[i1] += A(f, 0);
    angsum[i2] += A(f, 1);
    angsum[i3] += A(f, 2);
  }

  // 3) K = 2*pi - sum angles
  const double two_pi = 2.0 * M_PI;
  NumericVector k(Vn);
  for (int i = 0; i < Vn; ++i) {
    k[i] = two_pi - angsum[i];
  }

  // 4) Boundary correction: subtract pi at boundary vertices
  IntegerVector bverts = outline_vertices(face);
  for (int idx = 0; idx < bverts.size(); ++idx) {
    int vi = bverts[idx] - 1; // 0-based index
    if (vi >= 0 && vi < Vn) k[vi] -= M_PI;
  }

  return k;
}
