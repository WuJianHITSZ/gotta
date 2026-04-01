// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace Rcpp;

inline long long pack_edge(int i, int j) {
  return ( (static_cast<long long>(i) << 32) |
           (static_cast<unsigned int>(j)) );
}

static IntegerVector sorted_unique_vertices(const IntegerMatrix& face, int& maxv) {
  std::vector<int> vals;
  vals.reserve(static_cast<size_t>(face.nrow()) * static_cast<size_t>(face.ncol()));
  maxv = 0;
  for (int r = 0; r < face.nrow(); ++r) {
    for (int c = 0; c < face.ncol(); ++c) {
      int v = face(r, c);
      if (v > maxv) maxv = v;
      vals.push_back(v);
    }
  }
  std::sort(vals.begin(), vals.end());
  vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
  return IntegerVector(vals.begin(), vals.end());
}

static IntegerMatrix remap_face(const IntegerMatrix& face, const IntegerVector& vertex) {
  std::unordered_map<int,int> mp;
  mp.reserve(vertex.size() * 2);
  for (int i = 0; i < vertex.size(); ++i) mp[ vertex[i] ] = i + 1; // 1-based
  IntegerMatrix out(face.nrow(), face.ncol());
  for (int r = 0; r < face.nrow(); ++r)
    for (int c = 0; c < face.ncol(); ++c)
      out(r, c) = mp.at(face(r, c));
  return out;
}

// [[Rcpp::export]]
IntegerVector compute_bd(const IntegerMatrix& face_in, bool efficient = false) {
  if (face_in.ncol() != 3 || face_in.nrow() < 1)
    stop("face must be an n x 3 integer matrix with n >= 1.");

  int maxv_face = 0;
  IntegerVector vertex = sorted_unique_vertices(face_in, maxv_face);
  bool do_efficient = efficient || (vertex.size() < maxv_face);
  IntegerMatrix face = do_efficient ? remap_face(face_in, vertex) : IntegerMatrix(clone(face_in));
  if (do_efficient) maxv_face = vertex.size();

  // Oriented adjacency counts A(i,j)
  std::unordered_map<long long, int> A;
  A.reserve(static_cast<size_t>(face.nrow()) * 6);
  for (int r = 0; r < face.nrow(); ++r) {
    int a = face(r,0), b = face(r,1), c = face(r,2);
    ++A[pack_edge(a,b)];
    ++A[pack_edge(b,c)];
    ++A[pack_edge(c,a)];
  }

  // boundary_table in COLUMN-MAJOR order (like R's which(..., arr.ind=TRUE)):
  // iterate j = 1..max, then i = 1..max, include if A(i,j) exists AND A(i,j)-A(j,i) > 0
  std::vector< std::pair<int,int> > boundary_table;
  boundary_table.reserve(A.size());
  for (int j = 1; j <= maxv_face; ++j) {
    for (int i = 1; i <= maxv_face; ++i) {
      auto it = A.find(pack_edge(i,j));
      if (it == A.end()) continue;
      int a = it->second;
      auto rev = A.find(pack_edge(j,i));
      int b = (rev == A.end()) ? 0 : rev->second;
      if ((a - b) > 0) boundary_table.emplace_back(i, j);
    }
  }

  if (boundary_table.empty()) return IntegerVector(0);

  // Walk exactly like the R loop
  std::vector<int> boundary(boundary_table.size() + 2, NA_INTEGER);
  boundary[0] = boundary_table[0].first;

  auto seen_contains = [&](int val, size_t upto)->bool{
    for (size_t k = 0; k < upto; ++k) if (boundary[k] == val) return true;
    return false;
  };

  size_t count = 0;
  for (count = 1; count <= boundary_table.size(); ++count) {
    int current = boundary[count - 1];

    int pointer_idx = -1;
    for (size_t k = 0; k < boundary_table.size(); ++k) {
      if (boundary_table[k].first == current) { pointer_idx = (int)k; break; }
    }
    if (pointer_idx < 0) break;

    boundary[count] = boundary_table[pointer_idx].second;
    if (seen_contains(boundary[count], count)) break;
  }

  // boundary[seq_len(count)+1]
  IntegerVector out( (int)count );
  for (size_t i = 0; i < count; ++i) out[i] = boundary[i + 1];

  if (do_efficient) {
    for (int i = 0; i < out.size(); ++i) {
      int idx = out[i];
      if (idx < 1 || idx > vertex.size()) stop("map-back index out of range.");
      out[i] = vertex[idx - 1];
    }
  }
  return out;
}
