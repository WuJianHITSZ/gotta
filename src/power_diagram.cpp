// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "compute_bd.h"

using namespace Rcpp;

/*** ========================= Utilities ========================= ***/
static inline double cross2_xy(double ax, double ay, double bx, double by) {
  return ax*by - ay*bx;
}

/*** ==================== compute_bd (boundary) ==================== ***/
// Encode a directed edge (i -> j)
static inline long long edge_key(int i, int j, int maxN){
  return (static_cast<long long>(i) * (static_cast<long long>(maxN) + 1LL))
       + static_cast<long long>(j);
}

/*** ================== compute_vertex_ring (one-ring) ================== ***/
// [[Rcpp::export]]
List compute_vertex_ring(IntegerMatrix face) {
  const int n = face.nrow();
  if (face.ncol() != 3) stop("compute_vertex_ring expects a triangular face matrix with 3 columns.");

  IntegerVector global_bd = compute_bd(face, false);
  std::unordered_set<int> bdset; for (int v : global_bd) bdset.insert(v);

  int nv = 0;
  for (int i = 0; i < n; ++i) for (int c = 0; c < 3; ++c) if (face(i,c) > nv) nv = face(i,c);

  std::vector<std::vector<int>> rows_by_v(nv+1);
  for (int r = 0; r < n; ++r) {
    rows_by_v[face(r,0)].push_back(r);
    rows_by_v[face(r,1)].push_back(r);
    rows_by_v[face(r,2)].push_back(r);
  }

  List rings(nv);
  for (int idx = 1; idx <= nv; ++idx) {
    // facet = rows where idx appears exactly once
    std::vector<int> rows;
    for (int r : rows_by_v[idx]) {
      int cnt = (face(r,0)==idx) + (face(r,1)==idx) + (face(r,2)==idx);
      if (cnt == 1) rows.push_back(r);
    }
    IntegerMatrix facet((int)rows.size(), 3);
    for (int i = 0; i < (int)rows.size(); ++i) {
      int rr = rows[i];
      facet(i,0) = face(rr,0); facet(i,1) = face(rr,1); facet(i,2) = face(rr,2);
    }

    IntegerVector ring;
    if (bdset.find(idx) != bdset.end()) {
      // boundary vertex
      IntegerVector tmp = compute_bd(facet, true);
      std::vector<int> vr(tmp.begin(), tmp.end());
      int pos = -1;
      for (int k = 0; k < (int)vr.size(); ++k) if (vr[k]==idx){pos=k;break;}
      if (pos == -1) {
        ring = tmp;
      } else if (pos == 0 || pos == (int)vr.size()-1) {
        std::vector<int> filtered; for (int v : vr) if (v!=idx) filtered.push_back(v);
        ring = wrap(filtered);
      } else {
        std::vector<int> rot;
        for (size_t t = pos+1; t < vr.size(); ++t) rot.push_back(vr[t]);
        for (int t = 0; t <= pos; ++t)      rot.push_back(vr[t]);
        ring = wrap(rot);
      }
    } else {
      // interior vertex
      IntegerVector tmp = compute_bd(facet, false);
      if (tmp.size()==0) ring = tmp;
      else {
        IntegerVector closed(tmp.size()+1);
        for (int i = 0; i < tmp.size(); ++i) closed[i]=tmp[i];
        closed[tmp.size()] = tmp[0];
        ring = closed;
      }
    }

    // clean NAs (defensive)
    std::vector<int> cleaned; cleaned.reserve(ring.size());
    for (int i = 0; i < ring.size(); ++i) if (!IntegerVector::is_na(ring[i])) cleaned.push_back(ring[i]);
    rings[idx-1] = wrap(cleaned);
  }
  return rings;
}

/*** ===================== calculate_face_normal ===================== ***/
// signed 2D triangle area (right-hand). Return vector of per-face z (scalar).
// [[Rcpp::export]]
NumericVector calculate_face_normal(IntegerMatrix face, NumericMatrix uv){
  const int nf = face.nrow();
  NumericVector out(nf);
  for (int r = 0; r < nf; ++r) {
    int i = face(r,0)-1, j = face(r,1)-1, k = face(r,2)-1;
    double ax = uv(j,0) - uv(i,0), ay = uv(j,1) - uv(i,1);
    double bx = uv(k,0) - uv(i,0), by = uv(k,1) - uv(i,1);
    out[r] = cross2_xy(ax,ay,bx,by);
  }
  return out;
}

/*** ========================= face_dual_uv ========================= ***/
// p: 3x3 matrix rows are vertices, cols (x,y,z). Returns length-2 dp.
// [[Rcpp::export]]
NumericVector face_dual_uv(NumericMatrix p){
  if (p.nrow()!=3 || p.ncol()!=3) stop("face_dual_uv expects 3x3 matrix");
  double a = p(0,1)*(p(1,2)-p(2,2)) + p(1,1)*(p(2,2)-p(0,2)) + p(2,1)*(p(0,2)-p(1,2));
  double b = p(0,2)*(p(1,0)-p(2,0)) + p(1,2)*(p(2,0)-p(0,0)) + p(2,2)*(p(0,0)-p(1,0));
  double c  = p(0,0)*(p(1,1)-p(2,1)) + p(1,0)*(p(2,1)-p(0,1)) + p(2,0)*(p(0,1)-p(1,1));
  double ux = -a / c / 2.0;
  double uy = -b / c / 2.0;
  return NumericVector::create(ux, uy);
}

/*** ====================== compute_connectivity ====================== ***/
// Returns a dgCMatrix (vvif) via Matrix::sparseMatrix(i=..., j=..., x=ff)
// vvif[i,j] = face index where directed edge (i,j) occurs.
// [[Rcpp::export]]
SEXP compute_connectivity(IntegerMatrix face){
  const int nf = face.nrow();
  if (face.ncol() != 3) stop("face must be nf x 3");

  IntegerVector fi(nf), fj(nf), fk(nf), ff(nf);
  int nv = 0;
  for (int r = 0; r < nf; ++r) {
    fi[r] = face(r,0);
    fj[r] = face(r,1);
    fk[r] = face(r,2);
    ff[r] = r + 1;              // face id
    nv = std::max(nv, face(r,0));
    nv = std::max(nv, face(r,1));
    nv = std::max(nv, face(r,2));
  }

  // i = c(fi,fj,fk); j = c(fj,fk,fi); x = c(ff,ff,ff)
  IntegerVector i(3*nf), j(3*nf), x(3*nf);
  i[ Range(0,        nf-1) ] = fi;
  i[ Range(nf,     2*nf-1) ] = fj;
  i[ Range(2*nf,   3*nf-1) ] = fk;

  j[ Range(0,        nf-1) ] = fj;
  j[ Range(nf,     2*nf-1) ] = fk;
  j[ Range(2*nf,   3*nf-1) ] = fi;

  x[ Range(0,        nf-1) ] = ff;
  x[ Range(nf,     2*nf-1) ] = ff;
  x[ Range(2*nf,   3*nf-1) ] = ff;

  Function sparseMatrix("sparseMatrix"); // from Matrix package
  return sparseMatrix(
    _["i"]    = i,
    _["j"]    = j,
    _["x"]    = x,
    _["dims"] = IntegerVector::create(nv, nv)
  );
}

/*** ===================== intersectRayPolygon ===================== ***/
static inline std::string point_key(double x, double y, double tol){
  long long kx = llround(x / tol);
  long long ky = llround(y / tol);
  return std::to_string(kx) + "_" + std::to_string(ky);
}

struct Hit2D { double x,y,t; };

static NumericMatrix intersect_core(const NumericVector &O,
                                    const NumericVector &D,
                                    const NumericMatrix &poly,
                                    double tol)
{
  int N = poly.nrow();
  if (N < 2) return NumericMatrix(0,2);
  double Dn = std::sqrt(D[0]*D[0] + D[1]*D[1]);
  if (Dn < tol) return NumericMatrix(0,2);
  NumericVector Dunit = NumericVector::create(D[0]/Dn, D[1]/Dn);

  std::vector<Hit2D> hits; hits.reserve(N);
  for (int i = 0; i < N; ++i) {
    NumericVector P = poly(i,_);
    NumericVector Q = poly((i+1<N)?(i+1):0, _);
    double Sx = Q[0]-P[0], Sy = Q[1]-P[1];
    double denom = cross2_xy(D[0],D[1], Sx,Sy);

    if (std::fabs(denom) > tol) {
      double wx = P[0]-O[0], wy = P[1]-O[1];
      double t = cross2_xy(wx,wy, Sx,Sy) / denom;
      double u = cross2_xy(wx,wy, D[0],D[1]) / denom;
      if (t >= -tol && u >= -tol && u <= 1.0+tol) {
        t = std::max(0.0, t);
        u = std::max(0.0, std::min(1.0, u));
        hits.push_back({O[0]+t*D[0], O[1]+t*D[1], t});
      }
    } else {
      double wx = P[0]-O[0], wy = P[1]-O[1];
      if (std::fabs(cross2_xy(wx,wy, D[0],D[1])) <= tol) {
        double tP = ( (P[0]-O[0])*Dunit[0] + (P[1]-O[1])*Dunit[1] );
        double tQ = ( (Q[0]-O[0])*Dunit[0] + (Q[1]-O[1])*Dunit[1] );
        double tmin = std::max(0.0, std::min(tP,tQ));
        double tmax = std::max(0.0, std::max(tP,tQ));
        if (tmax >= tmin - tol) {
          double x = O[0] + tmin*Dunit[0];
          double y = O[1] + tmin*Dunit[1];
          hits.push_back({x,y, tmin*Dn});
        }
      }
    }
  }
  if (hits.empty()) return NumericMatrix(0,2);

  std::unordered_set<std::string> seen; seen.reserve(hits.size()*2);
  std::vector<Hit2D> uniq; uniq.reserve(hits.size());
  for (const auto& h : hits) {
    std::string k = point_key(h.x,h.y,1e-12);
    if (seen.insert(k).second) uniq.push_back(h);
  }
  std::sort(uniq.begin(), uniq.end(), [](const Hit2D& a, const Hit2D& b){ return a.t < b.t; });

  NumericMatrix out((int)uniq.size(),2);
  for (int i=0;i<(int)uniq.size();++i){ out(i,0)=uniq[i].x; out(i,1)=uniq[i].y; }
  colnames(out) = CharacterVector::create("x","y");
  return out;
}

// [[Rcpp::export]]
NumericMatrix intersectRayPolygon(Nullable<NumericVector> ray,
                                      NumericMatrix poly,
                                      Nullable<NumericVector> rayOrigin = R_NilValue,
                                      Nullable<NumericVector> rayDirection = R_NilValue,
                                      double tol = 1e-12)
{
  if (poly.ncol()!=2) stop("poly must be N x 2");
  NumericVector O(2), D(2);
  if (ray.isNotNull()) {
    NumericVector r(ray); if (r.size()!=4) stop("ray must be numeric(4)");
    O[0]=r[0]; O[1]=r[1]; D[0]=r[2]; D[1]=r[3];
  } else {
    if (rayOrigin.isNull() || rayDirection.isNull()) stop("provide ray or (origin, direction)");
    NumericVector ro(rayOrigin), rd(rayDirection);
    if (ro.size()!=2 || rd.size()!=2) stop("rayOrigin/rayDirection must be numeric(2)");
    O=ro; D=rd;
  }
  return intersect_core(O,D,poly,tol);
}

/*** ========================= power_diagram ========================= ***/
// [[Rcpp::export]]
List power_diagram(IntegerMatrix face,
                       NumericMatrix uv,
                       Nullable<NumericVector> h_in = R_NilValue,
                       Nullable<NumericVector> dh_in = R_NilValue)
{
  const int nf_prev = face.nrow();
  const int n = uv.nrow();

  NumericVector h = h_in.isNotNull() ? NumericVector(h_in) : NumericVector(n, 0.0);
  NumericVector dh = dh_in.isNotNull() ? NumericVector(dh_in) : NumericVector(n, 0.0);

  // functions we call from R
  Function convhulln("convhulln"); // geometry
  Function chull("chull");         // base
  Function sparseMatrix("sparseMatrix"); // Matrix

  double cstep = 1.0;
  int iter = 0, max_iter = 50;

  auto lift_points = [&](NumericMatrix& out_pl){
    out_pl = NumericMatrix(n, 3);
    for (int i = 0; i < n; ++i) {
      double s2 = uv(i,0)*uv(i,0) + uv(i,1)*uv(i,1);
      out_pl(i,0)=uv(i,0); out_pl(i,1)=uv(i,1); out_pl(i,2)=s2 - h[i];
    }
  };

  // backtracking step to ensure enough lower-hull faces
  while (true) {
    for (int i=0;i<n;++i) h[i] -= cstep * dh[i];

    NumericMatrix pl;
    lift_points(pl);

    List cx = convhulln(pl, Named("output.options")="n");
    NumericMatrix fn = cx["normals"];
    int cnt = 0;
    for (int i=0;i<fn.nrow();++i) if (fn(i,2) < 0.0) ++cnt;

    if (cnt < nf_prev) {
      for (int i=0;i<n;++i) h[i] += cstep * dh[i];
      cstep *= 0.5;
      if (max(abs(dh)) == 0.0) break;
      if (++iter >= max_iter) break;
    } else {
      break;
    }
  }

  // final recompute
  NumericMatrix pl;
  lift_points(pl);
  List cx = convhulln(pl, Named("output.options")="n");
  IntegerMatrix face_new = cx["hull"]; // 1-based
  NumericMatrix fn = cx["normals"];
  std::vector<int> keep; keep.reserve(face_new.nrow());
  for (int i=0;i<fn.nrow();++i) if (fn(i,2) < 0.0) keep.push_back(i);
  IntegerMatrix face_filt((int)keep.size(), 3);
  for (int r=0;r<(int)keep.size();++r){
    face_filt(r,0)=face_new(keep[r],0);
    face_filt(r,1)=face_new(keep[r],1);
    face_filt(r,2)=face_new(keep[r],2);
  }

  // orient faces by 2D signed area on uv
  NumericVector areas = calculate_face_normal(face_filt, uv);
  for (int r=0; r<face_filt.nrow(); ++r) {
    if (areas[r] > 0) {
      int tmp = face_filt(r,1); face_filt(r,1) = face_filt(r,2); face_filt(r,2) = tmp;
    }
  }

  // one-ring per vertex
  List vr = compute_vertex_ring(face_filt);

  // dual points for each face (lower hull)
  NumericMatrix dp(face_filt.nrow(), 2);
  for (int i=0;i<face_filt.nrow();++i){
    IntegerVector tri = face_filt(i,_);
    NumericMatrix tri_pl(3,3);
    for (int t=0;t<3;++t) tri_pl(t, _) = pl(tri[t]-1, _);
    NumericVector dpi = face_dual_uv(tri_pl);
    dp(i,0)=dpi[0]; dp(i,1)=dpi[1];
  }

  // 2D convex hull (base::chull), reverse and close
//   IntegerVector K_core = chull(uv); // K_core degenerates in rect case
  IntegerVector K_core = compute_bd(face_filt, false);
  std::reverse(K_core.begin(), K_core.end());
  IntegerVector K(K_core.size()+1);
  for (int i=0;i<K_core.size();++i) K[i]=K_core[i];
  K[K_core.size()] = K_core[0];

  // big bounding box around dp with padding 1
  NumericVector mindp = NumericVector::create(dp(0,0), dp(0,1));
  NumericVector maxdp = NumericVector::create(dp(0,0), dp(0,1));
  for (int i=1;i<dp.nrow();++i){
    mindp[0]=std::min(mindp[0], dp(i,0));
    mindp[1]=std::min(mindp[1], dp(i,1));
    maxdp[0]=std::max(maxdp[0], dp(i,0));
    maxdp[1]=std::max(maxdp[1], dp(i,1));
  }
  double minx=mindp[0]-1, miny=mindp[1]-1, maxx=maxdp[0]+1, maxy=maxdp[1]+1;
  NumericMatrix box(5,2);
  box(0,0)=minx; box(0,1)=miny;
  box(1,0)=maxx; box(1,1)=miny;
  box(2,0)=maxx; box(2,1)=maxy;
  box(3,0)=minx; box(3,1)=maxy;
  box(4,0)=minx; box(4,1)=miny;

  // vb: boundary ray hits
  NumericMatrix vb((int)K.size()-1, 2);
  for (int i=0;i<K.size()-1;++i){
    int i1 = K[i]-1, i2 = K[i+1]-1;
    double vx = uv(i2,0)-uv(i1,0), vy = uv(i2,1)-uv(i1,1);
    NumericVector vec = NumericVector::create( vy, -vx ); // rotate 90°
    NumericVector mid = NumericVector::create( (uv(i2,0)+uv(i1,0))/2.0,
                                               (uv(i2,1)+uv(i1,1))/2.0 );
    NumericVector ray = NumericVector::create(mid[0], mid[1], vec[0], vec[1]);
    NumericMatrix isects = intersectRayPolygon(ray, box);
    if (isects.nrow() == 0) { vb(i,0)=mid[0]; vb(i,1)=mid[1]; }
    else { vb(i,0)=isects(0,0); vb(i,1)=isects(0,1); }
  }

  // dpe = [dp; vb]
  NumericMatrix dpe(dp.nrow()+vb.nrow(), 2);
  for (int i=0;i<dp.nrow();++i){ dpe(i,0)=dp(i,0); dpe(i,1)=dp(i,1); }
  for (int i=0;i<vb.nrow();++i){ dpe(dp.nrow()+i,0)=vb(i,0); dpe(dp.nrow()+i,1)=vb(i,1); }

  // vvif sparse connectivity
  SEXP vvif = compute_connectivity(face_filt); // dgCMatrix

  // build cells (simple dense access for clarity; can be optimized)
  List cell(n);
  Function asMat("as.matrix");
  NumericMatrix VV = asMat(vvif); // NOTE: for large meshes, replace with sparse access
  for (int i=1;i<=n;++i){
    IntegerVector vri = vr[i-1];
    int pb = NA_INTEGER;
    for (int k=0;k<K.size();++k){ if (K[k]==i){ pb=k+1; break; } } // 1-based

    if (pb != NA_INTEGER) {
      IntegerVector fr((int)vri.size()+1);
      fr[fr.size()-1] = face_filt.nrow() + pb;
      if (pb == 1) fr[0] = face_filt.nrow() + (K.size()-1);
      else         fr[0] = face_filt.nrow() + (pb - 1);
      for (int j=0;j<(int)vri.size()-1;++j) fr[j+1] = (int)VV(i-1, vri[j]-1);
      IntegerVector fr_rev(fr.size());
      for (int t=0;t<fr.size();++t) fr_rev[t]=fr[fr.size()-1-t];
      IntegerVector fr_rev_nz = fr_rev[ fr_rev != 0 ];
      cell[i-1] = fr_rev_nz;
    } else {
      IntegerVector fr((int)vri.size());
      for (int j=0;j<(int)vri.size();++j) fr[j] = (int)VV(i-1, vri[j]-1);
      IntegerVector fr_rev(fr.size());
      for (int t=0;t<fr.size();++t) fr_rev[t]=fr[fr.size()-1-t];
      IntegerVector fr_rev_nz = fr_rev[ fr_rev != 0 ];
      cell[i-1] = fr_rev_nz;
    }
  }

  return List::create(
    _["face"] = face_filt,
    _["uv"]   = uv,
    _["dp"]   = dp,
    _["dpe"]  = dpe,
    _["cell"] = cell,
    _["h"]    = h
  );
}
