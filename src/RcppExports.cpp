#include <Rcpp.h>
#include <R_ext/Rdynload.h>

using namespace Rcpp;

// Forward declarations
NumericVector calculate_gradient(NumericMatrix cp, List pd, Function sigma, bool more_accurate = false);
SEXP calculate_hessian(NumericMatrix cp, List pd, Function sigma);
IntegerVector compute_bd(const IntegerMatrix& face_in, bool efficient = false);
NumericMatrix compute_centroid(const NumericMatrix& cp, const List& pd);
List construct_sea(const NumericMatrix& vertex, const IntegerMatrix& face);
NumericVector discrete_gaussian_curvature(const NumericMatrix& vertex, const IntegerMatrix& face);
NumericMatrix disk_harmonic_map(const IntegerMatrix& face, const NumericMatrix& vertex, int fixed_point, std::string method = "Polthier");
NumericVector face_area(const IntegerMatrix& face, const NumericMatrix& vertex);
List polybool(SEXP pa, SEXP pb, std::string op, Rcpp::Nullable<LogicalVector> ha = R_NilValue, Rcpp::Nullable<LogicalVector> hb = R_NilValue, Rcpp::Nullable<NumericVector> ug_in = R_NilValue);
List power_diagram(IntegerMatrix face, NumericMatrix uv, Nullable<NumericVector> h_in = R_NilValue, Nullable<NumericVector> dh_in = R_NilValue);
NumericMatrix rect_harmonic_map(const IntegerMatrix& face, const NumericMatrix& vertex, SEXP corner, std::string method = "Polthier");
Rcpp::List si_create_linear(const Rcpp::NumericMatrix& xy, const Rcpp::NumericVector& z, const Rcpp::IntegerMatrix& tri, const double eps = 1e-12);
Rcpp::NumericVector si_predict_linear(const Rcpp::List& model, SEXP query, double inside_tol = 1e-12);
NumericVector vertex_area(const IntegerMatrix& face, const NumericMatrix& vertex, std::string type = "one_ring");

RcppExport SEXP _gotta_calculate_gradient(SEXP cpSEXP, SEXP pdSEXP, SEXP sigmaSEXP, SEXP more_accurateSEXP) {
    BEGIN_RCPP
    NumericMatrix cp(cpSEXP);
    List pd(pdSEXP);
    Function sigma(sigmaSEXP);
    bool more_accurate = as<bool>(more_accurateSEXP);
    return wrap(calculate_gradient(cp, pd, sigma, more_accurate));
    END_RCPP
}

RcppExport SEXP _gotta_calculate_hessian(SEXP cpSEXP, SEXP pdSEXP, SEXP sigmaSEXP) {
    BEGIN_RCPP
    NumericMatrix cp(cpSEXP);
    List pd(pdSEXP);
    Function sigma(sigmaSEXP);
    return calculate_hessian(cp, pd, sigma);
    END_RCPP
}

RcppExport SEXP _gotta_compute_bd(SEXP face_inSEXP, SEXP efficientSEXP) {
    BEGIN_RCPP
    IntegerMatrix face_in(face_inSEXP);
    bool efficient = as<bool>(efficientSEXP);
    return wrap(compute_bd(face_in, efficient));
    END_RCPP
}

RcppExport SEXP _gotta_compute_centroid(SEXP cpSEXP, SEXP pdSEXP) {
    BEGIN_RCPP
    NumericMatrix cp(cpSEXP);
    List pd(pdSEXP);
    return wrap(compute_centroid(cp, pd));
    END_RCPP
}

RcppExport SEXP _gotta_construct_sea(SEXP vertexSEXP, SEXP faceSEXP) {
    BEGIN_RCPP
    NumericMatrix vertex(vertexSEXP);
    IntegerMatrix face(faceSEXP);
    return wrap(construct_sea(vertex, face));
    END_RCPP
}

RcppExport SEXP _gotta_discrete_gaussian_curvature(SEXP vertexSEXP, SEXP faceSEXP) {
    BEGIN_RCPP
    NumericMatrix vertex(vertexSEXP);
    IntegerMatrix face(faceSEXP);
    return wrap(discrete_gaussian_curvature(vertex, face));
    END_RCPP
}

RcppExport SEXP _gotta_disk_harmonic_map(SEXP faceSEXP, SEXP vertexSEXP, SEXP fixed_pointSEXP, SEXP methodSEXP) {
    BEGIN_RCPP
    IntegerMatrix face(faceSEXP);
    NumericMatrix vertex(vertexSEXP);
    int fixed_point = as<int>(fixed_pointSEXP);
    std::string method = as<std::string>(methodSEXP);
    return wrap(disk_harmonic_map(face, vertex, fixed_point, method));
    END_RCPP
}

RcppExport SEXP _gotta_face_area(SEXP faceSEXP, SEXP vertexSEXP) {
    BEGIN_RCPP
    IntegerMatrix face(faceSEXP);
    NumericMatrix vertex(vertexSEXP);
    return wrap(face_area(face, vertex));
    END_RCPP
}

RcppExport SEXP _gotta_polybool(SEXP paSEXP, SEXP pbSEXP, SEXP opSEXP, SEXP haSEXP, SEXP hbSEXP, SEXP ug_inSEXP) {
    BEGIN_RCPP
    std::string op = as<std::string>(opSEXP);
    Nullable<LogicalVector> ha(haSEXP);
    Nullable<LogicalVector> hb(hbSEXP);
    Nullable<NumericVector> ug_in(ug_inSEXP);
    return wrap(polybool(paSEXP, pbSEXP, op, ha, hb, ug_in));
    END_RCPP
}

RcppExport SEXP _gotta_power_diagram(SEXP faceSEXP, SEXP uvSEXP, SEXP h_inSEXP, SEXP dh_inSEXP) {
    BEGIN_RCPP
    IntegerMatrix face(faceSEXP);
    NumericMatrix uv(uvSEXP);
    Nullable<NumericVector> h_in(h_inSEXP);
    Nullable<NumericVector> dh_in(dh_inSEXP);
    return wrap(power_diagram(face, uv, h_in, dh_in));
    END_RCPP
}

RcppExport SEXP _gotta_rect_harmonic_map(SEXP faceSEXP, SEXP vertexSEXP, SEXP cornerSEXP, SEXP methodSEXP) {
    BEGIN_RCPP
    IntegerMatrix face(faceSEXP);
    NumericMatrix vertex(vertexSEXP);
    std::string method = as<std::string>(methodSEXP);
    return wrap(rect_harmonic_map(face, vertex, cornerSEXP, method));
    END_RCPP
}

RcppExport SEXP _gotta_si_create_linear(SEXP xySEXP, SEXP zSEXP, SEXP triSEXP, SEXP epsSEXP) {
    BEGIN_RCPP
    NumericMatrix xy(xySEXP);
    NumericVector z(zSEXP);
    IntegerMatrix tri(triSEXP);
    double eps = as<double>(epsSEXP);
    return wrap(si_create_linear(xy, z, tri, eps));
    END_RCPP
}

RcppExport SEXP _gotta_si_predict_linear(SEXP modelSEXP, SEXP querySEXP, SEXP inside_tolSEXP) {
    BEGIN_RCPP
    List model(modelSEXP);
    double inside_tol = as<double>(inside_tolSEXP);
    return wrap(si_predict_linear(model, querySEXP, inside_tol));
    END_RCPP
}

RcppExport SEXP _gotta_vertex_area(SEXP faceSEXP, SEXP vertexSEXP, SEXP typeSEXP) {
    BEGIN_RCPP
    IntegerMatrix face(faceSEXP);
    NumericMatrix vertex(vertexSEXP);
    std::string type = as<std::string>(typeSEXP);
    return wrap(vertex_area(face, vertex, type));
    END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gotta_calculate_gradient", (DL_FUNC) &_gotta_calculate_gradient, 4},
    {"_gotta_calculate_hessian", (DL_FUNC) &_gotta_calculate_hessian, 3},
    {"_gotta_compute_bd", (DL_FUNC) &_gotta_compute_bd, 2},
    {"_gotta_compute_centroid", (DL_FUNC) &_gotta_compute_centroid, 2},
    {"_gotta_construct_sea", (DL_FUNC) &_gotta_construct_sea, 2},
    {"_gotta_discrete_gaussian_curvature", (DL_FUNC) &_gotta_discrete_gaussian_curvature, 2},
    {"_gotta_disk_harmonic_map", (DL_FUNC) &_gotta_disk_harmonic_map, 4},
    {"_gotta_face_area", (DL_FUNC) &_gotta_face_area, 2},
    {"_gotta_polybool", (DL_FUNC) &_gotta_polybool, 6},
    {"_gotta_power_diagram", (DL_FUNC) &_gotta_power_diagram, 4},
    {"_gotta_rect_harmonic_map", (DL_FUNC) &_gotta_rect_harmonic_map, 4},
    {"_gotta_si_create_linear", (DL_FUNC) &_gotta_si_create_linear, 4},
    {"_gotta_si_predict_linear", (DL_FUNC) &_gotta_si_predict_linear, 3},
    {"_gotta_vertex_area", (DL_FUNC) &_gotta_vertex_area, 3},
    {NULL, NULL, 0}
};

extern "C" void R_init_gotta(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
