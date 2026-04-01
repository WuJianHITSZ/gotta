# Generated manually to provide the package-side .Call wrappers needed for
# building this project without running compileAttributes().

calculate_gradient <- function(cp, pd, sigma, more_accurate = FALSE) {
    .Call(`_gotta_calculate_gradient`, cp, pd, sigma, more_accurate)
}

calculate_hessian <- function(cp, pd, sigma) {
    .Call(`_gotta_calculate_hessian`, cp, pd, sigma)
}

compute_bd <- function(face_in, efficient = FALSE) {
    .Call(`_gotta_compute_bd`, face_in, efficient)
}

compute_centroid <- function(cp, pd) {
    .Call(`_gotta_compute_centroid`, cp, pd)
}

construct_sea <- function(vertex, face) {
    .Call(`_gotta_construct_sea`, vertex, face)
}

discrete_gaussian_curvature <- function(vertex, face) {
    .Call(`_gotta_discrete_gaussian_curvature`, vertex, face)
}

disk_harmonic_map <- function(face, vertex, fixed_point, method = "Polthier") {
    .Call(`_gotta_disk_harmonic_map`, face, vertex, fixed_point, method)
}

face_area <- function(face, vertex) {
    .Call(`_gotta_face_area`, face, vertex)
}

polybool <- function(pa, pb, op, ha = NULL, hb = NULL, ug_in = NULL) {
    .Call(`_gotta_polybool`, pa, pb, op, ha, hb, ug_in)
}

power_diagram <- function(face, uv, h_in = NULL, dh_in = NULL) {
    .Call(`_gotta_power_diagram`, face, uv, h_in, dh_in)
}

rect_harmonic_map <- function(face, vertex, corner, method = "Polthier") {
    .Call(`_gotta_rect_harmonic_map`, face, vertex, corner, method)
}

si_create_linear <- function(xy, z, tri, eps = 1e-12) {
    .Call(`_gotta_si_create_linear`, xy, z, tri, eps)
}

si_predict_linear <- function(model, query, inside_tol = 1e-12) {
    .Call(`_gotta_si_predict_linear`, model, query, inside_tol)
}

vertex_area <- function(face, vertex, type = "one_ring") {
    .Call(`_gotta_vertex_area`, face, vertex, type)
}
