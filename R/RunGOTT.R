# Constants for shape types
SHAPE_FREE_BOUNDARY <- "free_boundary"
SHAPE_RECT <- "rect"
SHAPE_DISK <- "disk"

#' RunGOTT S3 generic
RunGOTT <- function(object, ...) {
  UseMethod("RunGOTT")
}

#' Default RunGOTT method for SpatialMesh objects
RunGOTT.SpatialMesh <- function(object, shape = SHAPE_FREE_BOUNDARY, aspect.ratio = 1,
                                reduction = NULL, verbose = FALSE, gott.job = NULL,
                                assay = NULL, is.pseudo.initial = FALSE) {

  if (verbose) message("[RunGOTT] Starting RunGOTT...")

  mesh <- object
  face <- mesh@face
  corner <- mesh@corner   # cached for boundary cleanup

  if (!is.null(gott.job)) {
    shape <- gott.job$shape
    aspect.ratio <- gott.job$aspect.ratio
    assay <- gott.job$vertex.job$assay
    reduction <- gott.job$vertex.job$reduction
  } else {
    gott.job <- GOTTJob(assay = assay, reduction = reduction,
                        shape = shape, aspect.ratio = aspect.ratio)
  }

  if (verbose) message(sprintf("[RunGOTT] Parameters: shape=%s, aspect.ratio=%s, assay=%s", shape, aspect.ratio, assay))

  if (!is.null(reduction)) {
    assay <- reduction
  }

  layout.name <- paste0(assay, ".", shape)   # where results will be written

  area_key <- paste0("area.", assay)
  area <- mesh@vertex_metadata[[area_key]]

  spatial_initial <- if (is.pseudo.initial) {
    mesh@layout$spatial_pseudo
  } else {
    mesh@layout$spatial_coords
  }

  if (shape != SHAPE_FREE_BOUNDARY) {
    if (verbose) message("[RunGOTT] Processing shape boundary...")
    result <- process_shape_boundary(mesh, spatial_initial, shape, area,
                                     gott.job, layout.name, face, verbose)
    
  } else {
    if (verbose) message("[RunGOTT] Processing free boundary...")
    result <- process_free_boundary(mesh, spatial_initial, face, area,
                                    gott.job, layout.name, verbose)
  }
  mesh <- result$mesh
  
  mesh@layout[[layout.name]] <- result$uv_vertex
  mesh@cartogram[[layout.name]] <- result$cartogram

  if (verbose) message("[RunGOTT] Updating mesh metadata...")
  mesh <- update_mesh_metadata(mesh, face, uv_vertex = result$uv_vertex,
                               cartogram = result$cartogram, layout.name = layout.name)

  if (verbose) message("[RunGOTT] RunGOTT finished.")
  mesh
}

#' RunGOTT method for Seurat objects
RunGOTT.Seurat <- function(object, ...) {
  mesh <- SpatialMesh(object)
  mesh_updated <- RunGOTT(mesh, ...)
  SpatialMesh(object) <- mesh_updated
  object
}

#' Process shape with boundary constraints
process_shape_boundary <- function(mesh, spatial_initial, shape, area,
                                   gott.job, layout.name, face, verbose = FALSE) {
  spatial_shape <- paste0("spatial.", shape)
  layout_names <- names(mesh@layout)

  if (!spatial_shape %in% layout_names) {
    if (verbose) message("[process_shape_boundary] Initializing spatial coordinates (GOTT)...")
    cartogram_spatial <- GOTT(mesh, spatial_initial, gott.job = gott.job,
                              has.initialized = FALSE, verbose = verbose)
    uv_spatial <- cartogram_spatial$uv
    mesh@layout[[spatial_shape]] <- uv_spatial
  } else {
    if (verbose) message("[process_shape_boundary] Using existing spatial coordinates...")
    uv_spatial <- mesh@layout[[spatial_shape]]
  }

  if (verbose) message("[process_shape_boundary] Computing final cartogram (GOTT)...")
  cartogram <- GOTT(mesh, uv_spatial, population = area, gott.job = gott.job, verbose = verbose)
  uv_vertex <- cartogram$uv

  list(mesh = mesh, uv_vertex = uv_vertex, cartogram = cartogram)
}

#' Process free boundary case (disk construction)
process_free_boundary <- function(mesh, spatial_initial, face, area,
                                  gott.job, layout.name, verbose = FALSE) {
  if (verbose) message("[process_free_boundary] Constructing disk...")
  disk <- construct_sea(spatial_initial, face)
  boundary <- compute_bd(disk$face)

  mesh.disk <- new("SpatialMesh",
                   face = disk$face,
                   boundary = boundary,
                   corner = boundary[1L],
                   layout = list(spatial_coords = disk$vertex))

  n_vertices <- nrow(disk$vertex)
  n_area <- length(area)
  population <- numeric(n_vertices)

  if (n_area < n_vertices) {
    population[seq_len(n_area)] <- area
    mean_area <- mean(area)
    population[(n_area + 1L):n_vertices] <- mean_area
  } else {
    population <- area[seq_len(n_vertices)]
  }

  if (verbose) message("[process_free_boundary] Computing cartogram for free boundary (GOTT)...")
  cartogram <- GOTT(mesh.disk, disk$vertex, population = population,
                    shape = SHAPE_DISK, aspect.ratio = 1, has.initialized = TRUE, verbose = verbose)

  n_real <- min(n_area, n_vertices)
  uv_vertex <- cartogram$uv[seq_len(n_real), , drop = FALSE]

  cartogram$uv <- uv_vertex
  cartogram$cell <- cartogram$cell[seq_len(n_real)]
  cartogram$h <- cartogram$h[seq_len(n_real)]  # trim sea vertices

  mesh@layout$spatial.free_boundary <- disk$vertex[seq_len(n_real), , drop = FALSE]

  list(mesh = mesh, uv_vertex = uv_vertex, cartogram = cartogram)
}

#' Update mesh metadata with computed values
update_mesh_metadata <- function(mesh, face, uv_vertex, cartogram, layout.name) {
  face_valid <- face
  if (!is.integer(face_valid)) {
    storage.mode(face_valid) <- "integer"
  }
  if (nrow(face_valid) > 0L) {
    flag <- rowSums(face_valid > nrow(uv_vertex)) == 0L
    face_valid <- face_valid[flag, , drop = FALSE]
  }
  vertex_area_calc <- vertex_area(face_valid, uv_vertex) / 3   # each face contributes to 3 vertices
  face_area_calc <- face_area(face_valid, uv_vertex)

  mesh@vertex_metadata[[paste0("h.", layout.name)]] <- cartogram$h
  mesh@vertex_metadata[[paste0("area.", layout.name)]] <- vertex_area_calc
  mesh@face_metadata[[paste0("area.", layout.name)]] <- face_area_calc

  mesh
}

#' Core GOTT computation function
GOTT <- function(mesh, uv_initial, population = NULL, gott.job = NULL,
                 has.initialized = TRUE, aspect.ratio = 1,
                 shape = SHAPE_FREE_BOUNDARY, verbose = FALSE) {

  if (!is.null(gott.job)) {
    shape <- gott.job$shape
    aspect.ratio <- gott.job$aspect.ratio
  }

  boundary_index <- mesh@boundary
  corner <- mesh@corner
  face <- mesh@face

  if (!has.initialized) {
    if (verbose) message("[GOTT] Initializing UV coordinates...")
    uv_initial <- initialize_uv_coordinates(face, uv_initial, corner, shape, aspect.ratio)
  }

  area_initial <- vertex_area(face, uv_initial) / 3
  if (is.null(population)) {
    population <- area_initial
  } else {
    sum_area <- sum(area_initial)
    sum_pop <- sum(population)
    if (sum_pop > 0) {
      population <- population * (sum_area / sum_pop)
    }
  }

  boundary <- uv_initial[boundary_index, , drop = FALSE]
  if (verbose) message("[GOTT] Running discrete optimal transport...")
  cartogram <- discrete_optimal_transport(boundary, face, uv_initial,
                                          function(coord) 1, population)
  uv_vertex <- compute_centroid(boundary, cartogram)

  if (verbose) message("[GOTT] Cleaning up boundary...")
  uv_vertex <- cleanup_boundary(uv_vertex, boundary_index, corner, shape,
                                aspect.ratio, uv_initial)

  cartogram$uv <- uv_vertex
  cartogram$face <- NULL
  cartogram$dp <- NULL

  cartogram
}

#' Initialize UV coordinates based on shape
initialize_uv_coordinates <- function(face, uv_initial, corner, shape, aspect.ratio) {
  if (shape == SHAPE_RECT) {
    uv_initial <- rect_harmonic_map(face, uv_initial, corner)
  } else {
    uv_initial <- disk_harmonic_map(face, uv_initial, corner[1L])
  }
  uv_initial[, 2L] <- uv_initial[, 2L] * aspect.ratio
  uv_initial
}

#' Clean up boundary vertices based on shape type
cleanup_boundary <- function(uv_vertex, boundary_index, corner,
                             shape, aspect.ratio, uv_initial) {
  if (shape == SHAPE_RECT) {
    uv_vertex <- cleanup_rect_boundary(uv_vertex, boundary_index, corner)
  } else {
    uv_vertex <- cleanup_disk_boundary(uv_vertex, boundary_index, aspect.ratio)
  }
  uv_vertex
}

#' Clean up rectangular boundary
cleanup_rect_boundary <- function(uv_vertex, boundary_index, corner) {
  corner_indices <- find_corner_indices(boundary_index, corner)
  i1 <- corner_indices[1L]
  n_boundary <- length(boundary_index)
  boundary_index <- c(boundary_index[i1:n_boundary], boundary_index[1L:i1])
  corner_indices <- find_corner_indices(boundary_index, corner)
  i1 <- corner_indices[1L]
  i2 <- corner_indices[2L]
  i3 <- corner_indices[3L]
  i4 <- corner_indices[4L]
  i0 <- corner_indices[5L]

  L <- c(min(uv_vertex[, 1L]), min(uv_vertex[, 2L]))
  U <- c(max(uv_vertex[, 1L]), max(uv_vertex[, 2L]))

  uv_vertex[boundary_index[i1:i2], 2L] <- L[2L]
  uv_vertex[boundary_index[i2:i3], 1L] <- U[1L]
  uv_vertex[boundary_index[i3:i4], 2L] <- U[2L]
  uv_vertex[boundary_index[i4:i0], 1L] <- L[1L]

  uv_vertex
}

#' Find corner indices in boundary
find_corner_indices <- function(boundary_index, corner) {
  i1_first <- which(boundary_index == corner[1L])[1L]
  i1_last <- tail(which(boundary_index == corner[1L]), 1L)
  i2 <- which(boundary_index == corner[2L])[1L]
  i3 <- which(boundary_index == corner[3L])[1L]
  i4 <- which(boundary_index == corner[4L])[1L]

  c(i1_first, i2, i3, i4, i1_last)
}

#' Clean up disk boundary
cleanup_disk_boundary <- function(uv_vertex, boundary_index, aspect.ratio) {
  ar <- if (is.null(aspect.ratio)) 1 else aspect.ratio

  uv_boundary <- uv_vertex[boundary_index, , drop = FALSE]
  uv_boundary[, 2L] <- uv_boundary[, 2L] / ar

  cbl <- sqrt(uv_boundary[, 1L]^2 + uv_boundary[, 2L]^2)
  cbl[cbl == 0] <- 1
  uv_boundary <- uv_boundary / cbind(cbl, cbl)
  uv_boundary[, 2L] <- uv_boundary[, 2L] * ar

  uv_vertex[boundary_index, ] <- uv_boundary
  uv_vertex
}
