ComputeCellArea <- function(object, vertex.job = NULL, verbose = FALSE,
                         reduction = NULL, n.components = NULL, assay = DefaultAssay(object), features = NULL){
  if (verbose) message("Start")

  if (!is.null(vertex.job)) {
    assay <- vertex.job$assay %||% assay
    features <- vertex.job$features
    reduction <- vertex.job$reduction
    n.components <- vertex.job$n.components
  }
  if (!is.null(reduction)) {
    assay <- reduction
  }

  mesh <- SpatialMesh(object)
  face <- mesh@face
  spatial_coords <- mesh@layout$spatial_coords

  if (verbose) {
    message(sprintf("Building vertex coordinates (assay=%s, reduction=%s)", assay, reduction %||% "NA"))
  }
  vertex <- VertexAccessor(object, assay, features, reduction, n.components)
  if (verbose) {
    message(sprintf("Vertex dims after access: %s cells x %s components", nrow(vertex), ncol(vertex)))
  }
  if (nrow(spatial_coords) > nrow(vertex)) {
    n_missing <- nrow(spatial_coords) - nrow(vertex)
    if (verbose) message(sprintf("Filling holes in vertex coordinates (adding %s spots)", n_missing))
    class(vertex) <- "vertex"
    vertex <- FillHoles(vertex, face = face, array_coords = spatial_coords)
    if (verbose) {
      message(sprintf("Vertex dims after fill: %s cells x %s components", nrow(vertex), ncol(vertex)))
    }
  }

  if (methods::is(object, "Seurat") &&
      "seurat_clusters" %in% colnames(object@meta.data) &&
      !("seurat_clusters" %in% colnames(mesh@vertex_metadata))) {
    if (verbose) message("Propagating seurat_clusters to mesh vertex_metadata")
    n_cell <- ncol(object)
    if (nrow(mesh@vertex_metadata) == 0L) {
      mesh@vertex_metadata <- data.frame(matrix(NA, nrow = nrow(spatial_coords), ncol = 0))
    }
    clusters <- object@meta.data[colnames(object), "seurat_clusters", drop = TRUE]
    values <- rep(NA, nrow(mesh@vertex_metadata))
    idx <- seq_len(min(n_cell, length(values)))
    values[idx] <- clusters[idx]
    mesh@vertex_metadata[[1]] <- factor(values)
    colnames(mesh@vertex_metadata)[1] <- "seurat_clusters"
  }

  if (verbose) message("Computing areas and curvature")
  mesh@vertex_metadata[[paste0("area.", assay)]] <- vertex_area(face = face, vertex = vertex) / 3
  mesh@vertex_metadata[[paste0("curvature.", assay)]] <- discrete_gaussian_curvature(vertex, face)
  mesh@face_metadata[[paste0("area.", assay)]] <- face_area(face = face, vertex = vertex)

  if (verbose) message("Finished")
  SpatialMesh(object) <- mesh
  object
}
