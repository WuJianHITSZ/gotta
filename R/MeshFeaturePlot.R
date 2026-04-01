## S3 mesh feature plotting (v1.6.7)

# Internal helpers ---------------------------------------------------------
.get_mesh_layout <- function(mesh, layout.name, cell.limit = NULL) {
  layout <- mesh@layout[[layout.name]]
  if (is.null(layout)) {
    stop(sprintf("Layout '%s' not found in mesh@layout.", layout.name))
  }
  if (!is.null(cell.limit)) {
    keep <- seq_len(min(nrow(layout), cell.limit))
    layout <- layout[keep, , drop = FALSE]
  }
  layout
}

.get_mesh_face <- function(mesh, use_quad = FALSE) {
  if (!methods::is(mesh, "SpatialMesh")) {
    stop("Mesh face retrieval requires a SpatialMesh object.")
  }
  if (isTRUE(use_quad)) {
    face <- mesh@quad
  } else {
    face <- mesh@face
  }
  if (is.null(face)) {
    stop("Mesh face information is missing.")
  }
  face
}

.filter_faces <- function(face, max_index) {
  face[rowSums(face > max_index) == 0, , drop = FALSE]
}

.resolve_vertex_feature <- function(mesh, features, values_override = NULL, limit = NULL) {
  # Prefer explicit override, then cached vertex metadata
  if (!is.null(values_override)) {
    values <- values_override
  } else if (features %in% colnames(mesh@vertex_metadata)) {
    values <- mesh@vertex_metadata[[features]]
  } else {
    stop(sprintf("Feature '%s' not found in mesh vertex_metadata.", features))
  }
  values <- as.vector(values)
  if (!is.null(limit)) {
    values <- values[seq_len(min(length(values), limit))]
  }
  values
}

.resolve_face_feature <- function(mesh, features, mask = NULL) {
  # Pull face-level feature with optional row mask
  if (!(features %in% colnames(mesh@face_metadata))) {
    stop(sprintf("Feature '%s' not found in mesh face_metadata.", features))
  }
  values <- mesh@face_metadata[[features]]
  if (!is.null(mask)) {
    values <- values[mask]
  }
  values
}

.build_polygon_df <- function(face, fill = NULL) {
  n_face <- nrow(face)
  group <- rep(seq_len(n_face), each = ncol(face))
  df <- data.frame(
    vertex = as.vector(t(face)),
    group = group
  )
  if (!is.null(fill)) {
    df$fill <- rep(fill, each = ncol(face))
  }
  df
}

# Generics -----------------------------------------------------------------
MeshFeaturePlot <- function(object, ...) UseMethod("MeshFeaturePlot")
SurfFeaturePlot <- function(object, ...) UseMethod("SurfFeaturePlot")

MeshFeaturePlot.SpatialMesh <- function(object, features, layout.name = "spatial_coords", title = NULL,
                                        pt.size = 2, show.quad = FALSE, shape = NULL,
                                        feature.values = NULL, cell.limit = NULL) {
  if (!is.null(shape)) {
    layout.name <- paste0(layout.name, ".", shape)
  }

  layout <- .get_mesh_layout(object, layout.name, cell.limit)
  color <- .resolve_vertex_feature(object, features, values_override = feature.values, limit = nrow(layout))

  face <- .filter_faces(.get_mesh_face(object, use_quad = show.quad), nrow(layout))

  plot <- triplot(face = face, layout = layout) +
    geom_point(aes(x = layout[, 2], y = layout[, 1], color = color), size = pt.size) +
    labs(
      title = title,
      x = paste(layout.name, "2", sep = "_"),
      y = paste(layout.name, "1", sep = "_"),
      color = features
    )

  if (is.numeric(color[1])) {
    plot <- plot + scale_color_viridis_c(option = "turbo")
  }

  plot
}

SurfFeaturePlot.SpatialMesh <- function(object, features = NULL,
                                        layout.name = "spatial_coords",
                                        title = NULL,
                                        show.quad = FALSE,
                                        shape = "free_boundary",
                                        assay = NULL,
                                        feature.values = NULL,
                                        cell.limit = NULL) {

  if (!is.null(assay)) {
    layout.name <- paste0(assay, ".", shape)
  }
  if (is.null(features)) {
    features <- paste0("area.", layout.name)
  }

  layout <- .get_mesh_layout(object, layout.name, cell.limit)

  face_original <- .get_mesh_face(object, use_quad = show.quad)
  flag <- rowSums(face_original > nrow(layout)) == 0
  face <- face_original[flag, , drop = FALSE]

  if ((features %in% colnames(object@face_metadata)) && !isTRUE(show.quad)) {
    fill <- .resolve_face_feature(object, features, mask = flag)
  } else {
    color <- .resolve_vertex_feature(object, features, values_override = feature.values, limit = nrow(layout))
    fill <- rowMeans(matrix(color[face], nrow = nrow(face), ncol = ncol(face)))
  }

  plot <- trisurf(face = face, layout = layout, fill = fill)

  plot +
    labs(
      title = title,
      x = paste(layout.name, "2", sep = "_"),
      y = paste(layout.name, "1", sep = "_"),
      fill = features
    )
}

# Seurat methods -----------------------------------------------------------
MeshFeaturePlot.Seurat <- function(object, features, layout.name = "spatial_coords", title = NULL,
                                   pt.size = 2, show.quad = FALSE, shape = NULL, ...) {
  mesh <- SpatialMesh(object)
  values_override <- NULL
  if (!(features %in% colnames(mesh@vertex_metadata)) && !(features %in% colnames(mesh@face_metadata))) {
    fetched <- FetchData(object, vars = features)
    if (ncol(fetched) != 1) {
      stop(sprintf("Expected one column when fetching '%s' from Seurat.", features))
    }
    values_override <- fetched[, 1, drop = TRUE]
  }
  MeshFeaturePlot(mesh, features = features, layout.name = layout.name, title = title,
                  pt.size = pt.size, show.quad = show.quad, shape = shape,
                  feature.values = values_override, cell.limit = ncol(object), ...)
}

SurfFeaturePlot.Seurat <- function(object, features = NULL,
                                   layout.name = "spatial_coords",
                                   title = NULL,
                                   show.quad = FALSE,
                                   shape = "free_boundary",
                                   assay = NULL, ...) {
  mesh <- SpatialMesh(object)
  values_override <- NULL
  feature_name <- if (is.null(features) && !is.null(assay)) paste0("area.", paste0(assay, ".", shape)) else features
  feature_name <- if (is.null(feature_name)) paste0("area.", layout.name) else feature_name

  if (!(feature_name %in% colnames(mesh@vertex_metadata)) && !(feature_name %in% colnames(mesh@face_metadata))) {
    fetched <- FetchData(object, vars = feature_name)
    if (ncol(fetched) != 1) {
      stop(sprintf("Expected one column when fetching '%s' from Seurat.", feature_name))
    }
    values_override <- fetched[, 1, drop = TRUE]
  }

  SurfFeaturePlot(mesh, features = feature_name, layout.name = layout.name, title = title,
                  show.quad = show.quad,
                  shape = shape, assay = assay, feature.values = values_override, cell.limit = ncol(object), ...)
}

# Shared plotting helpers --------------------------------------------------
triplot <- function(face, layout) {
  data_frame_face <- .build_polygon_df(face)

  ggplot() +
    geom_polygon(data = data_frame_face, aes(x = layout[vertex, 2], y = layout[vertex, 1], group = group), fill = NA, color = "black") +
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

trisurf <- function(face, layout, fill) {
  data_frame_face <- .build_polygon_df(face, fill = fill)
  data_frame_face <- subset(data_frame_face, vertex != 0)

  plot <- ggplot() +
    geom_polygon(data = data_frame_face, aes(x = layout[vertex, 2], y = layout[vertex, 1], group = group, fill = fill), color = "black") +
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

  if (is.numeric(fill[1])) {
    plot <- plot + scale_fill_viridis_c(option = "turbo")
  }

  plot
}
