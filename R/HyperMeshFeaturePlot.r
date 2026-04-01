## S3 hyper mesh feature plotting (v1.6.7)

# Generics -----------------------------------------------------------------
HyperMeshFeaturePlot <- function(object, ...) UseMethod("HyperMeshFeaturePlot")
HyperSurfFeaturePlot <- function(object, ...) UseMethod("HyperSurfFeaturePlot")

# SpatialMesh methods ------------------------------------------------------
HyperMeshFeaturePlot.SpatialMesh <- function(object, features, layout.name = "spatial_coords",
                                             title = NULL, pt.size = 2,
                                             feature.values = NULL, cell.limit = NULL) {
  layout.name <- paste0(layout.name, ".hyperview")
  layout <- .get_mesh_layout(object, layout.name, cell.limit)
  color_values <- .resolve_vertex_feature(object, features, values_override = feature.values, limit = nrow(layout))

  face <- .filter_faces(.get_mesh_face(object, use_quad = FALSE), nrow(layout))
  face <- face - 1L

  color_levels <- unique(color_values)
  palette_fn <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
  color_map <- palette_fn(length(color_levels))
  colors <- color_map[match(color_values, color_levels)]

  mesh_plot <- plot_ly()
  mesh_plot <- add_trace(
    mesh_plot,
    type = "mesh3d",
    x = layout[, 1],
    y = layout[, 2],
    z = layout[, 3],
    i = face[, 1],
    j = face[, 2],
    k = face[, 3],
    color = I("rgba(255, 0, 0, 0.7)"),
    opacity = 0.7,
    lighting = list(ambient = 0.8)
  )

  mesh_plot <- add_trace(
    mesh_plot,
    type = "scatter3d",
    mode = "markers",
    x = layout[, 1],
    y = layout[, 2],
    z = layout[, 3],
    marker = list(
      size = pt.size,
      color = colors,
      colorscale = "Viridis",
      opacity = 0.8
    )
  )

  if (!is.null(title)) {
    mesh_plot <- layout(mesh_plot, title = title)
  }

  mesh_plot
}

HyperSurfFeaturePlot.SpatialMesh <- function(object, features, layout.name = "spatial_coords",
                                             title = NULL, pt.size = 2,
                                             feature.values = NULL, cell.limit = NULL) {
  layout.name <- paste0(layout.name, ".hyperview")
  layout <- .get_mesh_layout(object, layout.name, cell.limit)
  color <- .resolve_vertex_feature(object, features, values_override = feature.values, limit = nrow(layout))

  face <- .filter_faces(.get_mesh_face(object, use_quad = FALSE), nrow(layout))
  face <- face - 1L

  surface_plot <- plot_ly()
  surface_plot <- add_mesh(
    surface_plot,
    x = layout[, 1],
    y = layout[, 2],
    z = layout[, 3],
    i = face[, 1],
    j = face[, 2],
    k = face[, 3],
    intensity = color,
    colorscale = "Viridis",
    showscale = TRUE
  )

  if (is.null(title)) {
    surface_plot <- layout(surface_plot, title = "Gaussian Curvature")
  } else {
    surface_plot <- layout(surface_plot, title = title)
  }

  surface_plot
}

# Seurat methods -----------------------------------------------------------
HyperMeshFeaturePlot.Seurat <- function(object, features, layout.name = "spatial_coords",
                                        title = NULL, pt.size = 2, ...) {
  mesh <- SpatialMesh(object)
  values_override <- NULL
  if (!(features %in% colnames(mesh@vertex_metadata)) && !(features %in% colnames(mesh@face_metadata))) {
    fetched <- FetchData(object, vars = features)
    if (ncol(fetched) != 1) {
      stop(sprintf("Expected one column when fetching '%s' from Seurat.", features))
    }
    values_override <- fetched[, 1, drop = TRUE]
  }
  HyperMeshFeaturePlot(mesh, features = features, layout.name = layout.name, title = title, pt.size = pt.size,
                       feature.values = values_override, cell.limit = ncol(object), ...)
}

HyperSurfFeaturePlot.Seurat <- function(object, features, layout.name = "spatial_coords",
                                        title = NULL, pt.size = 2, ...) {
  mesh <- SpatialMesh(object)
  values_override <- NULL
  if (!(features %in% colnames(mesh@vertex_metadata)) && !(features %in% colnames(mesh@face_metadata))) {
    fetched <- FetchData(object, vars = features)
    if (ncol(fetched) != 1) {
      stop(sprintf("Expected one column when fetching '%s' from Seurat.", features))
    }
    values_override <- fetched[, 1, drop = TRUE]
  }
  HyperSurfFeaturePlot(mesh, features = features, layout.name = layout.name, title = title, pt.size = pt.size,
                       feature.values = values_override, cell.limit = ncol(object), ...)
}
