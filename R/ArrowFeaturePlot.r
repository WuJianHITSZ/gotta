#' Arrow Plot (Helper Function)
#' 
#' Plots spots at given layouts and overlays arrows indicating trajectory vectors.
#' 
#' @param mesh A SpatialMesh object
#' @param layout.name Character. Name of layout in mesh@layout or mesh@cartogram
#' @param trajectory_name Character. Name of trajectory in mesh@layout
#' @param feature.name Character. Name of feature in mesh@vertex_metadata to use for point color (optional)
#' @param arrow.length Scalar multiplier for arrow length. If NULL, arrows are auto-scaled (default: NULL)
#' @param arrow.size Size of arrow heads (default: 0.3)
#' @param point.size Size of spot points (default: 1)
#' @param point.color Color of spot points. Can be a single color string or a numeric vector of length N (default: "black")
#' @param arrow.color Color of arrows. Can be a single color string or a numeric vector of length N (default: "red")
#' @param title Plot title (default: NULL)
#' @param xlab X-axis label (default: "X")
#' @param ylab Y-axis label (default: "Y")
#' @param coord.fixed Logical. Whether to use fixed aspect ratio (default: TRUE)
#' @param is.show.arrow Logical. Whether to show arrows (default: TRUE)
#' @return ggplot object
#' @noRd
ArrowPlot <- function(mesh,
                     layout.name,
                     trajectory_name,
                     feature.name = NULL,
                     arrow.length = NULL,
                     arrow.size = 0.1,
                     point.size = 4,
                     point.color = "black",
                     arrow.color = "black",
                     title = NULL,
                     xlab = "X",
                     ylab = "Y",
                     coord.fixed = TRUE,
                     is.show.arrow = TRUE,
                     cell.limit = NULL) {
  
  # Validate mesh
  if (!methods::is(mesh, "SpatialMesh")) {
    stop("mesh must be a SpatialMesh object")
  }
  
  # Store original layout.name for axis labels
  layout.name.original <- layout.name
  
  # Handle layout.name: can be a string (lookup in mesh) or a matrix (use directly)
  if (is.character(layout.name) && length(layout.name) == 1) {
    # String input: get layout from mesh using layout.name
    layout <- GetLayout(mesh, layout.name)
  } else if (is.matrix(layout.name) || is.data.frame(layout.name)) {
    # Matrix input: use directly as layout
    layout <- layout.name
    layout.name <- NULL  # Clear layout.name so we use xlab/ylab defaults
  } else {
    stop("layout.name must be a single character string or a matrix/data.frame")
  }
  
  # Handle trajectory_name: can be a string (lookup in mesh) or a matrix (use directly)
  if (is.character(trajectory_name) && length(trajectory_name) == 1) {
    # String input: get trajectory from mesh@layout using trajectory_name
    if (!trajectory_name %in% names(mesh@layout)) {
      stop(sprintf("Trajectory '%s' not found in mesh@layout", trajectory_name))
    }
    trajectory <- mesh@layout[[trajectory_name]]
  } else if (is.matrix(trajectory_name) || is.data.frame(trajectory_name)) {
    # Matrix input: use directly as trajectory
    trajectory <- trajectory_name
  } else {
    stop("trajectory_name must be a single character string or a matrix/data.frame")
  }

  # reverse the trajectory
  trajectory <- -trajectory
  
  # Validate layout and trajectory matrices
  if (!is.matrix(layout) && !is.data.frame(layout)) {
    stop("layout must be a matrix or data.frame")
  }
  if (is.data.frame(layout)) {
    layout <- as.matrix(layout)
  }
  if (ncol(layout) != 2) {
    stop("layout must have 2 columns")
  }
  
  if (!is.matrix(trajectory) && !is.data.frame(trajectory)) {
    stop("trajectory must be a matrix or data.frame")
  }
  if (is.data.frame(trajectory)) {
    trajectory <- as.matrix(trajectory)
  }
  if (ncol(trajectory) != 2) {
    stop("trajectory must have 2 columns")
  }
  
  if (nrow(layout) != nrow(trajectory)) {
    stop(sprintf("layout and trajectory must have the same number of rows (got %d and %d)", 
                 nrow(layout), nrow(trajectory)))
  }
  
  # If feature.name is provided, get point.color from mesh@vertex_metadata
  if (!is.null(feature.name)) {
    if (!is.character(feature.name) || length(feature.name) != 1) {
      stop("feature.name must be a single character string")
    }
    if (!feature.name %in% colnames(mesh@vertex_metadata)) {
      stop(sprintf("Feature '%s' not found in mesh@vertex_metadata", feature.name))
    }
    point.color <- mesh@vertex_metadata[[feature.name]]
  }
  
  # Apply cell.limit filtering if specified
  if (!is.null(cell.limit)) {
    n_keep <- min(nrow(layout), cell.limit)
    layout <- layout[seq_len(n_keep), , drop = FALSE]
    trajectory <- trajectory[seq_len(n_keep), , drop = FALSE]
    if (length(point.color) > 1) {
      point.color <- point.color[seq_len(n_keep)]
    }
    if (length(arrow.color) > 1) {
      arrow.color <- arrow.color[seq_len(n_keep)]
    }
  }
  
  n_spots <- nrow(layout)
  
  # Validate and process point.color
  point.color.is.vector <- length(point.color) > 1
  point.color.is.numeric <- is.numeric(point.color) && point.color.is.vector
  point.color.is.categorical <- point.color.is.vector && (is.factor(point.color) || is.character(point.color))
  
  if (point.color.is.vector) {
    if (length(point.color) != n_spots) {
      stop("point.color vector must have length equal to number of spots (", n_spots, ")")
    }
    # Convert factor to character for processing
    if (is.factor(point.color)) {
      point.color <- as.character(point.color)
    }
  }
  
  # Validate and process arrow.color
  arrow.color.is.vector <- length(arrow.color) > 1
  arrow.color.is.numeric <- is.numeric(arrow.color) && arrow.color.is.vector
  arrow.color.is.categorical <- arrow.color.is.vector && (is.factor(arrow.color) || is.character(arrow.color))
  
  if (arrow.color.is.vector) {
    if (length(arrow.color) != n_spots) {
      stop("arrow.color vector must have length equal to number of spots (", n_spots, ")")
    }
    # Convert factor to character for processing
    if (is.factor(arrow.color)) {
      arrow.color <- as.character(arrow.color)
    }
  }
  
  # Note: For categorical colors, we let ggplot2 use default discrete scales
  # (matching MeshFeaturePlot behavior which doesn't specify a scale for categorical)
  
  # Calculate arrow endpoints
  if (is.null(arrow.length)) {
    # Auto-scale arrow length based on layout range
    coord_range <- max(c(range(layout[, 1], na.rm = TRUE),
                         range(layout[, 2], na.rm = TRUE)))
    trajectory_magnitude <- sqrt(rowSums(trajectory^2, na.rm = TRUE))
    max_trajectory <- max(trajectory_magnitude[trajectory_magnitude > 0], na.rm = TRUE)
    if (max_trajectory > 0) {
      arrow.length <- coord_range * 0.05 / max_trajectory
    } else {
      arrow.length <- 1
    }
  }
  
  # Scale trajectory vectors
  trajectory_scaled <- trajectory * arrow.length
  
  # Create data frame for plotting (swap x-y axes to match MeshFeaturePlot)
  plot_data <- data.frame(
    x = layout[, 2],  # Column 2 becomes x
    y = layout[, 1],  # Column 1 becomes y
    xend = layout[, 2] + trajectory_scaled[, 2],  # Swap trajectory axes too
    yend = layout[, 1] + trajectory_scaled[, 1]
  )
  
  # Add color mappings to plot_data if vector (numeric or categorical)
  if (point.color.is.numeric) {
    plot_data$point_color <- point.color
  } else if (point.color.is.categorical) {
    plot_data$point_color <- point.color
  }
  
  if (arrow.color.is.numeric) {
    plot_data$arrow_color <- arrow.color
  } else if (arrow.color.is.categorical) {
    plot_data$arrow_color <- arrow.color
  }
  
  # Build plot based on which colors are vectors
  # Use color for arrows and fill for points to allow separate scales
  p <- ggplot(plot_data)
  
  # Determine if we have vector colors (numeric or categorical)
  point.is.vector <- point.color.is.numeric || point.color.is.categorical
  arrow.is.vector <- arrow.color.is.numeric || arrow.color.is.categorical
  
  if (arrow.is.vector && point.is.vector) {
    # Both are vectors - use color for arrows, fill for points
    if (is.show.arrow) {
      if (arrow.color.is.numeric) {
        p <- p +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = arrow_color),
                       arrow = arrow(length = unit(arrow.size, "cm"), type = "closed"),
                       linewidth = 0.5) +
          scale_color_viridis_c(option = "turbo")
      } else {
        # Categorical arrow color - use default discrete scale
        p <- p +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = arrow_color),
                       arrow = arrow(length = unit(arrow.size, "cm"), type = "closed"),
                       linewidth = 0.5)
      }
    }
    # Use color aesthetic to match MeshFeaturePlot (not fill)
    p <- p +
      geom_point(aes(x = x, y = y, color = point_color),
                 size = point.size)
    if (point.color.is.numeric) {
      p <- p + scale_color_viridis_c(option = "turbo")
    }
    # For categorical, use default discrete scale (matching MeshFeaturePlot)
  } else if (arrow.is.vector) {
    # Only arrow color is vector
    if (is.show.arrow) {
      if (arrow.color.is.numeric) {
        p <- p +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = arrow_color),
                       arrow = arrow(length = unit(arrow.size, "cm"), type = "closed"),
                       linewidth = 0.5) +
          scale_color_viridis_c(option = "turbo")
      } else {
        # Categorical arrow color - use default discrete scale
        p <- p +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = arrow_color),
                       arrow = arrow(length = unit(arrow.size, "cm"), type = "closed"),
                       linewidth = 0.5)
      }
    }
    p <- p +
      geom_point(aes(x = x, y = y),
                 size = point.size,
                 color = point.color)
  } else if (point.is.vector) {
    # Only point color is vector
    # Use color aesthetic to match MeshFeaturePlot (not fill)
    p <- p +
      geom_point(aes(x = x, y = y, color = point_color),
                 size = point.size)
    if (point.color.is.numeric) {
      p <- p + scale_color_viridis_c(option = "turbo")
    }
    # For categorical, use default discrete scale (matching MeshFeaturePlot)
    if (is.show.arrow) {
      p <- p +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                     color = arrow.color,
                     arrow = arrow(length = unit(arrow.size, "cm"), type = "closed"),
                     linewidth = 0.5)
    }
  } else {
    # Both are single colors
    p <- p +
      geom_point(aes(x = x, y = y),
                 size = point.size,
                 color = point.color)
    if (is.show.arrow) {
      p <- p +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                     color = arrow.color,
                     arrow = arrow(length = unit(arrow.size, "cm"), type = "closed"),
                     linewidth = 0.5)
    }
  }
  
  # Generate axis labels from layout.name (match MeshFeaturePlot style)
  # If layout.name.original is a character, use it; otherwise use xlab/ylab defaults
  if (is.character(layout.name.original) && length(layout.name.original) == 1) {
    x_label <- paste(layout.name.original, "2", sep = "_")
    y_label <- paste(layout.name.original, "1", sep = "_")
  } else {
    x_label <- xlab
    y_label <- ylab
  }
  
  # Add labels, theme, and styling (match MeshFeaturePlot)
  # Set legend name from feature.name if available (matching MeshFeaturePlot's color = features)
  labs_list <- list(
    title = title,
    x = x_label,
    y = y_label
  )
  if (!is.null(feature.name) && point.is.vector) {
    labs_list$color <- feature.name
  }
  p <- p +
    do.call(labs, labs_list) +
    scale_y_reverse() +  # Reverse y-axis to match MeshFeaturePlot
    coord_fixed(ratio = 1) +  # Always use fixed aspect ratio
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  return(p)
}

#' Arrow Feature Plot
#' 
#' S3 generic function that plots spots at given layouts and overlays arrows 
#' indicating trajectory vectors. Wrapper around ArrowPlot that automatically 
#' sets layout and trajectory based on source.assay, target.assay, and shape parameters.
#' 
#' @param object A SpatialMesh or Seurat object
#' @param source.assay Character. Source assay name for the layout (default: "spatial")
#' @param target.assay Character. Target assay name for the feature
#' @param shape Character. Shape type (default: "free_boundary")
#' @param features Character. Feature name in mesh@vertex_metadata to use for point color. If NULL, uses h.<target.assay>.<shape> (default: NULL)
#' @param cell.limit Integer. Maximum number of cells to plot (default: NULL)
#' @param ... Additional arguments passed to ArrowPlot (arrow.length, arrow.size, point.size, arrow.color, title, xlab, ylab, coord.fixed, is.show.arrow)
#' @return ggplot object
#' @export
ArrowFeaturePlot <- function(object, ...) {
  UseMethod("ArrowFeaturePlot")
}

#' Arrow Feature Plot method for SpatialMesh objects
#' 
#' @param object A SpatialMesh object
#' @param source.assay Character. Source assay name for the layout (default: "spatial")
#' @param target.assay Character. Target assay name for the feature
#' @param shape Character. Shape type (default: "free_boundary")
#' @param features Character. Feature name in mesh@vertex_metadata to use for point color. If NULL, uses h.<target.assay>.<shape> (default: NULL)
#' @param alignment.job GOTTJob or AlignmentJob object. If provided, parses source.assay, target.assay, and shape from it (default: NULL)
#' @param cell.limit Integer. Maximum number of cells to plot (default: NULL)
#' @param ... Additional arguments passed to ArrowPlot
#' @return ggplot object
#' @export
ArrowFeaturePlot.SpatialMesh <- function(object,
                                        source.assay = "spatial",
                                        target.assay,
                                        shape = "free_boundary",
                                        features = NULL,
                                        alignment.job = NULL,
                                        cell.limit = NULL,
                                        ...) {
  
  mesh <- object
  
  # Parse alignment.job if provided
  if (!is.null(alignment.job)) {
    # Check if alignment.job is a GOTTJob
    if (inherits(alignment.job, "GOTTJob")) {
      source.assay <- "spatial"
      target.assay <- alignment.job$vertex.job$assay
      shape <- alignment.job$shape
    } 
    # Check if alignment.job is an AlignmentJob
    else if (inherits(alignment.job, "AlignmentJob")) {
      source.assay <- alignment.job$source.assay
      target.assay <- alignment.job$target.assay
      shape <- alignment.job$shape
    } else {
      stop("alignment.job must be a GOTTJob or AlignmentJob object")
    }
  }
  
  # Validate inputs
  if (!is.character(source.assay) || length(source.assay) != 1) {
    stop("source.assay must be a single character string")
  }
  if (!is.character(target.assay) || length(target.assay) != 1) {
    stop("target.assay must be a single character string")
  }
  if (!is.character(shape) || length(shape) != 1) {
    stop("shape must be a single character string")
  }
  if (!is.null(features) && (!is.character(features) || length(features) != 1)) {
    stop("features must be a single character string or NULL")
  }
  
  # Determine layout.name based on source.assay and shape (same logic as FindTrajectory)
  if (source.assay == "spatial" && shape == "free_boundary") {
    layout.name <- "spatial_coords"
  } else if (source.assay == "spatial") {
    layout.name <- paste0("spatial.", shape)
  } else {
    layout.name <- paste0(source.assay, ".", shape)
  }
  
  # Determine trajectory_name
  trajectory_name <- paste0("arrow.", source.assay, "2", target.assay, ".", shape)
  
  # Determine feature.name
  if (is.null(features)) {
    feature.name <- paste0("h.", target.assay, ".", shape)
  } else {
    feature.name <- features
  }
  
  # Call ArrowPlot with mesh, layout.name, trajectory_name, feature.name, and cell.limit
  ArrowPlot(mesh = mesh,
            layout.name = layout.name,
            trajectory_name = trajectory_name,
            feature.name = feature.name,
            cell.limit = cell.limit,
            ...)
}

#' Arrow Feature Plot method for Seurat objects
#' 
#' @param object A Seurat object
#' @param source.assay Character. Source assay name for the layout (default: "spatial")
#' @param target.assay Character. Target assay name for the feature
#' @param shape Character. Shape type (default: "free_boundary")
#' @param features Character. Feature name in mesh@vertex_metadata to use for point color. If NULL, uses h.<target.assay>.<shape> (default: NULL)
#' @param alignment.job GOTTJob or AlignmentJob object. If provided, parses source.assay, target.assay, and shape from it (default: NULL)
#' @param cell.limit Integer. Maximum number of cells to plot (default: NULL, uses ncol(object))
#' @param ... Additional arguments passed to ArrowPlot
#' @return ggplot object
#' @export
ArrowFeaturePlot.Seurat <- function(object,
                                   source.assay = "spatial",
                                   target.assay,
                                   shape = "free_boundary",
                                   features = NULL,
                                   alignment.job = NULL,
                                   cell.limit = NULL,
                                   ...) {
  
  # Extract SpatialMesh from Seurat object
  mesh <- SpatialMesh(object)
  
  # Apply default cell.limit from Seurat object if not specified
  if (is.null(cell.limit)) {
    cell.limit <- ncol(object)
  }
  
  # If features is provided, try to get it from mesh@vertex_metadata first,
  # then from Seurat object via FetchData if not found
  if (!is.null(features)) {
    if (!features %in% colnames(mesh@vertex_metadata)) {
      # Feature not in vertex_metadata, try FetchData from Seurat object
      mesh <- FetchData2Mesh(mesh, object, features)
    }
  }
  
  # Call the SpatialMesh method
  ArrowFeaturePlot.SpatialMesh(mesh,
                               source.assay = source.assay,
                               target.assay = target.assay,
                               shape = shape,
                               features = features,
                               alignment.job = alignment.job,
                               cell.limit = cell.limit,
                               ...)
}

