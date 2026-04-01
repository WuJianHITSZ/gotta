# -------------------------------------------------------------------------
# 1. Core Plotting Logic (Helper Function)
# -------------------------------------------------------------------------
#' Internal helper function for plotting power diagram cells
#' @param layout_list List of matrices, each containing polygon coordinates
#' @param fill_values Vector of fill values (one per cell)
#' @param border_color Color for polygon borders
#' @return ggplot object
.CartFeaturePlot_Core <- function(layout_list, 
                                   fill_values, 
                                   border_color = "black") {
  
  n_cells <- length(layout_list)
  if (n_cells == 0) stop("Layout list is empty. Cannot plot.")
  
  # --- Pre-compute row counts (single pass) ---
  row_counts <- vapply(layout_list, nrow, integer(1))
  total_rows <- sum(row_counts)
  
  if (total_rows == 0) stop("All cells have empty coordinates. Cannot plot.")
  
  # --- Pre-allocate matrices for memory efficiency ---
  # Instead of do.call(rbind), pre-allocate and fill
  layout_mat <- matrix(0, nrow = total_rows, ncol = 2)
  group_vec <- integer(total_rows)
  
  # --- Fill matrices in single pass ---
  row_idx <- 1L
  for (i in seq_len(n_cells)) {
    n_rows <- row_counts[i]
    if (n_rows > 0) {
      end_idx <- row_idx + n_rows - 1L
      layout_mat[row_idx:end_idx, ] <- layout_list[[i]]
      group_vec[row_idx:end_idx] <- i
      row_idx <- end_idx + 1L
    }
  }
  
  # --- Build data frame efficiently ---
  poly_data <- data.frame(
    x = layout_mat[, 2],  # V2 -> x
    y = layout_mat[, 1],  # V1 -> y
    group = group_vec,
    fill = fill_values[group_vec]
  )
  
  # --- Filter NA values before plotting (reduces memory) ---
  valid_rows <- !is.na(poly_data$fill)
  if (!all(valid_rows)) {
    poly_data <- poly_data[valid_rows, , drop = FALSE]
  }
  
  if (nrow(poly_data) == 0) {
    warning("No valid data points after filtering. Returning empty plot.")
    return(ggplot() + theme_void())
  }
  
  # --- Build Plot ---
  g <- ggplot(data = poly_data) +
    geom_polygon(aes(x = x, y = y, group = group, fill = fill), 
                 color = border_color,
                 linewidth = 0.2) +
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right")
  
  # --- Adaptive Color Scale ---
  if (is.numeric(fill_values)) {
    g <- g + scale_fill_viridis_c(option = "turbo")
  } else {
    g <- g + scale_fill_discrete()
  }
  
  return(g)
}

# -------------------------------------------------------------------------
# 2. S3 Generic CartFeaturePlot
# -------------------------------------------------------------------------
#' CartFeaturePlot S3 generic
#' @param object Object to plot (SpatialMesh or Seurat)
#' @param features Feature name(s) to display
#' @param layout.name Name of the layout to use
#' @param title Plot title
#' @param ... Additional arguments
CartFeaturePlot <- function(object, 
                            features = NULL, 
                            layout.name = NULL, 
                            title = NULL,
                            cell.limit = NULL,
                            ...) {
  UseMethod("CartFeaturePlot")
}

# -------------------------------------------------------------------------
# 3. Default Method for SpatialMesh
# -------------------------------------------------------------------------
#' Default CartFeaturePlot method for SpatialMesh objects
#' Features are primarily sought in vertex_metadata slot
CartFeaturePlot.SpatialMesh <- function(object, 
                                         features = NULL, 
                                         layout.name = NULL, 
                                         title = NULL,
                                         cell.limit = NULL,
                                         ...) {
  
  mesh <- object
  
  # --- B. Auto-detect layout.name if NULL ---
  if (is.null(layout.name)) {
    if (is.null(mesh@cartogram) || length(mesh@cartogram) == 0) {
      stop("No power diagrams found in mesh@cartogram. Cannot auto-detect layout.name.")
    }
    layout.name <- names(mesh@cartogram)[1]
    message(paste0("Auto-detected layout.name: '", layout.name, "'"))
  }
  
  # --- C. Extract shape from layout.name ---
  # Extract shape as the string after the last period, or last two parts if last is "boundary"
  layout_parts <- strsplit(layout.name, ".", fixed = TRUE)[[1]]
  if (length(layout_parts) > 0) {
    if (length(layout_parts) > 1 && layout_parts[length(layout_parts)] == "boundary") {
      # If last part is "boundary", take last two parts (e.g., "free_boundary")
      shape <- paste(layout_parts[(length(layout_parts)-1):length(layout_parts)], collapse = ".")
    } else {
      # Otherwise, take the last part (e.g., "disk")
      shape <- layout_parts[length(layout_parts)]
    }
  } else {
    shape <- NULL
  }
  
  # --- D. Parameter Validation & Defaults ---
  if (is.null(features)) {
    features <- paste0("area.", layout.name)
  }
  
  # --- E. Validate power diagram exists ---
  if (is.null(mesh@cartogram) || !layout.name %in% names(mesh@cartogram)) {
    stop(paste0("Power diagram for layout '", layout.name, "' not found."))
  }
  
  cartogram <- mesh@cartogram[[layout.name]]
  
  # --- F. Validate power diagram and extract coordinates ---
  if (is.null(cartogram$uv)) {
    stop(paste0("Power diagram for layout '", layout.name, "' does not contain 'uv' coordinates."))
  }
  
  # --- G. Compute Boundary (Single Computation) ---
  # Only compute boundary if shape is not "free_boundary"
  skip_boundary_clipping <- !is.null(shape) && shape == "free_boundary"
  
  if (!skip_boundary_clipping) {
    boundary_indices <- mesh@boundary
    boundary_coords <- cartogram$uv[boundary_indices, , drop = FALSE]
  }
  
  # --- H. Process Power Diagram Cells (Memory-Efficient Loop) ---
  n_cells <- length(cartogram$cell)
  if (n_cells == 0) {
    stop("Power diagram contains no cells.")
  }
  
  # Apply cell.limit filtering if specified
  if (!is.null(cell.limit)) {
    n_cells <- min(n_cells, cell.limit)
  }
  
  # Pre-extract dpe for efficiency (avoid repeated list access)
  dpe <- cartogram$dpe
  
  # Pre-allocate list for better memory management
  layout_list <- vector("list", length = n_cells)
  
  # Process cells with optional boundary clipping
  for (i in seq_len(n_cells)) {
    cell_indices <- cartogram$cell[[i]]
    
    if (length(cell_indices) == 0) {
      layout_list[[i]] <- matrix(numeric(0), ncol = 2)
      next
    }
    
    # Extract cell coordinates
    cell_coords <- dpe[cell_indices, , drop = FALSE]
    
    # Clip cell with boundary (skip if shape == "free_boundary")
    if (skip_boundary_clipping) {
      # Use cell coordinates directly without clipping
      layout_list[[i]] <- cell_coords
    } else {
      # Clip cell with boundary
      pb_result <- polybool(cell_coords, boundary_coords, "and")
      
      # Extract clipped polygon (handle nested list structure)
      if (length(pb_result) > 0 && length(pb_result[[1]]) > 0) {
        layout_list[[i]] <- pb_result[[1]][[1]]
      } else {
        layout_list[[i]] <- matrix(numeric(0), ncol = 2)
      }
    }
  }
  
  # --- I. Fetch Feature Data from vertex_metadata ---
  if (!features %in% colnames(mesh@vertex_metadata)) {
    stop(paste0("Feature '", features, "' not found in mesh@vertex_metadata."))
  }
  fill_values <- mesh@vertex_metadata[[features]]
  fill_values <- as.vector(fill_values)
  
  # Apply cell.limit to fill_values if specified
  if (!is.null(cell.limit)) {
    fill_values <- fill_values[seq_len(min(length(fill_values), cell.limit))]
  }
  
  # --- Truncate cell list to match feature length if needed ---
  n_fill <- length(fill_values)
  if (n_fill < n_cells) {
    message(paste0("Feature length (", n_fill, 
                   ") is smaller than number of cells (", n_cells, 
                   "). Truncating to first ", n_fill, " cells."))
    layout_list <- layout_list[seq_len(n_fill)]
    n_cells <- n_fill
  } else if (n_fill > n_cells) {
    warning(paste0("Feature length (", n_fill, 
                   ") is larger than number of cells (", n_cells, 
                   "). Using first ", n_cells, " feature values."))
    fill_values <- fill_values[seq_len(n_cells)]
  }
  
  # --- J. Build Plot Using Core Function ---
  plot <- .CartFeaturePlot_Core(layout_list = layout_list, 
                                 fill_values = fill_values)
  
  # --- K. Add Labels ---
  plot <- plot +
    labs(title = title,
         x = paste0(layout.name, "_2"),
         y = paste0(layout.name, "_1"),
         fill = features)
  
  return(plot)
}

# -------------------------------------------------------------------------
# 4. Method for Seurat Objects
# -------------------------------------------------------------------------
#' CartFeaturePlot method for Seurat objects
#' First seeks feature in spatial mesh vertex_metadata, then tries FetchData()
CartFeaturePlot.Seurat <- function(object, 
                                    features = NULL, 
                                    layout.name = NULL, 
                                    title = NULL,
                                    cell.limit = NULL,
                                    ...) {
  
  # --- A. Get SpatialMesh ---
  mesh <- SpatialMesh(object)
  
  # Apply default cell.limit from Seurat object if not specified
  if (is.null(cell.limit)) {
    cell.limit <- ncol(object)
  }
  
  # --- B. Auto-detect layout.name if NULL ---
  if (is.null(layout.name)) {
    if (is.null(mesh@cartogram) || length(mesh@cartogram) == 0) {
      stop("No power diagrams found in mesh@cartogram. Cannot auto-detect layout.name.")
    }
    layout.name <- names(mesh@cartogram)[1]
    message(paste0("Auto-detected layout.name: '", layout.name, "'"))
  }
  
  # --- C. Extract shape from layout.name ---
  layout_parts <- strsplit(layout.name, ".", fixed = TRUE)[[1]]
  if (length(layout_parts) > 0) {
    if (length(layout_parts) > 1 && layout_parts[length(layout_parts)] == "boundary") {
      shape <- paste(layout_parts[(length(layout_parts)-1):length(layout_parts)], collapse = ".")
    } else {
      shape <- layout_parts[length(layout_parts)]
    }
  } else {
    shape <- NULL
  }
  
  # --- D. Parameter Validation & Defaults ---
  if (is.null(features)) {
    features <- paste0("area.", layout.name)
  }
  
  # --- E. Validate power diagram exists ---
  if (is.null(mesh@cartogram) || !layout.name %in% names(mesh@cartogram)) {
    stop(paste0("Power diagram for layout '", layout.name, "' not found."))
  }
  
  cartogram <- mesh@cartogram[[layout.name]]
  
  # --- F. Validate power diagram and extract coordinates ---
  if (is.null(cartogram$uv)) {
    stop(paste0("Power diagram for layout '", layout.name, "' does not contain 'uv' coordinates."))
  }
  
  # --- G. Compute Boundary (Single Computation) ---
  skip_boundary_clipping <- !is.null(shape) && shape == "free_boundary"
  
  if (!skip_boundary_clipping) {
    boundary_indices <- mesh@boundary
    boundary_coords <- cartogram$uv[boundary_indices, , drop = FALSE]
  }
  
  # --- H. Process Power Diagram Cells (Memory-Efficient Loop) ---
  n_cells <- length(cartogram$cell)
  if (n_cells == 0) {
    stop("Power diagram contains no cells.")
  }
  
  # Apply cell.limit filtering if specified
  if (!is.null(cell.limit)) {
    n_cells <- min(n_cells, cell.limit)
  }
  
  dpe <- cartogram$dpe
  layout_list <- vector("list", length = n_cells)
  
  for (i in seq_len(n_cells)) {
    cell_indices <- cartogram$cell[[i]]
    
    if (length(cell_indices) == 0) {
      layout_list[[i]] <- matrix(numeric(0), ncol = 2)
      next
    }
    
    cell_coords <- dpe[cell_indices, , drop = FALSE]
    
    if (skip_boundary_clipping) {
      layout_list[[i]] <- cell_coords
    } else {
      pb_result <- polybool(cell_coords, boundary_coords, "and")
      if (length(pb_result) > 0 && length(pb_result[[1]]) > 0) {
        layout_list[[i]] <- pb_result[[1]][[1]]
      } else {
        layout_list[[i]] <- matrix(numeric(0), ncol = 2)
      }
    }
  }
  
  # --- I. Fetch Feature Data: Try vertex_metadata first, then FetchData ---
  fill_values <- NULL
  if (features %in% colnames(mesh@vertex_metadata)) {
    fill_values <- mesh@vertex_metadata[[features]]
    fill_values <- as.vector(fill_values)
  } else {
    # Try FetchData from Seurat object
    fill_data <- FetchData(object, vars = features)
    if (ncol(fill_data) == 0) {
      stop(paste0("Feature '", features, "' not found in mesh@vertex_metadata or in Seurat object."))
    }
    fill_values <- fill_data[[1]]  # Extract as vector directly
  }
  
  # Apply cell.limit to fill_values if specified
  if (!is.null(cell.limit)) {
    fill_values <- fill_values[seq_len(min(length(fill_values), cell.limit))]
  }
  
  # --- J. Truncate cell list to match feature length if needed ---
  n_fill <- length(fill_values)
  if (n_fill < n_cells) {
    message(paste0("Feature length (", n_fill, 
                   ") is smaller than number of cells (", n_cells, 
                   "). Truncating to first ", n_fill, " cells."))
    layout_list <- layout_list[seq_len(n_fill)]
    n_cells <- n_fill
  } else if (n_fill > n_cells) {
    warning(paste0("Feature length (", n_fill, 
                   ") is larger than number of cells (", n_cells, 
                   "). Using first ", n_cells, " feature values."))
    fill_values <- fill_values[seq_len(n_cells)]
  }
  
  # --- K. Build Plot Using Core Function ---
  plot <- .CartFeaturePlot_Core(layout_list = layout_list, 
                                 fill_values = fill_values)
  
  # --- L. Add Labels ---
  plot <- plot +
    labs(title = title,
         x = paste0(layout.name, "_2"),
         y = paste0(layout.name, "_1"),
         fill = features)
  
  return(plot)
}
