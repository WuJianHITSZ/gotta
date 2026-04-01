#' Retrieve SpatialMesh from a Seurat object
SpatialMesh <- function(object) {
  if (methods::is(object, "SpatialMesh")) {
    return(object)
  }
  if (!methods::is(object, "Seurat")) {
    stop("SpatialMesh expects a Seurat or SpatialMesh object.")
  }
  mesh <- object@tools[["spatial.mesh"]]
  if (is.null(mesh) || !methods::is(mesh, "SpatialMesh")) {
    stop("SpatialMesh not found in object@tools[['spatial.mesh']].")
  }
  mesh
}

#' Replacement function for SpatialMesh (enables SpatialMesh(object) <- mesh syntax)
`SpatialMesh<-` <- function(object, value) {
  if (!methods::is(value, "SpatialMesh")) {
    stop("`value` must be a SpatialMesh object.")
  }
  if (methods::is(object, "SpatialMesh")) {
    return(value)
  }
  if (!methods::is(object, "Seurat")) {
    stop("SpatialMesh<- expects a Seurat or SpatialMesh object.")
  }
  object@tools[["spatial.mesh"]] <- value
  object
}

#' Fetch Data from Seurat Object to SpatialMesh
#' 
#' Fetches feature data from a Seurat object and adds it to the vertex_metadata
#' of a SpatialMesh object, handling dimension mismatches appropriately.
#' 
#' @param mesh A SpatialMesh object
#' @param object A Seurat object
#' @param feature.name Character. Name of the feature to fetch from the Seurat object
#' @return SpatialMesh object with the feature added to vertex_metadata
#' @export
FetchData2Mesh <- function(mesh, object, feature.name) {
  # Validate inputs
  if (!methods::is(mesh, "SpatialMesh")) {
    stop("mesh must be a SpatialMesh object")
  }
  if (!is.character(feature.name) || length(feature.name) != 1) {
    stop("feature.name must be a single character string")
  }
  
  # Fetch data from Seurat object
  fetched_data <- FetchData(object, vars = feature.name)
  if (ncol(fetched_data) == 0 || nrow(fetched_data) == 0) {
    stop(sprintf("Feature '%s' not found in Seurat object.", feature.name))
  }
  
  # Extract feature values as vector
  feature_values <- fetched_data[, 1, drop = TRUE]
  
  # Create full-length vector matching mesh@vertex_metadata dimensions
  n_cell <- length(feature_values)
  n_metadata <- nrow(mesh@vertex_metadata)
  values <- rep(NA, n_metadata)
  idx <- seq_len(min(n_cell, n_metadata))
  values[idx] <- feature_values[idx]
  
  # Add to mesh@vertex_metadata
  mesh@vertex_metadata[[feature.name]] <- values
  
  return(mesh)
}

#' Fetch Layout from SpatialMesh
#' 
#' Searches for a layout by name in cartogram first, then in layout.
#' @param mesh A SpatialMesh object
#' @param layout_name Character. Name of the layout to fetch
#' @return Matrix of coordinates (uv) if found in cartogram, or layout coordinates if found in layout
#' @noRd
GetLayout <- function(mesh, layout_name) {
  # First search in cartogram
  if (layout_name %in% names(mesh@cartogram)) {
    return(mesh@cartogram[[layout_name]]$uv)
  }
  
  # If not found, search in layout
  if (layout_name %in% names(mesh@layout)) {
    return(mesh@layout[[layout_name]])
  }
  
  # If not found in either, throw error
  stop("Layout '", layout_name, "' not found in mesh@cartogram or mesh@layout")
}
