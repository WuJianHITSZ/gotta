#' Run Trajectory
#' 
#' Wrapper function around GradientLayout that computes trajectory vectors
#' for visualizing feature gradients in layout space.
#' 
#' @param object A SpatialMesh or Seurat object
#' @param alignment.job GOTTJob or AlignmentJob object. If provided, parses source.assay, target.assay, and shape from it (default: NULL)
#' @param target.assay Character. Target assay name for the feature (default: NULL)
#' @param source.assay Character. Source assay name for the layout (default: "spatial")
#' @param shape Character. Shape type (default: "free_boundary")
#' @param sigma Numeric. Optional smoothing parameter for gradient computation
#' @param ... Unused
#' @return SpatialMesh or Seurat object with trajectory stored
#' @export
RunTrajectory <- function(object, alignment.job = NULL, target.assay = NULL,
                          source.assay = "spatial", shape = "free_boundary",
                          sigma = NULL, ...) {
  is_seurat <- inherits(object, "Seurat")
  mesh <- if (is_seurat) SpatialMesh(object) else object
  
  # Normalize inputs: prefer alignment.job, otherwise build it from target/source/shape.
  if (is.null(alignment.job)) {
    if (!is.null(target.assay) && !is.null(source.assay) && !is.null(shape)) {
      alignment.job <- AlignmentJob(shape = shape,
                                    source.assay = source.assay,
                                    target.assay = target.assay)
    } else {
      stop("alignment.job must be provided, or supply target.assay, source.assay, and shape")
    }
  }
  
  # Parse alignment.job
  if (inherits(alignment.job, "GOTTJob")) {
    source.assay <- "spatial"
    target.assay <- alignment.job$vertex.job$assay
    shape <- alignment.job$shape
  } else if (inherits(alignment.job, "AlignmentJob")) {
    source.assay <- alignment.job$source.assay
    target.assay <- alignment.job$target.assay
    shape <- alignment.job$shape
  } else {
    stop("alignment.job must be a GOTTJob or AlignmentJob object")
  }
  
  # Validate inputs
  if (!methods::is(mesh, "SpatialMesh")) {
    stop("mesh must be a SpatialMesh object")
  }
  if (!is.character(target.assay) || length(target.assay) != 1) {
    stop("target.assay must be a single character string")
  }
  if (!is.character(source.assay) || length(source.assay) != 1) {
    stop("source.assay must be a single character string")
  }
  if (!is.character(shape) || length(shape) != 1) {
    stop("shape must be a single character string")
  }
  
  # Determine layout.name based on source.assay and shape
  if (source.assay == "spatial" && shape == "free_boundary") {
    layout.name <- "spatial_coords"
  } else if (source.assay == "spatial") {
    layout.name <- paste0("spatial.", shape)
  } else {
    layout.name <- paste0(source.assay, ".", shape)
  }
  
  # Set feature.name
  if (source.assay == "spatial") {
    feature.name <- paste0("h.", target.assay, ".", shape)
  } else {
    feature.name <- paste0("h.", source.assay, "2", target.assay, ".", shape)
  }
  
  
  # Log message indicating feature.name and layout.name
  message(sprintf("[RunTrajectory] Computing trajectory with feature.name='%s' and layout.name='%s'", 
                  feature.name, layout.name))
  
  # Call gradient_layout
  result <- gradient_layout(mesh = mesh, 
                            feature_name = feature.name, 
                            layout_name = layout.name, 
                            sigma = sigma)
  
  # Store trajectory in mesh@layout with naming convention: arrow.<source.assay>2<target.assay>.<shape>
  trajectory_name <- paste0("arrow.", source.assay, "2", target.assay, ".", shape)
  mesh@layout[[trajectory_name]] <- result$trajectory


  # Choose gradient/Jacobian strategy for divergence/curl/distortion.
  if (source.assay == "spatial" && shape == "free_boundary") {
    message("[RunTrajectory] Using GradientLayout for spatial/free_boundary Jacobian metrics")
    layout.name <- paste0(target.assay, ".", shape)
    output <- gradient_layout(mesh, feature_name = feature.name, layout_name = layout.name)
  } else {
    xy_name <- paste0(source.assay, ".", shape)
    uv_name <- paste0(target.assay, ".", shape)
    message(sprintf("[RunTrajectory] Using JacobianLayout with xy_name='%s', uv_name='%s'", xy_name, uv_name))
    output <- jacobian_layout(mesh, xy_name = xy_name, uv_name = uv_name, sigma = sigma)
  }

  # Store Jacobian-derived metrics using standard naming.
  divergence_name <- paste0("divergence.", source.assay, "2", target.assay, ".", shape)
  mesh@vertex_metadata[[divergence_name]] <- output$divergence
  curl_name <- paste0("curl.", source.assay, "2", target.assay, ".", shape)
  mesh@vertex_metadata[[curl_name]] <- output$curl
  distortion_name <- paste0("distortion.", source.assay, "2", target.assay, ".", shape)
  mesh@vertex_metadata[[distortion_name]] <- output$distortion
  
  if (is_seurat) {
    SpatialMesh(object) <- mesh
    return(object)
  }
  return(mesh)
}
