GradientSpatial <- function(mesh, feature.name, sigma = NULL) {
  # Validate inputs
  if (!methods::is(mesh, "SpatialMesh")) {
    stop("mesh must be a SpatialMesh object")
  }
  if (!is.character(feature.name) || length(feature.name) != 1) {
    stop("feature.name must be a single character string")
  }
  
  # Check dependencies
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The 'imager' package is required. Please install it with install.packages('imager').")
  }
  
  # Get feature from vertex_metadata or layout
  if (feature.name %in% colnames(mesh@vertex_metadata)) {
    feature <- mesh@vertex_metadata[[feature.name]]
  } else {
    feature <- GetLayout(mesh, feature.name)
  }
  
  # Validate feature
  if (is.null(feature) || length(feature) == 0) {
    stop(sprintf("Feature or layout '%s' not found in mesh@vertex_metadata or mesh@layout", feature.name))
  }
  
  # Convert to matrix if vector
  if (is.vector(feature)) {
    feature <- as.matrix(feature)
  }
  
  # Get coordinates and boundary
  array_coords <- GetLayout(mesh, "array_coords")
  layout <- GetLayout(mesh, "spatial_coords")
  boundary <- mesh@boundary
  
  # Validate dimensions
  if (nrow(feature) != nrow(array_coords)) {
    stop(sprintf("Feature has %d rows but array_coords has %d rows. They must match.", 
                 nrow(feature), nrow(array_coords)))
  }
  
  # very important! no changing!
  Dy <- matrix(c(0, 0, 0, -1, 0, 1, 0, 0, 0), nrow = 3, byrow = TRUE) / 2
  Dx <- matrix(c(-1, -1, 0, 0, 0, 0, 0, 1, 1), nrow = 3, byrow = TRUE) / 2 / sqrt(3)
  
  # Pre-compute kernel cimg objects (reused for all columns)
  Dx_cimg <- imager::as.cimg(Dx)
  Dy_cimg <- imager::as.cimg(Dy)
  
  # Compute image dimensions and indices
  imsize <- apply(array_coords, 2, max)
  n_features <- ncol(feature)
  if (is.null(n_features)) n_features <- 1L
  
  indices <- array_coords[, 1L] + (array_coords[, 2L] - 1L) * imsize[1L]
  n_points <- length(indices)
  
  # Pre-allocate output matrices (more memory efficient than lists)
  image_x_final <- matrix(0, nrow = n_points, ncol = n_features)
  image_y_final <- matrix(0, nrow = n_points, ncol = n_features)
  
  # Process each feature column
  for (k in seq_len(n_features)) {
    feat_col <- if (n_features > 1L) feature[, k] else feature
    
    # Create image matrix and convert to cimg
    img_mat <- matrix(0, nrow = imsize[1L], ncol = imsize[2L])
    img_mat[indices] <- feat_col
    img_cimg <- imager::as.cimg(img_mat)
    
    # Compute gradients with correlation
    img_x_cimg <- imager::correlate(img_cimg, Dx_cimg, dirichlet = FALSE)
    img_y_cimg <- imager::correlate(img_cimg, Dy_cimg, dirichlet = FALSE)
    
    # Apply smoothing if requested
    if (!is.null(sigma) && length(sigma) > 0 && sigma > 0) {
      img_x_cimg <- imager::isoblur(img_x_cimg, sigma)
      img_y_cimg <- imager::isoblur(img_y_cimg, sigma)
    }
    
    # Extract values at indices
    img_x_mat <- as.matrix(img_x_cimg)
    img_y_mat <- as.matrix(img_y_cimg)
    image_x_final[, k] <- img_x_mat[indices]
    image_y_final[, k] <- img_y_mat[indices]
  }
  
  # Apply boundary conditions
  if (length(boundary) > 0L) {
    image_x_final[boundary, ] <- 0
    image_y_final[boundary, ] <- 0
  }
  
  # Combine x and y gradients
  trajectory <- cbind(image_x_final, image_y_final)
  
  return(list(trajectory = trajectory, layout = layout))
}


GradientLayout <- function(mesh, feature.name, layout.name, sigma = NULL) {
  # Validate inputs
  if (!methods::is(mesh, "SpatialMesh")) {
    stop("mesh must be a SpatialMesh object")
  }
  if (!is.character(feature.name) || length(feature.name) != 1) {
    stop("feature.name must be a single character string")
  }
  if (!is.character(layout.name) || length(layout.name) != 1) {
    stop("layout.name must be a single character string")
  }
  
  # Check dependencies
  if (!requireNamespace("Rlinsolve", quietly = TRUE)) {
    stop("The 'Rlinsolve' package is required. Please install it with install.packages('Rlinsolve').")
  }
  
  # Get feature and layout
  feature <- mesh@vertex_metadata[[feature.name]]
  if (is.null(feature)) {
    stop(sprintf("Feature '%s' not found in mesh@vertex_metadata", feature.name))
  }
  
  layout <- GetLayout(mesh, layout.name)
  boundary <- mesh@boundary
  
  # Calculate gradients using GradientSpatial
  # trajectory_feature: [df/dx, df/dy] (assuming feature is scalar)
  result_feature <- GradientSpatial(mesh, feature.name, sigma)
  trajectory_feature <- result_feature$trajectory
  
  # trajectory_layout: [du/dx, dv/dx, du/dy, dv/dy] (assuming layout is 2D)
  result_layout <- GradientSpatial(mesh, layout.name, sigma)
  trajectory_layout <- result_layout$trajectory
  
  # Validate dimensions match
  n_points <- nrow(trajectory_layout)
  if (nrow(trajectory_feature) != n_points) {
    stop(sprintf("Dimension mismatch: trajectory_feature has %d rows but trajectory_layout has %d rows",
                 nrow(trajectory_feature), n_points))
  }
  
  # Validate trajectory_layout has 4 columns (for 2D layout: du/dx, dv/dx, du/dy, dv/dy)
  if (ncol(trajectory_layout) != 4L) {
    stop(sprintf("trajectory_layout must have 4 columns (got %d). Layout must be 2D.", ncol(trajectory_layout)))
  }
  
  # Validate trajectory_feature has 2 columns (for scalar feature: df/dx, df/dy)
  if (ncol(trajectory_feature) != 2L) {
    stop(sprintf("trajectory_feature must have 2 columns (got %d). Feature must be scalar.", ncol(trajectory_feature)))
  }
  
  # Pre-allocate output matrices
  trajectory <- matrix(0, nrow = n_points, ncol = 2L)
  distortion <- numeric(n_points)
  divergence <- numeric(n_points)
  curl <- numeric(n_points)
  # Solve linear systems for each point
  # The chain rule: [df/dx]   [du/dx  dv/dx] [df/du]
  #                 [df/dy] = [du/dy  dv/dy] [df/dv]
  # We solve Jacobi * [df/du, df/dv]^T = [df/dx, df/dy]^T
  # where Jacobi = [du/dx  dv/dx; du/dy  dv/dy]
  for (i in seq_len(n_points)) {
    # Construct Jacobian matrix Jacobi (2x2) from trajectory_layout
    # trajectory_layout[i, ] = [du/dx, dv/dx, du/dy, dv/dy]
    Jacobi <- matrix(trajectory_layout[i, ], nrow = 2L, ncol = 2L, byrow = TRUE)
    
    # Get right-hand side vector b = [df/dx, df/dy]
    b <- matrix(trajectory_feature[i, ], ncol = 1L)
    
    # Solve linear system Jacobi * x = b using Conjugate Gradient Squared method
    res <- Rlinsolve::lsolve.cgs(Jacobi, b, verbose = FALSE)
    
    # Store solution [df/du, df/dv] and Jacobian determinant
    trajectory[i, ] <- res$x
    distortion[i] <- det(Jacobi)
    divergence[i] <- sum(diag(Jacobi))
    curl[i] <- Jacobi[1,2] - Jacobi[2,1]
  }
  
  # Apply boundary conditions
  if (length(boundary) > 0L) {
    trajectory[boundary, ] <- 0
    distortion[boundary] <- 0
    divergence[boundary] <- 0
    curl[boundary] <- 0
  }
  
  return(list(
    trajectory = trajectory,
    distortion = distortion,
    layout = layout,
    feature = feature,
    divergence = divergence,
    curl = curl
  ))
}

#' Compute Jacobian Matrix
#' 
#' Computes the Jacobian determinant, divergence, and curl for each point by solving linear systems
#' using gradient information from spatial coordinates and layout coordinates.
#' 
#' @param mesh A SpatialMesh object
#' @param xy_name Character. Name of the layout for xy coordinates
#' @param uv_name Character. Name of the layout for uv coordinates
#' @param sigma Numeric. Optional smoothing parameter for gradient computation
#' @return List containing:
#'   \item{distortion}{Numeric vector of Jacobian determinants (one per point)}
#'   \item{divergence}{Numeric vector of divergence values (one per point)}
#'   \item{curl}{Numeric vector of curl values (one per point)}
#' @export
JacobianLayout <- function(mesh, xy_name, uv_name, sigma = NULL) {
  
  # Get layout matrices
  xy_mat <- GetLayout(mesh, xy_name)
  uv_mat <- GetLayout(mesh, uv_name)
  array_coords <- GetLayout(mesh, "array_coords")
  
  # Get boundary
  boundary <- mesh@boundary
  
  # Compute gradients using GradientSpatial
  # GradientSpatial returns trajectory matrix with 4 columns for 2D layout: [dx/dx, dy/dx, dx/dy, dy/dy]
  result_xy <- GradientSpatial(mesh, xy_name, sigma)
  xy <- result_xy$trajectory
  
  result_uv <- GradientSpatial(mesh, uv_name, sigma)
  uv <- result_uv$trajectory
  
  # Check if Rlinsolve is available
  if (!requireNamespace("Rlinsolve", quietly = TRUE)) {
    stop("The 'Rlinsolve' package is required. Please install it with install.packages('Rlinsolve').")
  }
  
  # Initialize output vectors
  n_points <- nrow(xy)
  distortion <- numeric(n_points)
  divergence <- numeric(n_points)
  curl <- numeric(n_points)
  
  # Process each point
  for (i in seq_len(n_points)) {
    # Reshape xy[i, ] from 1x4 to 2x2 matrix A
    # xy[i, ] = [dx/dx, dy/dx, dx/dy, dy/dy]
    # Reshape to: [dx/dx  dy/dx]
    #             [dx/dy  dy/dy]
    A_vals <- xy[i, ]
    A <- matrix(A_vals, nrow = 2, ncol = 2, byrow = TRUE)
    
    # Reshape uv[i, ] from 1x4 to 2x2 matrix B
    # uv[i, ] = [du/dx, dv/dx, du/dy, dv/dy]
    # Reshape to: [du/dx  dv/dx]
    #             [du/dy  dv/dy]
    B_vals <- uv[i, ]
    B <- matrix(B_vals, nrow = 2, ncol = 2, byrow = TRUE)
    
    # Solve A*Jacobian = B for Jacobian
    # We solve for each row: A*Jacobian[1,]' = B[1,]' and A*Jacobian[2,]' = B[2,]'
    Jacobian <- matrix(0, nrow = 2, ncol = 2)
    
    # Solve for first row of Jacobian
    b1 <- matrix(B[1, ], ncol = 1)  # B[1,] as column vector
    res1 <- Rlinsolve::lsolve.cgs(A, b1, verbose = FALSE)
    Jacobian[1, ] <- res1$x
    
    # Solve for second row of Jacobian
    b2 <- matrix(B[2, ], ncol = 1)  # B[2,] as column vector
    res2 <- Rlinsolve::lsolve.cgs(A, b2, verbose = FALSE)
    Jacobian[2, ] <- res2$x
    
    # Compute determinant (distortion), divergence, and curl of Jacobian
    distortion[i] <- det(Jacobian)
    divergence[i] <- Jacobian[1, 1] + Jacobian[2, 2]  # trace of Jacobian
    curl[i] <- Jacobian[1, 2] - Jacobian[2, 1]
  }
  
  # Apply boundary conditions
  if (length(boundary) > 0L) {
    distortion[boundary] <- 0
    divergence[boundary] <- 0
    curl[boundary] <- 0
  }
  
  return(list(
    distortion = distortion,
    divergence = divergence,
    curl = curl
  ))
}

#' Find Trajectory
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
FindTrajectory <- function(object, alignment.job = NULL, target.assay = NULL,
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
  feature.name <- paste0("h.", target.assay, ".", shape)
  
  # Log message indicating feature.name and layout.name
  message(sprintf("[FindTrajectory] Computing trajectory with feature.name='%s' and layout.name='%s'", 
                  feature.name, layout.name))
  
  # Call GradientLayout
  result <- GradientLayout(mesh = mesh, 
                          feature.name = feature.name, 
                          layout.name = layout.name, 
                          sigma = sigma)
  
  # Store trajectory in mesh@layout with naming convention: arrow.<source.assay>2<target.assay>.<shape>
  trajectory_name <- paste0("arrow.", source.assay, "2", target.assay, ".", shape)
  mesh@layout[[trajectory_name]] <- result$trajectory


  # Choose gradient/Jacobian strategy for divergence/curl/distortion.
  if (source.assay == "spatial" && shape == "free_boundary") {
    message("[FindTrajectory] Using GradientLayout for spatial/free_boundary Jacobian metrics")
    layout.name <- paste0(target.assay, ".", shape)
    output <- GradientLayout(mesh, feature.name = feature.name, layout.name = layout.name)
  } else {
    xy_name <- paste0(source.assay, ".", shape)
    uv_name <- paste0(target.assay, ".", shape)
    message(sprintf("[FindTrajectory] Using JacobianLayout with xy_name='%s', uv_name='%s'", xy_name, uv_name))
    output <- JacobianLayout(mesh, xy_name = xy_name, uv_name = uv_name, sigma = sigma)
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
