#' Compute Geometric Optimal Transport Tableau & Alignment
#' 
#' Performs discrete optimal transport and boundary alignment.
#' @param face Face matrix (Mx3 or Mx4)
#' @param uv_initial Initial UV coordinates (Nx2)
#' @param sigma Function for computing target ratios
#' @param corner Corner indices (NULL for disk, 1+ for rect)
#' @param shape Shape type: "disk" or "rect"
#' @param aspect.ratio Aspect ratio for shape (default: 1)
#' @param verbose Logical. Whether to print progress messages (default: FALSE)
#' @return Power diagram list with aligned UV coordinates
#' @noRd
GOTTA <- function(face, uv_initial, sigma, corner = NULL, shape = "disk", aspect.ratio = 1, verbose = FALSE){
  # Compute boundary indices once
  if (verbose) message("[GOTTA] Computing boundary indices...")
  boundary_index <- compute_bd(face)
  
  # Compute initial area distribution
  if (verbose) message("[GOTTA] Computing initial area distribution...")
  area_initial <- vertex_area(face, uv_initial) / 3
  
  # Extract boundary coordinates
  boundary <- uv_initial[boundary_index, , drop = FALSE]
  
  # Compute discrete optimal transport
  if (verbose) message("[GOTTA] Running discrete optimal transport...")
  cartogram <- discrete_optimal_transport(boundary, face, uv_initial, sigma, area_initial)
  
  # Update centroids
  if (verbose) message("[GOTTA] Computing centroids...")
  uv_align <- compute_centroid(boundary, cartogram)
  
  # Clean up boundary based on shape
  if (verbose) message("[GOTTA] Cleaning up boundary...")
  uv_align <- cleanup_boundary(uv_align, boundary_index, corner, shape, aspect.ratio, uv_initial)
  
  # Update power diagram and clean up
  if (verbose) message("[GOTTA] Finalizing cartogram...")
  cartogram$uv <- uv_align
  cartogram$dp <- NULL
  cartogram$face <- NULL
  
  return(cartogram)
}

#' Run Geometric Optimal Transport Tableau & Alignment
#' 
#' Aligns source and target assays using optimal transport on topological structures.
#' @param object Seurat object with TopoStruct or SpatialMesh object
#' @param shape Shape type: "disk" or "rect" (default: "disk")
#' @param source.assay Source assay/reduction name
#' @param target.assay Target assay/reduction name
#' @param target.object Target Seurat object (NULL = use source object)
#' @param target.name Target name for cross-object alignment
#' @param is.spatial.initial Use spatial coordinates as initial (default: FALSE = use vertex coordinates)
#' @param is.hyper.target Use hyper-parameter target (default: FALSE)
#' @return Modified Seurat object with alignment results
#' @export
RunGOTTA <- function(object, ...) {
  UseMethod("RunGOTTA")
}

#' Default RunGOTTA method for SpatialMesh objects
#' @noRd
RunGOTTA.SpatialMesh <- function(object, shape = "disk", 
                                  source.assay = NULL, 
                                  target.assay = NULL, 
                                  target.object = NULL, 
                                  target.name = NULL,
                                  is.spatial.initial = FALSE,
                                  is.hyper.target = FALSE,
                                  alignment.job = NULL,
                                  verbose = FALSE,
                                  ...) {
  if (!is.null(alignment.job)) {
    shape <- alignment.job$shape
    source.assay <- alignment.job$source.assay
    target.assay <- alignment.job$target.assay
    target.object <- alignment.job$target.object
    target.name <- alignment.job$target.name
    is.spatial.initial <- alignment.job$is.spatial.initial
    is.hyper.target <- alignment.job$is.hyper.target
  }
  if (verbose) message("[RunGOTTA] Starting RunGOTTA...")
  if (verbose) message(sprintf("[RunGOTTA] Parameters: shape=%s, source.assay=%s, target.assay=%s, target.object=%s, target.name=%s, is.spatial.initial=%s, is.hyper.target=%s", shape, source.assay, target.assay, target.object, target.name, is.spatial.initial, is.hyper.target))
  
  # ============================================================================
  # Step 1: Parameter Setup and Validation
  # ============================================================================
  if (verbose) message("[RunGOTTA] Setting up parameters and validating inputs...")
  mesh <- object
  
  # For SpatialMesh, target.object should also be a SpatialMesh or NULL
  if (is.null(target.object)) {
    mesh.target <- mesh
    align.assay <- paste0(source.assay, "2", target.assay)
    if (verbose) message("[RunGOTTA] Using same mesh as source and target")
  } else if (methods::is(target.object, "SpatialMesh")) {
    mesh.target <- target.object
    if(is.null(target.assay)){
      target.assay <- source.assay
    }
    align.assay <- paste0(source.assay, "2", target.name)
    if (verbose) message("[RunGOTTA] Using separate target mesh")
  } else {
    stop("When object is SpatialMesh, target.object must be NULL or SpatialMesh")
  }
  
  # Build layout names
  layout_name <- paste0(align.assay, ".", shape)
  source.assay <- paste0(source.assay, ".", shape)
  spatial_shape <- paste0("spatial.", shape)
  layout_name.target <- paste0(target.assay, ".", shape)
  
  # Validate inputs
  if(!source.assay %in% names(mesh@cartogram)){
    stop("source.assay '", source.assay, "' not found in mesh@cartogram")
  }
  
  if(!layout_name.target %in% names(mesh.target@cartogram) && 
     !layout_name.target %in% names(mesh.target@layout)){
    stop("layout_name.target '", layout_name.target, "' not found in target mesh@cartogram or mesh@layout")
  }
  
  # ============================================================================
  # Step 2: Extract Target Data
  # ============================================================================
  if (verbose) message("[RunGOTTA] Extracting target data...")
  face.target <- mesh.target@face
  uv_spatial.target <- GetLayout(mesh.target, spatial_shape)
  uv_vertex.target <- GetLayout(mesh.target, layout_name.target)

  # ============================================================================
  # Step 3: Compute Area Ratios for Interpolation
  # ============================================================================
  if (verbose) message("[RunGOTTA] Computing area ratios for interpolation...")
  # Compute vertex area (use vertex area if is.hyper.target, else use UV vertex data)
  area_key <- if (is.hyper.target) paste0("area.", target.assay) else paste0("area.", target.assay, ".", shape)
  area_vertex <- mesh.target@vertex_metadata[[area_key]]

  
  # Compute spatial coordinate area
  area_spatial <- vertex_area(face.target, uv_spatial.target) / 3
  # Normalize areas (preserve total area)
  sum_area_vertex <- sum(area_vertex)
  sum_area_spatial <- sum(area_spatial)
  area_vertex <- area_vertex / sum_area_vertex * sum_area_spatial
  
  # ============================================================================
  # Step 4: Build Scattered Interpolant Model
  # ============================================================================
  if (verbose) message("[RunGOTTA] Building scattered interpolant model...")
  area_ratio <- area_spatial / area_vertex
  mdl_lin <- si_create_linear(uv_vertex.target, area_ratio, face.target)
  
  # Create sigma function for GOTTA
  sigma <- function(coords) {
    si_predict_linear(mdl_lin, coords)
  }
  
  # ============================================================================
  # Step 5: Extract Source Data and Run Alignment
  # ============================================================================
  if (verbose) message("[RunGOTTA] Extracting source data and running alignment...")
  face <- mesh@face
  uv_spatial <- GetLayout(mesh, spatial_shape)
  uv_vertex <- GetLayout(mesh, source.assay)
  corner <- mesh@corner
  
  # Choose initial coordinates
  uv_initial <- if(is.spatial.initial) uv_spatial else uv_vertex
  if (verbose) message(sprintf("[RunGOTTA] Using %s as initial coordinates", if(is.spatial.initial) "spatial coordinates" else "vertex coordinates"))
  
  # Run GOTTA alignment
  if (verbose) message("[RunGOTTA] Running GOTTA alignment...")
  cartogram <- GOTTA(face, uv_initial, sigma, corner, shape, aspect.ratio = 1, verbose = verbose)
  
  # ============================================================================
  # Step 6: Store Results
  # ============================================================================
  if (verbose) message("[RunGOTTA] Storing results...")
  mesh@layout[[layout_name]] <- cartogram$uv
  mesh@cartogram[[layout_name]] <- cartogram
  if (verbose) message("[RunGOTTA] Updating mesh metadata...")
  mesh <- update_mesh_metadata(mesh, face, cartogram$uv, cartogram, layout_name)

  if (verbose) message("[RunGOTTA] RunGOTTA finished.")
  mesh
}

#' RunGOTTA method for Seurat objects
#' @noRd
RunGOTTA.Seurat <- function(object, 
                            shape = "disk", 
                            source.assay = NULL, 
                            target.assay = NULL, 
                            target.object = NULL, 
                            target.name = NULL,
                            is.spatial.initial = FALSE,
                            is.hyper.target = FALSE,
                            alignment.job = NULL,
                            verbose = FALSE,
                            ...) {
  
  if (is.null(alignment.job)) {
    if (verbose) message("[RunGOTTA.Seurat] Creating AlignmentJob from parameters...")
    alignment.job <- AlignmentJob(shape = shape,
                                  source.assay = source.assay,
                                  target.assay = target.assay,
                                  target.object = target.object,
                                  target.name = target.name,
                                  is.spatial.initial = is.spatial.initial,
                                  is.hyper.target = is.hyper.target)
  }
  
  if (verbose) message("[RunGOTTA.Seurat] Extracting SpatialMesh from Seurat object...")
  mesh <- SpatialMesh(object)
  
  # Handle target object - if it's a Seurat object, extract its mesh
  if (!is.null(target.object)) {
    if (methods::is(target.object, "Seurat")) {
      if (verbose) message("[RunGOTTA.Seurat] Extracting SpatialMesh from target Seurat object...")
      target.mesh <- SpatialMesh(target.object)
    } else if (methods::is(target.object, "SpatialMesh")) {
      if (verbose) message("[RunGOTTA.Seurat] Using provided SpatialMesh as target...")
      target.mesh <- target.object
    } else {
      stop("target.object must be NULL, Seurat, or SpatialMesh")
    }
  } else {
    # When target.object is NULL, use the same object as target
    if (verbose) message("[RunGOTTA.Seurat] Using same object as target...")
    target.mesh <- NULL
  }
  
  # Validate that source.assay exists in mesh@cartogram or mesh@layout
  # and target.assay exists in target.mesh@cartogram or target.mesh@layout
  if (!is.null(source.assay)) {
    if (verbose) message("[RunGOTTA.Seurat] Validating source assay...")
    source_assay_found <- any(c(
      paste0(source.assay, ".", shape) %in% names(mesh@cartogram),
      paste0(source.assay, ".", shape) %in% names(mesh@layout)
    ))
    if (!source_assay_found) {
      stop("source.assay '", source.assay, "' not found in source mesh cartogram or layout")
    }
  }
  if (!is.null(target.mesh) && !is.null(target.assay)) {
    if (verbose) message("[RunGOTTA.Seurat] Validating target assay...")
    target_assay_found <- any(c(
      paste0(target.assay, ".", shape) %in% names(target.mesh@cartogram),
      paste0(target.assay, ".", shape) %in% names(target.mesh@layout)
    ))
    if (!target_assay_found) {
      stop("target.assay '", target.assay, "' not found in target mesh cartogram or layout")
    }
  }
  
  # Update alignment.job with target.mesh for SpatialMesh method
  alignment.job$target.object <- target.mesh
  
  # Call the SpatialMesh method
  if (verbose) message("[RunGOTTA.Seurat] Calling RunGOTTA.SpatialMesh...")
  mesh_updated <- RunGOTTA(mesh, 
                           alignment.job = alignment.job,
                           verbose = verbose)
  
  # Assign updated mesh back to Seurat object
  if (verbose) message("[RunGOTTA.Seurat] Assigning updated mesh back to Seurat object...")
  SpatialMesh(object) <- mesh_updated
  if (verbose) message("[RunGOTTA.Seurat] RunGOTTA finished.")
  object
}

