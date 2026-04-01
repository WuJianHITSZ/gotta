## Find spatial neighbors with mesh construction (v1.6.7)

# Decode barcodes into array coords (Visium or plain "rowxcol")
ParseArrayCoords <- function(barcode, path = NULL, positions_df = NULL) {
  if (nchar(barcode[1]) >= 18) {
    barcode <- substr(barcode, 1, 18)
    if (!is.null(positions_df)) {
      positions_list <- positions_df
    } else {
      if (is.null(path) || !file.exists(path)) {
        stop("tissue_positions_list.csv not found. Provide positions.path or images with array coords.")
      }
      positions_list <- read.csv(path, header = FALSE)
    }
    coords <- positions_list[, c(3, 4)]
    colnames(coords) <- c("row", "col")
    rownames(coords) <- positions_list[, 1]
    flag <- rownames(coords) %in% barcode
    coords <- coords[flag, , drop = FALSE]
    ord <- match(barcode, rownames(coords))
    coords <- coords[ord, , drop = FALSE]
    coords$col <- floor(coords$col / 2) + ceiling(coords$row / 2)
    coords <- as.matrix(coords)
  } else {
    coords <- as.numeric(unlist(strsplit(barcode, "x"))) + 1
    coords <- matrix(coords, nrow = length(barcode), ncol = 2, byrow = TRUE)
  }
  coords
}

# Try to resolve array coords from image coordinates (Visium)
GetArrayCoordsFromImage <- function(object, images, barcode) {
  if (is.null(images) || is.null(object@images[[images]])) {
    return(NULL)
  }
  coords <- object@images[[images]]@coordinates
  if (is.null(coords)) {
    return(NULL)
  }
  col_candidates <- list(c("row", "col"), c("array_row", "array_col"))
  select <- NULL
  for (cand in col_candidates) {
    if (all(cand %in% colnames(coords))) {
      select <- cand
      break
    }
  }
  if (is.null(select)) {
    return(NULL)
  }
  coords <- as.matrix(coords[, select, drop = FALSE])
  rownames(coords) <- rownames(object@images[[images]]@coordinates)
  if (!all(barcode %in% rownames(coords))) {
    warning("Some barcodes missing in image coordinates; dropping missing.")
  }
  coords <- coords[barcode, , drop = FALSE]
  coords
}

# Resolve default path for tissue_positions_list.csv if available
ResolvePositionsPath <- function(path = NULL) {
  if (!is.null(path)) {
    return(path)
  }
  # Expanded candidates list to handle various working directories (Root, Package Root, Subdirs)
  # This makes it robust across RStudio, Positron, VS Code, and script execution contexts
  candidates <- c(
    file.path(getwd(), "inst", "extdata", "tissue_positions_list.csv"),          # WD is package root
    file.path(getwd(), "gotta", "inst", "extdata", "tissue_positions_list.csv"), # WD is project wrapper
    file.path("../inst", "extdata", "tissue_positions_list.csv"),                # WD is one level deep (e.g., functionality/)
    file.path("../../inst", "extdata", "tissue_positions_list.csv")              # WD is two levels deep
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  
  # Check installed package (try both possible case variations)
  for (pkg in c("GOTTA", "gotta")) {
    pkg_path <- system.file("extdata", "tissue_positions_list.csv", package = pkg)
    if (nzchar(pkg_path)) {
      return(pkg_path)
    }
  }
  NULL
}

# Resolve spatial coords from inputs (Staffli, image, or Visium-derived)
ParseSpatialCoords <- function(object, array_coords, col.names = NULL, images = NULL, is.visium = TRUE) {
  if (!is.null(col.names)) {
    coords <- as.matrix(object@tools[["Staffli"]]@meta.data[, col.names])
  } else if (!is.null(images)) {
    coords <- as.matrix(object@images[[images]]@coordss[, c("imagerow", "imagecol")])
  } else if (is.visium) {
    coords <- array_coords
    coords[, 2] <- coords[, 2] - ceiling(coords[, 1] / 2)
  } else {
    stop("Spatial coordss not provided and not Visium.")
  }
  coords
}

# Derive four corners evenly along boundary indices
BoundaryBuilder <- function(boundary) {
  step <- floor(length(boundary) / 4)
  corner <- boundary[seq(1, length(boundary), by = step)]
  corner[1:4]
}

MeshBuilder <- function(array_coords, mesh_type = "quad") {
  
  # 1. Input Handling
  # Ensure input is a matrix for consistent indexing
  if (is.data.frame(array_coords)) array_coords <- as.matrix(array_coords)
  
  # Adjust to 1-based indexing (assuming input is 0-based grid coordss)
  coord_adj <- array_coords + 1
  rows <- coord_adj[, 1]
  cols <- coord_adj[, 2]
  
  # 2. Establish Grid Neighbors
  # Linear Indexing: ID = row + (col - 1) * max_row
  max_row <- max(rows)
  global_ids <- rows + (cols - 1) * max_row
  
  # Define theoretical neighbor IDs
  target_id_bl <- global_ids + 1           # Bottom-Left
  target_id_tr <- global_ids + max_row     # Top-Right
  target_id_br <- global_ids + max_row + 1 # Bottom-Right
  
  # Find indices of these neighbors in the original array_coords list
  idx_tl <- 1:nrow(array_coords)
  idx_bl <- match(target_id_bl, global_ids)
  idx_tr <- match(target_id_tr, global_ids)
  idx_br <- match(target_id_br, global_ids)
  
  # Handle boundary wrapping (bottom edge shouldn't connect to next column top)
  is_bottom_edge <- (rows == max_row)
  idx_bl[is_bottom_edge] <- NA
  idx_br[is_bottom_edge] <- NA
  
  # 3. Build Raw Matrices (using original 1:N indices)
  
  # -- Quads --
  valid_quad <- !is.na(idx_bl) & !is.na(idx_br) & !is.na(idx_tr)
  raw_quad <- cbind(
    idx_tl[valid_quad], 
    idx_bl[valid_quad], 
    idx_br[valid_quad], 
    idx_tr[valid_quad]
  )
  
  # -- Faces (Triangles) --
  # Delta 1: TL-BL-BR
  valid_t1 <- !is.na(idx_bl) & !is.na(idx_br)
  faces_t1 <- cbind(idx_tl[valid_t1], idx_bl[valid_t1], idx_br[valid_t1])
  
  # Delta 2: BR-TR-TL
  valid_t2 <- !is.na(idx_br) & !is.na(idx_tr)
  faces_t2 <- cbind(idx_br[valid_t2], idx_tr[valid_t2], idx_tl[valid_t2])
  
  raw_face <- rbind(faces_t1, faces_t2)
  
  # 3b. Filter faces for quad meshes: only keep faces that are entirely contained in at least one quad
  if (mesh_type == "quad" && nrow(raw_quad) > 0 && nrow(raw_face) > 0) {
    # For each face, check if it is contained in at least one quad
    # A face is contained in a quad if at least 3 of its vertices match the quad's vertices
    flag <- apply(raw_face, 1, function(v) {
      # Check which quad vertices are in the face vertices
      matches <- matrix(raw_quad %in% v, nrow = nrow(raw_quad))
      # Count how many vertices of each quad match the face (should be >= 3 for containment)
      quad_match_counts <- rowSums(matches)
      # Return TRUE if at least one quad contains this face (has >= 3 matching vertices)
      any(quad_match_counts >= 3)
    })
    # Keep only faces that are contained in at least one quad
    raw_face <- raw_face[flag, , drop = FALSE]
    if (nrow(raw_face) == 0) {
      warning("After filtering faces to those contained in quads, no faces remain.")
    }
  }
  
  # 4. Filtering and Remapping
  
  # Determine which spots to keep based on mesh_type
  if (mesh_type == "quad") {
    if (nrow(raw_quad) == 0) stop("No valid quads found.")
    kept_indices <- unique(as.vector(raw_quad))
  } else {
    if (nrow(raw_face) == 0) stop("No valid faces found.")
    kept_indices <- unique(as.vector(raw_face))
  }
  
  # Sort kept indices to maintain relative order (optional but cleaner)
  kept_indices <- sort(kept_indices)
  
  # Messaging
  n_total <- nrow(array_coords)
  n_kept  <- length(kept_indices)
  message(sprintf("MeshBuilder: Filtered out %d spots. Keeping %d spots used in %s mesh.", 
                  n_total - n_kept, n_kept, mesh_type))
  
  # Create a Lookup Map: Old Index -> New Index
  # Initialize with NA, fill kept slots with 1:n_kept
  id_map <- rep(NA, n_total)
  id_map[kept_indices] <- 1:n_kept
  
  # 5. Apply Map to Matrices and Clean Up
  
  # Remap Quads
  new_quad <- matrix(id_map[raw_quad], ncol = 4)
  # Keep only rows where ALL vertices survived the filter
  # (If mesh_type="face", a quad might lose a point that was part of a face but not a quad, though unlikely given geometry)
  new_quad <- new_quad[complete.cases(new_quad), , drop = FALSE]
  
  # Remap Faces
  new_face <- matrix(id_map[raw_face], ncol = 3)
  # Keep only rows where ALL vertices survived
  # (If mesh_type="quad", faces consisting of "loose" triangles not part of a quad will be dropped here)
  new_face <- new_face[complete.cases(new_face), , drop = FALSE]
  
  # 6. Return Output
  return(list(
    face = new_face,
    quad = new_quad,
    kept_indices = kept_indices
  ))
}

# Main entry: construct SpatialMesh and attach to object
FindSpatialNeighbors <- function(object, images = NULL, col.names = NULL, is.visium = TRUE, mesh.type = "face",
                                 positions.path = NULL, verbose = FALSE) {
  barcode <- colnames(object)
  if (!nchar(barcode[1]) >= 18) {
    is.visium <- FALSE
  }

  # 1) Parse coords
  array_coords <- NULL
  # if (is.visium) {
  #   array_coords <- GetArrayCoordsFromImage(object, images, barcode)
  # }
  if (is.null(array_coords)) {
    positions.path <- ResolvePositionsPath(positions.path)
    array_coords <- ParseArrayCoords(barcode, positions.path)
  }
  spatial_coords <- ParseSpatialCoords(object, array_coords, col.names = col.names, images = images, is.visium = is.visium)

  # 2) Initial mesh and filtering
  if (verbose) message("Step 1: Filtering outliers...")
  mesh_init <- MeshBuilder(array_coords, mesh_type = mesh.type)
  idx_keep <- mesh_init$kept_indices
  array_coords <- array_coords[idx_keep, , drop = FALSE]
  spatial_coords <- spatial_coords[idx_keep, , drop = FALSE]
  object <- object[, idx_keep]

  # 3) Fill holes in array coords, then interpolate spatial coords using filled mesh faces (S3 FillHoles)
  if (verbose) message("Step 2: Processing (Hole filling & Interpolation)...")
  class(array_coords) <- "array_coords"
  filled_array_coords <- FillHoles(array_coords, verbose = verbose)

  # build mesh on filled coords to drive vertex interpolation
  mesh_after_fill <- MeshBuilder(filled_array_coords, mesh_type = mesh.type)
  mesh_faces <- if (!is.null(mesh_after_fill$face) && nrow(mesh_after_fill$face) > 0) mesh_after_fill$face else mesh_after_fill$quad

  class(spatial_coords) <- "vertex"
  filled_spatial_coords <- FillHoles(spatial_coords, face = mesh_faces, array_coords = filled_array_coords, verbose = verbose)

  # 4) Final mesh on filled coordss (reuse build)
  if (verbose) message("Step 3: Generating final topology...")
  final_mesh <- mesh_after_fill
  if (length(final_mesh$kept_indices) != nrow(filled_array_coords)) {
    filled_array_coords <- filled_array_coords[final_mesh$kept_indices, , drop = FALSE]
    filled_spatial_coords <- filled_spatial_coords[final_mesh$kept_indices, , drop = FALSE]
  }

  # 5) Boundary and corner computation
  boundary <- numeric()
  corner <- numeric()
  topo_face <- if (!is.null(final_mesh$face)) final_mesh$face else final_mesh$quad
  if (exists("compute_bd")) {
    boundary <- compute_bd(topo_face)
    if (length(boundary) >= 4) {
      corner <- BoundaryBuilder(boundary)
    } else {
      corner <- boundary
    }
  }

  # 6) Prepare metadata holders
  layout <- list(spatial_coords = filled_spatial_coords, array_coords = filled_array_coords, spatial_pseudo = filled_array_coords)
  face_metadata <- data.frame(matrix(data = NA, nrow = nrow(final_mesh$face), ncol = 0))
  vertex_metadata <- data.frame(matrix(data = NA, nrow = nrow(filled_spatial_coords), ncol = 0))

  # Visium shift for array coords
  if (is.visium) {
    layout$spatial_pseudo[, 2] <- layout$spatial_pseudo[, 2] - ceiling(layout$spatial_pseudo[, 1] / 2)
  }

  # 7) Construct SpatialMesh and return object
  object@tools[["spatial.mesh"]] <- new("SpatialMesh",
                                        layout = layout,
                                        face_metadata = face_metadata,
                                        vertex_metadata = vertex_metadata,
                                        face = final_mesh$face,
                                        quad = final_mesh$quad,
                                        boundary = boundary,
                                        corner = corner)

  if (verbose) {
    message(sprintf("Mesh built: %s vertices, %s faces, %s quads",
                    nrow(layout$spatial_coords), nrow(final_mesh$face), nrow(final_mesh$quad)))
  }
  object
}
