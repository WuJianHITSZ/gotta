FillHoles <- function(x, ...) {
  UseMethod("FillHoles")
}

#' Fill missing array coordinates (default method)
FillHoles.array_coords <- function(x, verbose = FALSE, ...) {
  array_coords <- x
  stopifnot(is.matrix(array_coords), ncol(array_coords) == 2)

  if (verbose) {
    message("Step 1: Filling geometric holes...")
    message(sprintf("[FillHoles.array_coords] Input coords: %s rows", nrow(array_coords)))
  }

  image_size <- apply(array_coords, 2, max)
  img <- matrix(FALSE, nrow = image_size[1], ncol = image_size[2])

  index_old <- (array_coords[, 2] - 1L) * nrow(img) + array_coords[, 1]
  img[index_old] <- TRUE

  img_filled <- fillHull(img)

  index_new <- which(img_filled)
  index_holes <- setdiff(index_new, index_old)
  n_added <- length(index_holes)
  if (n_added > 0) {
    holes_cols <- ((index_holes - 1L) %/% nrow(img_filled)) + 1L
    holes_rows <- ((index_holes - 1L) %%  nrow(img_filled)) + 1L
    array_coords <- rbind(array_coords, cbind(holes_rows, holes_cols))
  }

  if (verbose) {
    message(sprintf("[FillHoles.array_coords] Added %s hole pixels; output coords: %s rows", n_added, nrow(array_coords)))
  }

  array_coords
}

#' Fill missing vertex rows using mesh adjacency (vertex method)
FillHoles.vertex <- function(x, face, array_coords, verbose = FALSE, ...) {
  vertex <- x
  nCell <- max(face)
  vertex_new <- matrix(data = mean(vertex), nrow = nCell, ncol = ncol(vertex))
  vertex_new[seq_len(nrow(vertex)), ] <- vertex

  adjacency_matrix <- sparseMatrix(i = c(t(face)), j = c(t(face[, c(2, 3, 1)])), x = TRUE, dims = c(nCell, nCell))
  index_holes <- (nrow(vertex) + 1):nCell

  if (verbose) {
    message(sprintf("[FillHoles.vertex] nCell=%s, observed=%s, holes=%s", nCell, nrow(vertex), length(index_holes)))
  }

  if (verbose) message("Initializing hole values via KNN...")
  index_nn <- knn1(array_coords[seq_len(nrow(vertex)), ], array_coords[index_holes, ], seq_len(nrow(vertex)))
  vertex_new[index_holes, ] <- vertex[index_nn, ]

  # Precompute neighbor structure for holes; sparse multiply to average quickly
  neighbor_weights <- adjacency_matrix[, index_holes, drop = FALSE]
  neighbor_counts <- Matrix::colSums(neighbor_weights)
  neighbor_counts[neighbor_counts == 0] <- 1  # avoid divide-by-zero if any isolated

  if (verbose) message("Building adjacency matrix...")

  if (verbose) message(sprintf("Smoothing holes for %d iterations...", 10))
  for (i in seq_len(10)) {
    averaged <- as.matrix(t(neighbor_weights) %*% vertex_new)
    averaged <- averaged / neighbor_counts
    vertex_new[index_holes, ] <- averaged
  }

  if (verbose) {
    message(sprintf("[FillHoles.vertex] Filled %s missing vertices; final dims: %s x %s",
                    length(index_holes), nrow(vertex_new), ncol(vertex_new)))
  }

  vertex_new
}
