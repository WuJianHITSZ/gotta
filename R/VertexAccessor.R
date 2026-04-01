VertexAccessor <- function(object, assay = DefaultAssay(object), features = NULL,
                           reduction = NULL, n.components = NULL) {
  if (!is.null(reduction)) {
    if (!reduction %in% names(object@reductions)) {
      stop(sprintf("Reduction '%s' not found in object@reductions.", reduction))
    }
    vertex <- as.matrix(object@reductions[[reduction]]@cell.embeddings)
    max_comp <- ncol(vertex)
    if (!is.null(n.components)) {
      if (n.components < 2) {
        stop("n.components must be at least 2.")
      }
      use_comp <- min(n.components, max_comp)
      vertex <- vertex[, seq_len(use_comp), drop = FALSE]
      if (n.components > max_comp) {
        warning(sprintf("Requested n.components=%s exceeds available=%s; using %s.", n.components, max_comp, max_comp))
      }
    }
  } else {
    if (!assay %in% names(object@assays)) {
      stop(sprintf("Assay '%s' not found in object@assays.", assay))
    }
    # Ensure the desired assay is active for feature retrieval
    DefaultAssay(object) <- assay
    vertex <- t(as.matrix(object@assays[[assay]]@counts))

    if (!is.null(features)) {
      available <- rownames(object)
      present <- features %in% available
      if (sum(present) < 2) {
        stop("At least two requested features must exist in the assay.")
      }
      vertex <- vertex[, features[present], drop = FALSE]
      if (any(!present)) {
        warning(sprintf("%s features not found in assay '%s'.", sum(!present), assay))
      }
    }
  }

  # Pad to 3D if needed for downstream geometry routines
  if (ncol(vertex) == 2) {
    vertex <- cbind(vertex, rep(0, nrow(vertex)))
    colnames(vertex)[3] <- "pad_0"
  }

  storage.mode(vertex) <- "double"
  vertex
}
