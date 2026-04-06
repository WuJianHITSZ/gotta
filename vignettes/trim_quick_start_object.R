# Trim a Seurat object down to the minimal payload needed by
# vignettes/quick_start.Rmd.
#
# The quick-start pipeline only needs:
# - RegionLoupe in meta.data
# - original_y, original_x in meta.data or Staffli@meta.data
# - reduction "rnapca" (first 10 dimensions by default)
# - optional Staffli tool if coordinates are stored there
#
# Everything else is aggressively removed to make a lightweight
# reproducible toy object for GitHub distribution.

suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(SeuratObject)
})


trim_quick_start_object <- function(
  source_rds,
  output_rds,
  keep_reduction = "rnapca",
  keep_dims = 10,
  keep_meta = c("original_y", "original_x", "RegionLoupe", "seurat_clusters"),
  keep_tools = c("Staffli"),
  keep_assay = NULL,
  max_cells = NULL,
  sampling_method = "spatial_crop",
  stratify_by = "RegionLoupe",
  seed = 1L,
  keep_spatial_mesh = FALSE,
  compress = "xz"
) {
  object <- readRDS(source_rds)

  if (!methods::is(object, "Seurat")) {
    stop("source_rds must contain a Seurat object.")
  }

  if (!keep_reduction %in% names(object@reductions)) {
    stop(sprintf("Reduction '%s' not found in object@reductions.", keep_reduction))
  }

  if (!("RegionLoupe" %in% colnames(object@meta.data))) {
    stop("quick_start.Rmd requires 'RegionLoupe' in object@meta.data.")
  }
  if (!has_quick_start_coords(object)) {
    stop(
      "quick_start.Rmd requires original_x/original_y in object@meta.data ",
      "or object@tools[['Staffli']]@meta.data."
    )
  }

  if (is.null(keep_assay)) {
    keep_assay <- SeuratObject::DefaultAssay(object)
  }
  if (!keep_assay %in% names(object@assays)) {
    stop(sprintf("Assay '%s' not found in object@assays.", keep_assay))
  }

  cells_to_keep <- colnames(object)
  if (!is.null(max_cells) && max_cells < length(cells_to_keep)) {
    set.seed(seed)
    cells_to_keep <- choose_cells(
      object = object,
      max_cells = max_cells,
      sampling_method = sampling_method,
      stratify_by = stratify_by
    )
    object <- subset(object, cells = cells_to_keep)
  }

  keep_meta <- unique(c(keep_meta, "original_y", "original_x", "RegionLoupe"))
  keep_meta <- intersect(keep_meta, colnames(object@meta.data))
  object@meta.data <- object@meta.data[, keep_meta, drop = FALSE]

  object <- slim_reduction(
    object = object,
    reduction = keep_reduction,
    keep_dims = keep_dims
  )

  object <- drop_other_reductions(object, keep_reduction = keep_reduction)
  object <- drop_other_assays(object, keep_assay = keep_assay)
  object <- slim_assay_payload(object, assay = keep_assay)
  object <- drop_graphs_neighbors_commands(object)
  object <- drop_images_misc(object)
  object <- trim_tools(
    object,
    keep_tools = keep_tools,
    keep_spatial_mesh = keep_spatial_mesh
  )

  SeuratObject::DefaultAssay(object) <- keep_assay

  saveRDS(object, output_rds, compress = compress)

  invisible(object)
}


choose_cells <- function(object, max_cells, sampling_method = "spatial_crop",
                         stratify_by = "RegionLoupe") {
  all_cells <- colnames(object)

  if (identical(sampling_method, "spatial_crop") && has_quick_start_coords(object)) {
    return(choose_cells_spatial(object, max_cells))
  }

  if (!(stratify_by %in% colnames(object@meta.data))) {
    return(sample(all_cells, size = max_cells, replace = FALSE))
  }

  groups <- as.character(object@meta.data[[stratify_by]])
  groups[is.na(groups)] <- "NA"
  split_cells <- split(all_cells, groups)

  group_sizes <- vapply(split_cells, length, integer(1))
  target_sizes <- stats::setNames(
    pmax(1L, floor(group_sizes / sum(group_sizes) * max_cells)),
    names(group_sizes)
  )
  target_sizes <- stats::setNames(
    pmin(as.numeric(target_sizes), as.numeric(group_sizes)),
    names(group_sizes)
  )

  current_total <- sum(target_sizes)
  if (current_total < max_cells) {
    deficit <- max_cells - current_total
    room <- group_sizes - target_sizes
    while (deficit > 0 && any(room > 0)) {
      for (nm in names(room)[room > 0]) {
        target_sizes[nm] <- target_sizes[nm] + 1L
        room[nm] <- room[nm] - 1L
        deficit <- deficit - 1L
        if (deficit == 0) break
      }
    }
  }

  target_sizes <- target_sizes[names(split_cells)]
  target_sizes[is.na(target_sizes)] <- 1L

  selected <- unlist(
    lapply(seq_along(split_cells), function(i) {
      cells <- split_cells[[i]]
      n_keep <- min(length(cells), as.integer(target_sizes[[i]]))
      sample(cells, size = n_keep, replace = FALSE)
    }),
    use.names = FALSE
  )

  selected
}


choose_cells_spatial <- function(object, max_cells) {
  coords <- get_quick_start_coords(object)
  coords <- coords[colnames(object), , drop = FALSE]

  center <- colMeans(coords, na.rm = TRUE)
  dist2 <- (coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2
  keep_idx <- order(dist2, na.last = NA)[seq_len(min(max_cells, nrow(coords)))]
  rownames(coords)[keep_idx]
}


slim_reduction <- function(object, reduction, keep_dims) {
  red <- object@reductions[[reduction]]
  embed <- red@cell.embeddings

  keep_dims <- min(keep_dims, ncol(embed))
  red@cell.embeddings <- embed[, seq_len(keep_dims), drop = FALSE]

  if ("feature.loadings" %in% slotNames(red)) {
    red@feature.loadings <- matrix(0, nrow = 0, ncol = 0)
  }
  if ("feature.loadings.projected" %in% slotNames(red)) {
    red@feature.loadings.projected <- matrix(0, nrow = 0, ncol = 0)
  }
  if ("jackstraw" %in% slotNames(red)) {
    red@jackstraw <- new("JackStrawData")
  }
  if ("misc" %in% slotNames(red)) {
    red@misc <- list()
  }
  if ("stdev" %in% slotNames(red) && length(red@stdev) >= keep_dims) {
    red@stdev <- red@stdev[seq_len(keep_dims)]
  }

  object@reductions[[reduction]] <- red
  object
}


drop_other_reductions <- function(object, keep_reduction) {
  object@reductions <- object@reductions[keep_reduction]
  object
}


drop_other_assays <- function(object, keep_assay) {
  object@assays <- object@assays[keep_assay]
  object
}


slim_assay_payload <- function(object, assay) {
  assay_obj <- object@assays[[assay]]
  cell_names <- colnames(object)
  feature_names <- rownames(assay_obj)

  if (length(feature_names) == 0) {
    stop(sprintf("Assay '%s' has no features to retain.", assay))
  }

  keep_feature <- feature_names[1]

  if ("counts" %in% slotNames(assay_obj)) {
    counts <- assay_obj@counts[keep_feature, , drop = FALSE]
    counts[] <- 0
    assay_obj@counts <- counts
  }
  if ("data" %in% slotNames(assay_obj)) {
    data <- assay_obj@data[keep_feature, , drop = FALSE]
    data[] <- 0
    assay_obj@data <- data
  }
  if ("scale.data" %in% slotNames(assay_obj)) {
    assay_obj@scale.data <- matrix(
      0,
      nrow = 0,
      ncol = length(cell_names),
      dimnames = list(character(0), cell_names)
    )
  }
  if ("meta.features" %in% slotNames(assay_obj)) {
    assay_obj@meta.features <- assay_obj@meta.features[keep_feature, , drop = FALSE]
  }
  if ("var.features" %in% slotNames(assay_obj)) {
    assay_obj@var.features <- character(0)
  }
  if ("misc" %in% slotNames(assay_obj)) {
    assay_obj@misc <- list()
  }

  object@assays[[assay]] <- assay_obj
  object
}


drop_graphs_neighbors_commands <- function(object) {
  if ("graphs" %in% slotNames(object)) {
    object@graphs <- list()
  }
  if ("neighbors" %in% slotNames(object)) {
    object@neighbors <- list()
  }
  if ("commands" %in% slotNames(object)) {
    object@commands <- list()
  }
  object
}


drop_images_misc <- function(object) {
  if ("images" %in% slotNames(object)) {
    object@images <- list()
  }
  if ("misc" %in% slotNames(object)) {
    object@misc <- list()
  }
  object
}


has_quick_start_coords <- function(object) {
  !is.null(get_quick_start_coords(object))
}


get_quick_start_coords <- function(object) {
  meta_has_coords <- all(c("original_x", "original_y") %in% colnames(object@meta.data))
  if (meta_has_coords) {
    coords <- as.matrix(object@meta.data[, c("original_x", "original_y"), drop = FALSE])
    return(coords)
  }

  if (!("Staffli" %in% names(object@tools))) {
    return(NULL)
  }

  staffli <- object@tools[["Staffli"]]
  if (!methods::is(staffli, "Staffli")) {
    return(NULL)
  }

  if (!all(c("original_x", "original_y") %in% colnames(staffli@meta.data))) {
    return(NULL)
  }

  as.matrix(staffli@meta.data[, c("original_x", "original_y"), drop = FALSE])
}


trim_tools <- function(object, keep_tools = c("Staffli"), keep_spatial_mesh = FALSE) {
  keep_tools <- unique(keep_tools)
  if (isTRUE(keep_spatial_mesh)) {
    keep_tools <- unique(c(keep_tools, "spatial.mesh"))
  }

  current_tools <- names(object@tools)
  keep_tools <- intersect(keep_tools, current_tools)
  object@tools <- object@tools[keep_tools]
  object
}


# Example ---------------------------------------------------------------
# source_rds <- "C:/path/to/full-object.rds"
# output_rds <- "C:/path/to/quick-start-toy.rds"
#
# trim_quick_start_object(
#   source_rds = source_rds,
#   output_rds = output_rds,
#   keep_reduction = "rnapca",
#   keep_dims = 10,
#   max_cells = 2000,
#   stratify_by = "RegionLoupe",
#   keep_tools = c("Staffli"),
#   keep_spatial_mesh = FALSE
# )
