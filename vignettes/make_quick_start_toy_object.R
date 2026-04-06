source("C:/Users/jian/Documents/R/gotta-20260404-124151/vignettes/trim_quick_start_object.R")

source_rds <- "C:/Users/jian/Documents/R/result/example-object-sma-mouse-heart-3.rds"
output_rds <- "C:/Users/jian/Documents/R/result/example-object-sma-mouse-heart-3-toy.rds"
package_rds <- "C:/Users/jian/Documents/R/gotta-20260404-124151/inst/extdata/example-object-sma-mouse-heart-3-toy.rds"

trim_quick_start_object(
  source_rds = source_rds,
  output_rds = output_rds,
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
)

dir.create(dirname(package_rds), recursive = TRUE, showWarnings = FALSE)
file.copy(output_rds, package_rds, overwrite = TRUE)

message("Saved toy object to: ", output_rds)
message("Copied toy object to: ", package_rds)
