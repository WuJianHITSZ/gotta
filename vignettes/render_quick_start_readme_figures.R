library(gotta)
library(Matrix) # fill holes
library(patchwork)
library(ggplot2)

toy_path <- normalizePath(
  "C:/Users/jian/Documents/R/gotta-20260404-124151/inst/extdata/example-object-sma-mouse-heart-3-toy.rds",
  mustWork = TRUE
)
fig_dir <- normalizePath(
  "C:/Users/jian/Documents/R/gotta-20260404-124151/man/figures",
  mustWork = FALSE
)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(toy_path)

vertex.job <- VertexJob(reduction = "rnapca", n.components = 10)
gott.job <- GOTTJob(vertex.job = vertex.job, shape = "free_boundary")

object <- FindSpatialNeighbors(object, col.names = c("original_y", "original_x"), verbose = TRUE)
object <- ComputeCellArea(object, vertex.job = vertex.job, verbose = TRUE)

spatial_plot <- MeshFeaturePlot(object, layout.name = "spatial_coords", features = "RegionLoupe") +
  SurfFeaturePlot(object, layout.name = "spatial_coords", features = "area.rnapca")
ggsave(
  filename = file.path(fig_dir, "quick-start-spatial-combined.png"),
  plot = spatial_plot,
  width = 12,
  height = 6,
  dpi = 150
)

object <- RunGOTT(object, gott.job = gott.job, is.pseudo.initial = TRUE, verbose = TRUE)

gott_plot <- MeshFeaturePlot(object, layout.name = "rnapca.free_boundary", features = "RegionLoupe") +
  SurfFeaturePlot(object, layout.name = "rnapca.free_boundary")
ggsave(
  filename = file.path(fig_dir, "quick-start-gott-combined.png"),
  plot = gott_plot,
  width = 12,
  height = 6,
  dpi = 150
)
