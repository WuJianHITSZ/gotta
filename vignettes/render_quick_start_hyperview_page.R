library(gotta)
library(Matrix) # fill holes
library(patchwork)
library(htmlwidgets)

repo_root <- "C:/Users/jian/Documents/R/gotta-20260404-124151"
toy_path <- file.path(repo_root, "inst", "extdata", "example-object-sma-mouse-heart-3-toy.rds")
docs_dir <- file.path(repo_root, "docs")
page_dir <- file.path(docs_dir, "quick-start-hyperview")
page_path <- file.path(page_dir, "index.html")
page_libs <- file.path(page_dir, "libs")
preview_path <- file.path(repo_root, "man", "figures", "quick-start-hyperview-preview.png")
edge_path <- "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"

dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(preview_path), recursive = TRUE, showWarnings = FALSE)

object <- readRDS(toy_path)

vertex.job <- VertexJob(reduction = "rnapca", n.components = 10)
gott.job <- GOTTJob(vertex.job = vertex.job, shape = "free_boundary")

object <- FindSpatialNeighbors(object, col.names = c("original_y", "original_x"), verbose = TRUE)
object <- ComputeCellArea(object, vertex.job = vertex.job, verbose = TRUE)
object <- RunGOTT(object, gott.job = gott.job, is.pseudo.initial = TRUE, verbose = TRUE)
object <- RunHyperView(object, vertex.job = vertex.job, layout.name = "rnapca.free_boundary")

hyperview_widget <- HyperMeshFeaturePlot(
  object,
  layout.name = "rnapca.free_boundary",
  features = "RegionLoupe"
)

hyperview_widget <- plotly::layout(
  hyperview_widget,
  title = list(text = "gotta Quick Start: HyperView"),
  scene = list(
    xaxis = list(title = "HyperView 1"),
    yaxis = list(title = "HyperView 2"),
    zaxis = list(title = "HyperView 3")
  )
)

htmlwidgets::saveWidget(
  hyperview_widget,
  file = page_path,
  selfcontained = FALSE,
  libdir = page_libs,
  title = "gotta Quick Start HyperView"
)

if (!file.exists(edge_path)) {
  stop("Microsoft Edge not found; cannot render HyperView preview PNG.")
}

html_uri <- paste0("file:///", normalizePath(page_path, winslash = "/"))
cmd <- sprintf(
  '"%s" --headless --disable-gpu --hide-scrollbars --window-size=1400,1000 --screenshot="%s" "%s"',
  edge_path,
  normalizePath(preview_path, winslash = "/", mustWork = FALSE),
  html_uri
)
status <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = FALSE)
if (!identical(status, 0L) || !file.exists(preview_path)) {
  stop("Failed to render HyperView preview PNG.")
}
