# gotta

**Geometric Optimal Transport Tableau & Alignment**

`gotta` is an R package for spatial transcriptomics analysis using geometric optimal transport methods. It provides tools for constructing spatial meshes, computing optimal transport cartograms, aligning assay layouts, estimating trajectories, and visualizing mesh-based spatial analyses.

## Features

- **GOTT (Geometric Optimal Transport Tableau)**: Transform spatial coordinates into cartogram layouts that preserve geometric properties while revealing structure in high-dimensional data
- **GOTTA (GOTT Alignment)**: Align different assay layouts (e.g., spatial to PCA) using optimal transport
- **Spatial Mesh Construction**: Build topological mesh representations of spatial transcriptomics data
- **Trajectory Analysis**: Estimate and visualize cellular trajectories in spatial context
- **HyperView Analysis**: Compute and visualize curvature and other geometric properties
- **Rich Visualization**: Multiple plotting functions for mesh-based and cartogram visualizations

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("WuJianHITSZ/gotta")

# Or install from tarball
install.packages("gotta_1.6.9.tar.gz", repos = NULL, type = "source")
```

### System Requirements

- R >= 4.1.0
- C++17 compiler
- System dependencies for RcppEigen

## Quick Start

```r
library(gotta)
library(Seurat)

# Load your spatial transcriptomics object
object <- your_seurat_object

# Find spatial neighbors
object <- FindSpatialNeighbors(object, col.names = c("y", "x"))

# Run GOTT to create cartogram layout
object <- RunGOTT(object, shape = "disk", reduction = "pca", is.pseudo.initial = TRUE)

# Visualize the cartogram
CartFeaturePlot(object, layout.name = "pca.disk")

# Run GOTTA to align spatial and PCA layouts
alignment.job <- AlignmentJob(shape = "disk", source.assay = "spatial", target.assay = "pca")
object <- RunGOTTA(object, alignment.job = alignment.job)

# Visualize trajectory
object <- RunTrajectory(object, target.assay = "pca", source.assay = "spatial", shape = "disk")
ArrowFeaturePlot(object, features = "area.pca.disk")
```

## Main Functions

### Core Analysis

| Function | Description |
|----------|-------------|
| `RunGOTT()` | Run Geometric Optimal Transport Tableau transformation |
| `RunGOTTA()` | Run GOTT Alignment between assays |
| `RunTrajectory()` | Estimate and compute trajectories |
| `RunHyperView()` | Compute curvature and hyper-surface properties |
| `ComputeCellArea()` | Compute cell/vertex areas |

### Spatial Mesh

| Function | Description |
|----------|-------------|
| `SpatialMesh()` | Create or extract spatial mesh from object |
| `FindSpatialNeighbors()` | Find spatial neighbors for mesh construction |
| `FindSpatialCorners()` | Identify corner vertices in spatial layout |
| `FillHoles()` | Fill holes in spatial data |

### Visualization

| Function | Description |
|----------|-------------|
| `MeshFeaturePlot()` | Plot features on spatial mesh |
| `CartFeaturePlot()` | Plot features on cartogram layout |
| `SurfFeaturePlot()` | Surface plot of features |
| `ArrowFeaturePlot()` | Arrow plot showing trajectories/deformations |
| `HyperMeshFeaturePlot()` | Hyper-surface mesh visualization |
| `HyperSurfFeaturePlot()` | Hyper-surface feature visualization |

### Job Classes

| Class | Description |
|-------|-------------|
| `VertexJob` | Job configuration for vertex-based operations |
| `GOTTJob` | Job configuration for GOTT transformation |
| `AlignmentJob` | Job configuration for alignment operations |

## Supported Shapes

- **Disk**: Circular boundary constraint for cartogram transformation
- **Rectangle**: Rectangular boundary constraint
- **Free boundary**: Unconstrained transformation

## Vignettes

The package includes several demonstration vignettes:

- `demo_GOTT_disk` - GOTT transformation with disk shape
- `demo_GOTT_rect` - GOTT transformation with rectangle shape
- `demo_GOTT_free_boundary` - GOTT with free boundary
- `demo_GOTTA_disk` - GOTTA alignment with disk shape
- `demo_GOTTA_rect` - GOTTA alignment with rectangle shape
- `demo_GOTT_disk_findtrajectory` - Trajectory analysis workflow

View vignettes with:
```r
browseVignettes("gotta")
```

## Dependencies

### Required
- R (>= 4.1.0)
- Rcpp, RcppEigen
- ggplot2, plotly
- Matrix, geometry
- RTriangle, Rlinsolve
- SeuratObject
- shiny

### Suggested
- Seurat
- STutility
- knitr, rmarkdown

## Architecture

The package uses C++ for performance-critical computations including:
- Power diagram construction
- Harmonic mapping
- Gradient and Hessian calculations
- Discrete Gaussian curvature
- Scattered interpolation

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use `gotta` in your research, please cite:

```
Wu, J. (2026). gotta: Geometric Optimal Transport Tableau & Alignment. 
R package version 1.6.9. https://github.com/WuJianHITSZ/gotta
```

## Issues & Contributions

Please report issues and submit pull requests on the [GitHub repository](https://github.com/WuJianHITSZ/gotta).

## Author

**Jian Wu**  
Email: jianwu@example.com
