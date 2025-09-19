
library(plotly)

face <- object@tools[["mesh_Spatial"]]@face - 1
layout <- object@reductions[["pca"]]@cell.embeddings[,1:3]

color_map <- c( "green","cyan","yellow","red","purple", "blue")  # Assign colors to factor levels
colors <- color_map[as.numeric(factor(object$RegionLoupe, levels = unique(object$RegionLoupe)))]


# Create 3D polygon (polygon mesh) using plotly
mesh_plot <- plot_ly()

mesh_plot <- add_trace(
  mesh_plot,
  type = "mesh3d",
  x = layout[,1],
  y = layout[,2],
  z = layout[,3],
  i = face[,1],  # Indices of layout for each triangle (0-based)
  j = face[,2],
  k = face[,3],
  color = I("rgba(255, 0, 0, 0.7)"),  # Fill color (red with alpha = 0.7)
  opacity = 0.7,  # Opacity of the polygon surface
  lighting = list(ambient = 0.8)  # Adjust lighting (ambient intensity)
)

mesh_plot <- add_trace(
  mesh_plot,
  type = "scatter3d",
  mode = "markers",
  x = layout[,1],
  y = layout[,2],
  z = layout[,3],
  marker = list(
    size = 5,  # Marker size
    color = colors,  # Color based on z-coordinate
    colorscale = "Viridis",  # Color scale
    opacity = 0.8  # Marker opacity
  )
)

# Display the 3D polygon plot
print(mesh_plot)




# Plot 3D surface using plot_ly() with Delaunay mesh
surface_plot <- plot_ly()

surface_plot <- add_mesh(
  surface_plot,
  x = layout[,1],
  y = layout[,2],
  z = layout[,3],
  i = face[,1],  # Indices of layout for each triangle (0-based)
  j = face[,2],
  k = face[,3],
  intensity = object$nCount_MSI,  # Disable vertex intensity (no color fill)
  colorscale = "Viridis",  # Transparent color
  showscale = TRUE  # Hide color scale legend
)

surface_plot <- surface_plot %>% layout(
  title = "Gaussian Curvature"
)

print(surface_plot)

