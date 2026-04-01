FindSpatialCorners <- function(object, ...) {
  UseMethod("FindSpatialCorners")
}

FindSpatialCorners.SpatialMesh <- function(object, shape = NULL, feature = NULL, gott.job = NULL) {
  mesh <- object
  if (!is.null(gott.job) && is.null(shape)) {
    shape <- gott.job$shape
  }

  if (identical(shape, "free_boundary")) {
    warning("free_boundary shape: no corners need to be selected; returning object unchanged.")
    return(mesh)
  }

  corner_index <- mesh@corner
  corner_num <- if (identical(shape, "rect")) 4 else 1

  if (is.null(feature)) {
    if (ncol(mesh@vertex_metadata) == 0) {
      stop("No vertex_metadata available to use as a plotting feature.")
    }
    feature <- colnames(mesh@vertex_metadata)[1]
  }

  plot_Spatial <- MeshFeaturePlot(mesh, features = feature)
  SelectCorners(plot_Spatial, corner_num) # Produces data.frame corner in parent env

  spatial_coords <- mesh@layout[["spatial_coords"]]
  boundary_idx <- mesh@boundary
  if (length(boundary_idx) == 0L) {
    stop("No boundary information found in spatial.mesh@boundary.")
  }

  boundary_coords <- spatial_coords[boundary_idx, , drop = FALSE]
  corner <- as.matrix(corner)
  boundary_hits <- as.integer(knn(boundary_coords, corner, seq_len(nrow(boundary_coords))))
  corner_index <- boundary_idx[boundary_hits]

  mesh@corner <- corner_index
  mesh
}

FindSpatialCorners.Seurat <- function(object, shape = NULL, feature = "seurat_clusters", gott.job = NULL) {
  mesh <- SpatialMesh(object)
  mesh <- FindSpatialCorners(mesh, shape = shape, feature = feature, gott.job = gott.job)
  SpatialMesh(object) <- mesh
  object
}

# Define a function containing a Shiny app
SelectCorners <- function(plot_Spatial,corner_num) {
  # Define UI
  ui <- shiny::fluidPage(
    if(corner_num==1){
      shiny::titlePanel("Select one point on boundary!")
    }else if(corner_num==4){
      shiny::titlePanel("Select four boundary points in counter-clockwise order!")
    },
    
    shiny::mainPanel(
      plotly::plotlyOutput("plot")
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    corners <- shiny::reactiveValues(count = 0, x = numeric(), y = numeric())
    
    output$plot <- plotly::renderPlotly({
      # Convert ggplot object to plotly
      plotly::ggplotly(plot_Spatial)
    })
    
    shiny::observeEvent(plotly::event_data("plotly_click"), {
      # Extract corners spatial_coordss
      click <- plotly::event_data("plotly_click")
      if (!is.null(click)) {
        corners$count <- corners$count + 1
        
        # Store spatial_coordss
        corners$x[corners$count] <- click$x
        corners$y[corners$count] <- click$y
        
        if (corners$count == corner_num) {
          # Output corner spatial_coordss to a higher enviroment, notice x,y have been interchanged (and a sign)
          assign("corner", data.frame(x = -corners$y, y = corners$x), pos = 1)
          shiny::stopApp()
        }
        
        # Re-render plot to include newly clicked point
        output$plot <- plotly::renderPlotly({
          # Create ggplot object
          plot_marked <- plot_Spatial +
            geom_point(data = data.frame(x = corners$x, y = corners$y), aes(x, -y), color = "black", size = 3)
          
          # Convert ggplot object to plotly
          plotly::ggplotly(plot_marked)
        })
      }
    })
  }
  
  # Run the Shiny application
  app <- shiny::shinyApp(ui = ui, server = server)
  shiny::runApp(app)
}
