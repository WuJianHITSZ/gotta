MeshFeaturePlot <- function(object, features, reduction = "coordinate", title = NULL, pt.size = 2, show.quad = FALSE){
  
  if(show.quad){
    face <- object@graphs[["mesh_Spatial"]]@quad
  }else{
    face <- object@graphs[["mesh_Spatial"]]@face
  }

  layout <- object@reductions[[reduction]]@cell.embeddings
  color <- object@meta.data[[features]]
  
  plot <- triplot(face = face, layout = layout) +
    geom_point(aes(x = layout[,2], y = layout[,1], color = color), size = pt.size) +
    labs(title = title,
         x = paste(reduction,"2",sep = "_"),
         y = paste(reduction,"1",sep = "_"),
         color = features)
  
  if(is.numeric(color[1])){
    plot <- plot + scale_color_viridis_c(option = "turbo")
  }
  
  return(plot)
}

SurfFeaturePlot <- function(object, features, reduction = "coordinate", title = NULL, pt.size = 2, show.quad = FALSE){
  
  if(show.quad){
    face <- object@graphs[["mesh_Spatial"]]@quad
  }else{
    face <- object@graphs[["mesh_Spatial"]]@face
  }
  
  layout <- object@reductions[[reduction]]@cell.embeddings
  color <- object@meta.data[[features]]
  fill <- rowMeans(matrix(data = color[face], nrow = nrow(face), ncol = ncol(face)))
  
  plot <- trisurf(face = face, layout = layout, fill = fill) +
    labs(title = title,
         x = paste(reduction,"2",sep = "_"),
         y = paste(reduction,"1",sep = "_"),
         fill = features)
  
  return(plot)
}

triplot <- function(face, layout){
  
  data_frame_face <- data.frame(vertex = face)
  col_names <- colnames(data_frame_face)
  data_frame_face$group <- rownames(data_frame_face)
  
  data_frame_face <- pivot_longer(data_frame_face, cols = col_names, names_to = NULL, values_to = "vertex")
  
  plot <- ggplot() +
    geom_polygon(data = data_frame_face, aes(x = layout[vertex,2], y = layout[vertex,1], group = group), fill = NA, color="black") +
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  return(plot)
}

trisurf <- function(face, layout, fill){
  
  data_frame_face <- data.frame(vertex = face)
  col_names <- colnames(data_frame_face)
  data_frame_face$group <- rownames(data_frame_face)
  data_frame_face$fill <- fill
  
  data_frame_face <- pivot_longer(data_frame_face, cols = col_names, names_to = NULL, values_to = "vertex")
  
  plot <- ggplot() +
    geom_polygon(data = data_frame_face, aes(x = layout[vertex,2], y = layout[vertex,1], group = group, fill = fill), color="black") +
    scale_fill_viridis_c(option = "turbo") +  
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  return(plot)
}

FindFaceInTissue <- function(face, tissue_spot_index, row){
  
  face_spot_in_tissue <- matrix(face %in% tissue_spot_index, nrow = nrow(face), ncol = ncol(face))
  face <- face[apply(face_spot_in_tissue, 1, sum)==ncol(face),]
  
  face_row <- matrix(row[face], nrow = nrow(face), ncol = ncol(face))
  face_row_range <- apply(face_row, 1, max) - apply(face_row, 1, min)
  face <- face[face_row_range==1,]
  
  return(face)
}

RowCol2Face <- function(row, col){
  
  row <- row + 1
  col <- col + 1
  nrow <- max(row)
  ncol <- max(col)
  tissue_spot_index <- (col - 1) * nrow + row
  slice_spot_index <- 1:(nrow*ncol)
  slice_spot_row <- matrix(rep(1:nrow, each=ncol), nrow = nrow, byrow = TRUE)
  slice_spot_row <- c(slice_spot_row)
  
  ## face index in slice spots
  quad <- cbind(slice_spot_index, slice_spot_index+1, slice_spot_index+nrow+1,  slice_spot_index+nrow);
  delta1 <- FindFaceInTissue(face = quad[,1:3], tissue_spot_index, slice_spot_row)
  delta2 <- FindFaceInTissue(face = quad[,c(3,4,1)], tissue_spot_index, slice_spot_row)
  face <- rbind(delta1,delta2)
  
  ## face index in tissue spots
  face <- matrix(match(face, tissue_spot_index), nrow = nrow(face), ncol = ncol(face))
  
  return(face)
}

RowCol2Quad <- function(row, col){
  
  row <- row + 1
  col <- col + 1
  nrow <- max(row)
  ncol <- max(col)
  tissue_spot_index <- (col - 1) * nrow + row
  slice_spot_index <- 1:(nrow*ncol)
  slice_spot_row <- matrix(rep(1:nrow, each=ncol), nrow = nrow, byrow = TRUE)
  slice_spot_row <- c(slice_spot_row)
  
  ## face index in slice spots
  quad <- cbind(slice_spot_index, slice_spot_index+1, slice_spot_index+nrow+1,  slice_spot_index+nrow);
  quad <- FindFaceInTissue(face = quad, tissue_spot_index, slice_spot_row)
  
  ## face index in tissue spots
  quad <- matrix(match(quad, tissue_spot_index), nrow = nrow(quad), ncol = ncol(quad))
  
  return(quad)
}


FindSpatialNeighbors <- function(object, images = NULL){
  if(!is.null(images)){
    coordinate <- as.matrix(object@images[[images]]@coordinates[,c("imagerow","imagecol")])
    colnames(coordinate) <- c("coordinate_1","coordinate_2")
    object@reductions[["coordinate"]] <- CreateDimReducObject(
      embeddings = coordinate,
      assay = "Spatial",
      key = "coordinate_",
    )    
    row <- object@images[[images]]@coordinates$row
    col <- object@images[[images]]@coordinates$col
    col <- floor(col/2) + ceiling(row/2)
  }else{
    coordinate <- object@reductions[["coordinate"]]@cell.embeddings
    print("you are here")
    row <- coordinate[,1]
    col <- coordinate[,2]
  }
  
  face <- RowCol2Face(row = row, col = col)
  quad <- RowCol2Quad(row = row, col = col)
  
  ## kick out isolated points
  vertex_reduc <- sort(unique(c(face)))
  face <- matrix(match(face, vertex_reduc), nrow = nrow(face), ncol = ncol(face))
  quad <- matrix(match(quad, vertex_reduc), nrow = nrow(quad), ncol = ncol(quad))
  object <- object[,vertex_reduc]
  
  boundary <- ComputeBoundary(face = face, ncell = ncol(object))
  
  setClass("Mesh", slots = list(face = "matrix", boundary = "numeric", corners = "numeric", quad = "matrix"))
  object@graphs[["mesh_Spatial"]] <- new("Mesh", face = face, boundary = boundary, corners = boundary[1], quad = quad)    
 
  return(object)
}


ComputeBoundary <- function(face, ncell){
  adjacency_matrix <- sparseMatrix(i = c(t(face)), j = c(t(face[,c(2,3,1)])), x = 1, dims = c(ncell, ncell))
  
  boundary_list <- which((adjacency_matrix - t(adjacency_matrix))>0, arr.ind = TRUE)
  
  boundary <- vector(mode = "numeric", length = nrow(boundary_list))
  pointer <- boundary
  boundary[1] <- boundary_list[1,1]
  for(count in 1:length(boundary_list)){
    pointer[count] <- which(boundary_list[,1]==boundary[count])
    boundary[count+1] <- boundary_list[pointer[count],2]
    if(boundary[count+1] %in% boundary[c(1:count)]){
      break
    }
  }
  boundary <- boundary[c(1:count+1)]
  
  return(boundary)
}

OutputVertex <- function(object, assay = NULL, features = NULL, reduction.use = "coordinate", n.components = NULL){
  
  if(!is.null(assay)){
    # in newer versions of seurat, counts is stored in layers
    if("layers" %in% slotNames(object@assays[[assay]])){
      vertex <- t(as.matrix(object@assays[[assay]]@layers[["counts"]]))
      print("we are here")
    }else{
      vertex <- t(as.matrix(object@assays[[assay]]@counts))
    }
    
    if(!is.null(features)){
      DefaultAssay(object) <- assay
      features_assay <- rownames(object)
      flag <- (features_assay %in% features)
      if(sum(flag)<2){
        stop("Valid features cannot be less than two!")
      }
      vertex <- vertex[,flag]
      
      flag <- (features %in% features_assay)
      if(sum(flag)>0){
        warnings_text <- paste(sum(flag),"features are not found in the assay",assay)
        warnings(warnings_text)
      }
    }
  }else{
    vertex <- as.matrix(object@reductions[[reduction.use]]@cell.embeddings)
    if(!is.null(n.components)){
      if(n.components<2){
        stop("Valid features cannot be less than two!")
      }
      n.components <- min(n.components,ncol(vertex))
      vertex <- vertex[,c(1:n.components)]
      if(n.components>ncol(vertex)){
        warnings_text <- paste("n.components can not be larger than",ncol(vertex))
        warnings(warnings_text)
      }      
    }
  }
  
  return(vertex)
}

# Define a function containing a Shiny app
SelectCorners <- function(plot_spatial, corner_num) {
  # Define UI
  ui <- fluidPage(
    if(corner_num==1){
      titlePanel("Select one point on boundary!")
    }else if(corner_num==4){
      titlePanel("Select four boundary points in counter-clockwise order!")
    },
    
    mainPanel(
      plotlyOutput("plot")
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    corners <- reactiveValues(count = 0, x = numeric(), y = numeric())
    
    output$plot <- renderPlotly({
      # Convert ggplot object to plotly
      ggplotly(plot_spatial)
    })
    
    observeEvent(event_data("plotly_click"), {
      # Extract corners coordinates
      corner <- event_data("plotly_click")
      if (!is.null(corner)) {
        corners$count <- corners$count + 1
        
        # Store coordinates
        corners$x[corners$count] <- corner$x
        corners$y[corners$count] <- corner$y
        
        if (corners$count == corner_num) {
          # Output corner coordinates to a higher enviroment, notice x,y have been interchanged (and a sign)
          assign("corners", data.frame(x = -corners$y, y = corners$x), pos = 1)
          stopApp()
        }
        
        # Re-render plot to include newly clicked point
        output$plot <- renderPlotly({
          # Create ggplot object
          plot_marked <- plot_spatial +
            geom_point(data = data.frame(x = corners$x, y = corners$y), aes(x, -y), color = "black", size = 3)
          
          # Convert ggplot object to plotly
          ggplotly(plot_marked)
        })
      }
    })
  }
  
  # Run the Shiny application
  app <- shinyApp(ui = ui, server = server)
  runApp(app)
}

RunGOTT <- function(object, sample_id = NULL, reduction.save = "GOTT", shape = "disk", aspect_ratio = 1, assay = NULL, features = NULL, reduction.use = "coordinate", dim.use = NULL, select_points = FALSE){
  # judge whether selecting corners is needed
  corner_num <- 1
  if(shape=="rect"){
    corner_num <- 4
  }
  corner_index <- object@graphs[["mesh_Spatial"]]@corners
  if(length(corner_index)!=corner_num){
    select_points <- TRUE
  }
  
  # select corners or fixed point from a spatial plot via a shiny app
  if(select_points){
    plot_spatial <- MeshFeaturePlot(object, features = "nCount_Spatial")
    
    SelectCorners(plot_spatial,corner_num) # a data frame called corners has been produced
    
    # find index of corners/fixed point via knn classification
    coordinate <- as.matrix(object@reductions$coordinate@cell.embeddings)
    corners <- as.matrix(corners)
    corner_index <- as.integer(knn(coordinate,corners,c(1:nrow(coordinate))))
    object@graphs[["mesh_Spatial"]]@corners <- corner_index
  }  
  
  # write varibale to matlab: barcode, vertex, sample_id, and corner_index
  vertex <- OutputVertex(object, assay = assay, features = features, reduction.use = reduction.use, n.components = dim.use)
  barcode <- colnames(object)
  
  print("writing matlab files...")
  writeMat("gott_input.mat", vertex = vertex, barcode = barcode, sample_id = sample_id, index = corner_index, shape = shape, aspect_ratio = aspect_ratio)
  
  ## run matlab script here, most time-costly!!
  print("trying to boost matlab...")
  system("/Applications/MATLAB_R2022a.app/bin/matlab -nodisplay -r \"run('interface_matlab_r_gott.m'); exit\"")
  
  ## read uv from matlab and write it to seurat object reduction
  print("reading results from matlab...")
  gott_output <- readMat("gott_output.mat")
  
  object@meta.data$area <- gott_output[["area"]]
  uv <- gott_output[["uv"]]
  
  colnames(uv) <- c(paste(reduction.save,"1",sep = "_"),paste(reduction.save,"2",sep = "_"))
  rownames(uv) <- colnames(object)
  
  object@reductions[[reduction.save]] <- CreateDimReducObject(
    embeddings = uv,
    assay = "Spatial",
    key = paste(reduction.save,"",sep = "_"),
  )
  
  return(object)
}

GetCoordinateBarcode <- function(coordinate){
  slice_coordinate_list <- read.csv("tissue_positions_list.csv", header = FALSE)
  slice_coordinate <- slice_coordinate_list[,c("V3","V4")]
  colnames(slice_coordinate) <- c("row","col")
  slice_coordinate$col <- floor(slice_coordinate$col/2)
  
  tissue_spot_index <- knn(slice_coordinate,coordinate,c(1:nrow(slice_coordinate)))
  barcode <- slice_coordinate_list$V1[tissue_spot_index]
  
  return(barcode)
}

CreateSpatialObject <- function(count_matrix, coordinate){
  coordinate <- round(coordinate)
  # colnames(coordinate) <- c("col","row")
  coordinate <- coordinate[,c(2,1)]

  barcode <- GetCoordinateBarcode(coordinate)
  row.names(coordinate) <- barcode
  row.names(count_matrix) <- barcode

  order <- match(sort(barcode),barcode)
  coordinate <- coordinate[order,]
  count_matrix <- count_matrix[order,]

  object <- CreateSeuratObject(counts = as.sparse(t(count_matrix)), assay = "Spatial")

  coordinate <- as.matrix(coordinate)
  colnames(coordinate) <- c("coordinate_1","coordinate_2")
  object@reductions[["coordinate"]] <- CreateDimReducObject(
    embeddings = coordinate,
    assay = "Spatial",
    key = "coordinate_",
  )

  return(object)
}