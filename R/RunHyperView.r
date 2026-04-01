RunHyperView <- function(object, gott.job=NULL, vertex.job=NULL, layout.name=NULL, alpha = 0.99){
    if(!is.null(gott.job)){
        vertex.job <- gott.job$vertex.job
        shape <- gott.job$shape
    }
    if(!is.null(vertex.job)){
        assay <- vertex.job$assay
        features <- vertex.job$features
        reduction <- vertex.job$reduction
        n.components <- vertex.job$n.components
    }
    if(!is.null(reduction)){
        assay <- reduction
    }
    if(is.null(layout.name)){
        layout.name <- paste0(assay,".",shape)
    }
    mesh <- SpatialMesh(object)
    face <- mesh@face

    spatial_coords <- mesh@layout$spatial_coords
    
    vertex <- VertexAccessor(object, assay, features, reduction, n.components)
    if(nrow(spatial_coords)>nrow(vertex)){
        class(vertex) <- "vertex"
        vertex <- FillHoles(vertex, face = face, array_coords = spatial_coords)
    }


    uv <- mesh@layout[[layout.name]]

    population <- face_area(face,vertex) / 3

    mixed_data_opt <- cbind(uv * alpha, vertex * (1 - alpha))
    
    # 2. PCA
    pca_res_opt <- prcomp(mixed_data_opt, center = TRUE, rank. = 3)
    uvw <- pca_res_opt$x[, 1:3]
    
    # 3. Calculate Scale Factor
    # Note: In the loop, we calculate area, then scale, then calculate again.
    # The MATLAB script repeats this logic outside the loop.
    current_areas_opt <- face_area(face, uvw)
    scale_factor <- sum(population) / sum(current_areas_opt)
    
    # 4. Apply Scaling
    uvw <- uvw * sqrt(scale_factor)

    # 5. Store Results
    hyperview.name <- paste0(layout.name, ".hyperview")

    mesh@layout[[hyperview.name]] <- uvw
    mesh@vertex_metadata[[paste0("area.",hyperview.name)]] <- vertex_area(face, uvw) / 3
    mesh@face_metadata[[paste0("area.",hyperview.name)]] <- face_area(face, uvw)
    mesh@vertex_metadata[[paste0("curvature.",layout.name)]] <- discrete_gaussian_curvature(uvw, face)
    
    SpatialMesh(object) <- mesh

    return(object)
}
