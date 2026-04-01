# SpatialMesh S4 class definition
setClass("SpatialMesh", slots = list(layout = "list", 
                                face_metadata = "data.frame", 
                                vertex_metadata = "data.frame", 
                                face = "matrix", 
                                quad = "matrix", 
                                boundary = "integer", 
                                corner = "integer", 
                                cartogram = "list"),
         prototype = list(layout = list(),
                         face_metadata = data.frame(),
                         vertex_metadata = data.frame(),
                         face = matrix(nrow = 0, ncol = 0),
                         quad = matrix(nrow = 0, ncol = 0),
                         boundary = integer(0),
                         corner = integer(0),
                         cartogram = list()))

# VertexJob S3 constructor: represents vertex processing job
# Based on SimpleVertexFactory which returns: list(vertex = matrix, area = numeric vector)
# The vertex matrix is stored in reduction, and area vector is stored in vertex_metadata
VertexJob <- function(assay = NULL, reduction = NULL, features = NULL, 
                     n.components = NULL, face_area = numeric(0), 
                     vertex_area = numeric(0), vertex = matrix(nrow = 0, ncol = 0)) {
  # If reduction is not null, then assay = reduction
  if (!is.null(reduction)) {
    assay <- reduction
  }
                      
  structure(list(
    assay = assay,
    reduction = reduction,
    features = features,
    n.components = n.components,
    face_area = face_area,
    vertex_area = vertex_area,
    vertex = vertex
  ), class = "VertexJob")
}

# GOTTJob S3 constructor: represents GOTT processing job
# Uses composition: contains a VertexJob as a field
GOTTJob <- function(vertex.job = NULL, assay = NULL, reduction = NULL, 
                   features = NULL, n.components = NULL, face_area = numeric(0), 
                   vertex_area = numeric(0), vertex = matrix(nrow = 0, ncol = 0),
                   shape = "free_boundary", aspect.ratio = 1.0, 
                   layout.name = character(0)) {
  # Create or use provided VertexJob
  if (is.null(vertex.job)) {
    vertex.job <- VertexJob(assay = assay, reduction = reduction, features = features,
                           n.components = n.components, face_area = face_area,
                           vertex_area = vertex_area, vertex = vertex)
  }
  # Create GOTTJob structure with VertexJob as a field
  job <- structure(list(
    vertex.job = vertex.job,
    shape = shape,
    aspect.ratio = aspect.ratio,
    layout.name = layout.name
  ), class = "GOTTJob")
  return(job)
}

# AlignmentJob S3 constructor: represents alignment processing job
# Used for RunGOTTA to align source and target assays using optimal transport
AlignmentJob <- function(shape = "disk",
                        source.assay = NULL,
                        target.assay = NULL,
                        target.object = NULL,
                        target.name = NULL,
                        is.spatial.initial = FALSE,
                        is.hyper.target = FALSE) {
  # If target.assay is NULL, set it to source.assay
  if (is.null(target.assay) && !is.null(source.assay)) {
    target.assay <- source.assay
  }
  
  structure(list(
    shape = shape,
    source.assay = source.assay,
    target.assay = target.assay,
    target.object = target.object,
    target.name = target.name,
    is.spatial.initial = is.spatial.initial,
    is.hyper.target = is.hyper.target
  ), class = "AlignmentJob")
}

