## changeblog v1.1.1
# this is a script specilized in geometric computation
# try to compute face areas and vertex areas, the part of calling matlab to run gott is temporally commented

ComputeFaceArea <- function(face, vertex){
  # Assume: 
  #   face   is an n×3 integer matrix (each row = vertex indices of a triangle)
  #   vertex is an m×3 numeric matrix (rows = xyz coords of vertices)
  
  fi <- face[,1]
  fj <- face[,2]
  fk <- face[,3]
  
  vij <- vertex[fj, , drop=FALSE] - vertex[fi, , drop=FALSE]
  vjk <- vertex[fk, , drop=FALSE] - vertex[fj, , drop=FALSE]
  vki <- vertex[fi, , drop=FALSE] - vertex[fk, , drop=FALSE]
  
  # row-wise dot products
  a <- sqrt(rowSums(vij * vij))
  b <- sqrt(rowSums(vjk * vjk))
  c <- sqrt(rowSums(vki * vki))
  
  s <- (a + b + c) / 2.0
  face_area <- sqrt(s * (s - a) * (s - b) * (s - c))   # area of each triangle
  
  return(face_area)
}

# ---- main API ---------------------------------------------------------------

ComputeVertexArea <- function(face, vertex, type = c("one_ring", "mixed")) {
  type <- match.arg(type)
  fa <- ComputeFaceArea(face, vertex)  # per-face area (length = nrow(face))
  
  if (type == "one_ring") {
    he_dat <- compute_halfedge(face)           # he: 3m x 2, heif: 3m
    vals  <- fa[he_dat$heif]                   # face area of each halfedge's face
    idx   <- he_dat$he[, 1]                    # source vertex of each halfedge
    # accumulate (like MATLAB accumarray)
    va <- numeric(nrow(vertex))
    s  <- tapply(vals, idx, sum)
    va[as.integer(names(s))] <- s
    return(va)
  }
  
  # type == "mixed"
  i1 <- face[, 1]; i2 <- face[, 2]; i3 <- face[, 3]
  dvf12 <- vertex[i2, , drop = FALSE] - vertex[i1, , drop = FALSE]
  dvf23 <- vertex[i3, , drop = FALSE] - vertex[i2, , drop = FALSE]
  dvf31 <- vertex[i1, , drop = FALSE] - vertex[i3, , drop = FALSE]
  
  c1 <- vector_cot(dvf12, -dvf31)
  c2 <- vector_cot(dvf23, -dvf12)
  c3 <- vector_cot(dvf31, -dvf23)
  
  # row-wise squared lengths
  l2_12 <- rowSums(dvf12 * dvf12)
  l2_23 <- rowSums(dvf23 * dvf23)
  l2_31 <- rowSums(dvf31 * dvf31)
  
  vaf1 <- (l2_12 * c3 + l2_31 * c2) / 8
  vaf2 <- (l2_23 * c1 + l2_12 * c3) / 8
  vaf3 <- (l2_31 * c2 + l2_23 * c1) / 8
  
  # obtuse-angle handling (match MATLAB logic)
  ind1 <- (c1 < 0)
  vaf1[ind1] <- fa[ind1] / 2; vaf2[ind1] <- fa[ind1] / 4; vaf3[ind1] <- fa[ind1] / 4
  ind2 <- (c2 < 0)
  vaf1[ind2] <- fa[ind2] / 4; vaf2[ind2] <- fa[ind2] / 2; vaf3[ind2] <- fa[ind2] / 4
  ind3 <- (c3 < 0)
  vaf1[ind3] <- fa[ind3] / 4; vaf2[ind3] <- fa[ind3] / 4; vaf3[ind3] <- fa[ind3] / 2
  
  # accumarray(face(:), [vaf1; vaf2; vaf3])
  va <- numeric(nrow(vertex))
  idx <- c(i1, i2, i3)
  vals <- c(vaf1, vaf2, vaf3)
  s <- tapply(vals, idx, sum)
  va[as.integer(names(s))] <- s
  va
}

# ---- helpers ----------------------------------------------------------------

# per-face area using cross product (vectorized)
# face_area <- function(face, vertex) {
#   i1 <- face[, 1]; i2 <- face[, 2]; i3 <- face[, 3]
#   v1 <- vertex[i1, , drop = FALSE]
#   v2 <- vertex[i2, , drop = FALSE]
#   v3 <- vertex[i3, , drop = FALSE]
#   a  <- v2 - v1
#   b  <- v3 - v1
#   # cross product (row-wise) and its norm
#   cx <- a[, 2] * b[, 3] - a[, 3] * b[, 2]
#   cy <- a[, 3] * b[, 1] - a[, 1] * b[, 3]
#   cz <- a[, 1] * b[, 2] - a[, 2] * b[, 1]
#   0.5 * sqrt(cx * cx + cy * cy + cz * cz)
# }

# build halfedges (source, target) and map to face index
compute_halfedge <- function(face) {
  m  <- nrow(face)
  he <- rbind(
    cbind(face[, 1], face[, 2]),
    cbind(face[, 2], face[, 3]),
    cbind(face[, 3], face[, 1])
  )
  heif <- rep(seq_len(m), each = 3)
  list(he = he, heif = heif)
}

# cot(angle) between row-vectors v1 and v2
vector_cot <- function(v1, v2) {
  cs <- rowSums(v1 * v2)
  d1 <- rowSums(v1 * v1)
  d2 <- rowSums(v2 * v2)
  # guard tiny negatives from round-off in the sqrt
  disc <- pmax(d1 * d2 - cs * cs, 0)
  cs / sqrt(disc)
}
