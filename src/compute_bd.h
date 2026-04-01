#ifndef COMPUTE_BD_HPP
#define COMPUTE_BD_HPP

// Compute boundary loop from a triangle mesh (R-style semantics)
// - face: n x 3 integer matrix of 1-based vertex indices (R layout)
// - efficient: if TRUE (or if face indices are non-contiguous), the function
//              remaps face to 1..K, computes the loop, then maps back.
//
// This function replicates the R version:
//
//   adjacency_matrix <- sparseMatrix(i=c(t(face)),
//                                    j=c(t(face[,c(2,3,1)])),
//                                    x=1, dims=max(face) x max(face))
//   boundary_table <- which((adjacency_matrix - t(adjacency_matrix)) > 0,
//                           arr.ind=TRUE)   # NOTE: column-major order
//   boundary <- walk starting at boundary_table[1, 1], always taking the FIRST
//               match where boundary_table[,1] == current, stopping on revisit.
//   boundary <- boundary[seq_len(count) + 1]
//
// The C++ implementation matches that logic exactly, including:
//   - column-major enumeration (by j then i) when forming boundary_table
//   - the “first match” rule during traversal
//   - the final slice [seq_len(count)+1]
//   - auto “efficient” remap when there are gaps in vertex ids
//
// Return: IntegerVector of boundary vertex ids (1-based), same as R.
//
// Requirements:
//   - Link against Rcpp and compile the corresponding .cpp that defines this.
//   - Include this header in translation units that need the declaration.
//
// Example (C++):
//   #include "compute_bd.hpp"
//   using namespace Rcpp;
//   // ... obtain IntegerMatrix face ...
//   IntegerVector bd = compute_bd_cpp(face, false);
//
// Example (R):
//   Rcpp::sourceCpp("compute_bd.cpp")  # Implementation file
//   bd <- compute_bd_cpp(face, efficient = FALSE)

#include <Rcpp.h>

/// Declaration of the exported function implemented in compute_bd.cpp
/// @param face_in  n x 3 IntegerMatrix of 1-based indices (R-style)
/// @param efficient  if true (or auto-triggered), remap to contiguous ids
/// @return IntegerVector boundary loop matching the R implementation
Rcpp::IntegerVector compute_bd(const Rcpp::IntegerMatrix& face_in,
                                   bool efficient = false);

#endif // COMPUTE_BD_HPP
