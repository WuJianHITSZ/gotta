//
//  polybool.h
//  GOTTA_v1.6.1
//
//  Created by Jian Wu on 2025/11/12.
//


// polybool.h
// Header for polybool.cpp
// Exposes: polybool() callable from other C++ files
//
// Dependencies: Rcpp, ClipperLib (clipper.hpp)
// Created: 2025-11-12

#ifndef POLYBOOL_H
#define POLYBOOL_H

#include <Rcpp.h>
#include "clipper.hpp"

// The polybool() function performs boolean operations between polygons.
// Parameters:
//   pa, pb  — polygon(s) as R matrices (n×2) or lists of matrices
//   op      — string specifying operation: "or", "and", "notb"/"diff", "xor"
//   ha, hb  — optional LogicalVectors marking holes (same length as list elements)
//   ug_in   — optional NumericVector scaling factor (default 1e6)
// Returns:
//   An Rcpp::List with two elements:
//     $pc — list of resulting polygons (each an n×2 NumericMatrix)
//     $hc — LogicalVector indicating which polygons are holes
//
Rcpp::List polybool(SEXP pa,
                    SEXP pb,
                    std::string op,
                    Rcpp::Nullable<Rcpp::LogicalVector> ha = R_NilValue,
                    Rcpp::Nullable<Rcpp::LogicalVector> hb = R_NilValue,
                    Rcpp::Nullable<Rcpp::NumericVector> ug_in = R_NilValue);

#endif // POLYBOOL_H
