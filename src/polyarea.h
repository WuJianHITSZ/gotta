//
//  polyarea.h
//  GOTTA_v1.6.1
//
//  Created by Jian Wu on 2025/11/12.
//


// polyarea.h
// Header for polyarea.cpp
// Provides area computation for polygons using the shoelace formula.
//
// Created: 2025-11-12
// Dependencies: Rcpp
// Exposes:
//   - polyarea(List pa): compute areas of a list of polygons
//   - polyarea_matrix(const NumericMatrix& poly): compute area of a single polygon
//
// Usage:
//   #include "polyarea.h"
//   NumericVector A = polyarea(pa);
//   double a = polyarea_matrix(poly);

#ifndef POLYAREA_H
#define POLYAREA_H

#include <Rcpp.h>

// Compute the area of a list of polygons.
// Each element of the list must be a numeric matrix (n × 2).
// Returns a NumericVector of areas (same shape as input list if it had dim attributes).
Rcpp::NumericVector polyarea(Rcpp::List pa);

// Compute the area of a single polygon represented as an n × 2 matrix.
// Returns the absolute area as a double.
double polyarea_matrix(const Rcpp::NumericMatrix& poly);

#endif // POLYAREA_H
