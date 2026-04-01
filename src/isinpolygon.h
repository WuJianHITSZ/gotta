//
//  isinpolygon.h
//  GOTTA_v1.6.1
//
//  Created by Jian Wu on 2025/11/12.
//


// isinpolygon.h
// Header for isinpolygon.cpp
// Provides point-in-polygon testing using the ray-casting algorithm.
//
// Created: 2025-11-12
// Dependencies: Rcpp
// Exposes:
//   - isinpolygon(const NumericMatrix& polygon, const NumericMatrix& xy)
//
// Usage:
//   #include "isinpolygon.h"
//   LogicalVector inside = isinpolygon(polygon, xy);

#ifndef ISINPOLYGON_H
#define ISINPOLYGON_H

#include <Rcpp.h>

// Test whether points (xy) lie inside a polygon.
// Parameters:
//   polygon — NumericMatrix (m × 2) representing polygon vertices
//   xy      — NumericMatrix (n × 2) representing query points
// Returns:
//   LogicalVector (length n): TRUE if inside, FALSE if outside, NA if point has NA.
Rcpp::LogicalVector isinpolygon(const Rcpp::NumericMatrix& polygon,
                                const Rcpp::NumericMatrix& xy);

#endif // ISINPOLYGON_H
