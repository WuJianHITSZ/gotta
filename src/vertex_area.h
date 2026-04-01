//
//  vertex_area.h
//  GOTTA
//
//  Created by Jian Wu on 2025/11/4.
//


#ifndef VERTEX_AREA_H
#define VERTEX_AREA_H

#include <Rcpp.h>
#include <string>

// ---- exported function declarations ----
// (useful if other .cpp files want to call them directly)
Rcpp::NumericVector face_area(const Rcpp::IntegerMatrix& face,
                              const Rcpp::NumericMatrix& vertex);

Rcpp::NumericVector vertex_area(const Rcpp::IntegerMatrix& face,
                                const Rcpp::NumericMatrix& vertex,
                                std::string type = "one_ring");

#endif  // VERTEX_AREA_H
