## Copyright (C) 2017-present, Jiali Lin
## All rights reserved.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.


#' Calculates the covariance matrix sigma using a the squared exponential function
#' with diagnal covraince matrix.
#'
#' @param X1 vectors of data.
#' @param X2 vectors of data.
#'
#' @return Sigma a covariance matrix.
#'
#' @export
CalcSigma <- function(X1, X2, l = 1) {
  Sigma <- matrix(rep(0, length(X1) * length(X2)), nrow = length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i, j] <- exp(-0.5 * (abs(X1[i] - X2[j]) / l) ** 2)
    }
  }
  return(Sigma)
}
