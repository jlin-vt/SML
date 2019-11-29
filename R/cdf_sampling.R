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


#' Find CDF
#'
#' @param f a function of x.
#'
#' @return F cdf.
#'
#' @export
CDF = function(f){
  F <- function(x) {integrate(f, 0, x)$value}
  F <- Vectorize(F)
  return(F)
}

#' Find inverse CDF
#'
#' @param F a cdf of f(x).
#' @param l lower bound of domain.
#' @param u upper bound of domain.
#'
#' @return F.inv inverse cdf.
InverseCDF = function(F, l, u){
  F.inv <- function(y){uniroot(function(x){F(x) - y}, interval = c(l, u))$root}
  F.inv <- Vectorize(F.inv)
  return(F.inv)
}
