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

#' Rejection sampler for univariate data.
#'
#' @param theta.list linspace.
#' @param f target distribution.
#' @param g proposed distribution.
#' @param M constant determinig the accptance rate.
#'
#' @return sample samples
#'
#' @export
RejectionSampler <- function(theta.list, f, g, M){

  ## Define target and proposal
  f.theta <- f(theta.list)
  g.theta <- g(theta.list)

  ## Pre-allocation
  sample <- rep(NA, 1000)
  tag <- 0

  ## Main algorithm
  while (tag <= 1000){
    theta.i <- rexp(1, rate = 1)
    f.thetai <- sqrt(2 / pi) * exp(-theta.i ** 2 / 2)
    M.gthetai <- M * dexp(theta.i, rate = 1)
    prob.i <- f.thetai/M.gthetai
    u <- runif(1)
    if (u < prob.i){
      tag <- tag+1
      sample[tag] <- theta.i
    }
  }

  return(sample)
}
