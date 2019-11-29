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

#' Demo of polya urn process.
#
#   Assume F ~ DP(alpha, G0) with alpha = 5, G0 ~ N (0, 1).
#   Let theta.i ~ F, generate 10000 samples from F using the Polya urn scheme.
#'
#' @param G_0 Base measure.
#' @param N the number of customers (observaions).
#' @param alpha concentration parameter.
#'
#' @return balls clusetring results.
#'
#' @export
PolyaUrn <- function(G_0, N, alpha) {
  balls <- c()

  ## Main algorithm
  for (i in 1:N) {
    if (runif(1) < alpha / (alpha + length(balls))) {
      ## Add a new ball color.
      new_color <- G_0()
      balls <- c(balls, new_color)
    }
    else {
      ## Pick out a ball from the urn, and add back a ball of the same color.
      ball <- balls[sample(1:length(balls), 1)]
      balls <- c(balls, ball)
    }
  }

  return(balls)
}
