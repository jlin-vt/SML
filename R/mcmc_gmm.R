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


#' Metropolis Hasting samlinging MCMC for 1-d Gaussian Mixture Model.
#'
#' @param p target distribution.
#' @param N number of iterations need to run.
#'
#' @return xs samples.
#'
#' @export
MetropolisHastingSampler <- function(p, N){
  # Initial
  xs <- rep(0, N)
  xs[1] <- 0.3

  for (i in 1:N){
    # Proposal distribution
    xs_ast <- rnorm(1, 2, var[j])
    alpha <- p(xs_ast) / p(xs[i])

    # Accept the sample with prob = min(alpha,1)
    if (runif(1) < min(alpha, 1)){
      xs[i + 1] <- xs_ast
    } else {
      xs[i + 1] <- xs[i]
    }
  }

  return(xs)
}
