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


#' SMC method using SIR algorithm in linear Gaussian state-space model.
#'
#' @param y obeservations.
#' @param T time.
#' @param N size.
#' @param sigma noise.
#' @param phi scale parameter.
#'
#' @return a list of objects:
#'    mu.Xs mean sample.
#'    var.Xs variance sample.
#'
#' @export

SirStochasticVolatility <- function(y, T, N, gamma, sigma, phi){
  ## Pre-allocation
  N <- 1000
  Xs <- matrix(NA, N, Time)
  logweight <- matrix(NA, N, Time)
  weight.normalized <- matrix(NA, N, Time)
  mu.Xs <- rep(NA, Time)
  var.Xs <- rep(NA, Time)

  ## Main algorithm
  for (t in 1:Time){
    if (t == 1){
      Xs[, t] <- rnorm(N, mean = 0, sd = sigma / sqrt(1 - phi^2))
    } else {
      ## Propagate particles according to prior
      Xs[, t] <- rnorm(N, mean = phi * Xs[, t-1], sd = sigma)
    }
    logweight[, t] <- dnorm(y[t], mean = 0, sd = exp(gamma + Xs[, t]), log = T)

    weight.normalized[,t] <- exp(logweight[, t])/sum(exp(logweight[, t]))

    ## Resample
    Xs.resample <- sample(Xs[, t], size = N, replace = T, prob = weight.normalized[, t])
    Xs[, t] <- Xs.resample
    mu.Xs[t] <- mean(Xs.resample)
    var.Xs[t] <- var(Xs.resample)
  }

  return(list(mu.Xs = mu.Xs,
              var.Xs = var.Xs))

}
