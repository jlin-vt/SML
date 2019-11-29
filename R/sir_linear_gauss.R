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
#' @param S.v noise.
#' @param S.w noise.
#' @param phi scale parameter.
#'
#' @return a list of objects:
#'    mu.Xs mean sample.
#'    var.Xs variance sample.
#'
#' @export

SirLinearGauss <- function(y, T, N, S.v, S.w, phi){
  ## Pre-allocation
  N <- 1000
  Xs <- matrix(NA, N, T)
  logweight <- matrix(NA, N, T)
  weight.normalized <- matrix(NA, N, T)
  mu.Xs <- rep(NA, T)
  var.Xs <- rep(NA, T)

  ## Main algorithm
  for (t in 1:T){

    if (t == 1){
      Xs[, t] <- rnorm(N, mean =  0, sd = 1)
    }else{
      ## Propagate particles according to prior
      Xs[, t] <- rnorm(N, mean = phi * Xs[, t-1], sd = S.v)
    }
    logweight[, t] <- dnorm(y[t], mean = Xs[, t], sd = S.w, log = T)
    weight.normalized[, t] <- exp(logweight[, t]) / sum(exp(logweight[, t]))

    ## Resample
    Xs.resample <- sample(Xs[, t], size = N, replace = T, prob = weight.normalized[, t])
    Xs[, t] <- Xs.resample
    mu.Xs[t] <- mean(Xs.resample)
    var.Xs[t] <- var(Xs.resample)
  }

  return(list(mu.Xs = mu.Xs,
              var.Xs = var.Xs))

}
