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


#' Gibbs sampler for Bayesian variables selection with spike-and-slab prior.
#'
#' @param covariate n by p covariate matrix (does not include the column of ones).
#' @param Y continuous response, n by 1.
#' @param nMc number of MCMC iterations.
#' @param v1 p by 1 vector, the large s.d. of beta when gamma_j=1.
#' @param v0 p by 1 vector, the small s.d. of beta when gamma_j=0.
#' @param w the probablity that gamma_j=1 for each variable.
#' @param gamma0 initial value of gamma, contains zeros and ones.
#' @param beta0 initial value of beta.
#' @param sigma0 initial value of sigma2.
#' @param a prior parameters for sigma^2. (Inverse Gamma parameters).
#' @param b prior parameters for sigma^2. (Inverse Gamma parameters).
#'
#' @return a list of objects:
#'    sample.beta All posterior samples of beta.
#'    sample.gamma all posterior samples of gamma.
#'    sample.sigma2 All posterior samples of sigma2.
#'    logodds The log (conditional) posterior odds used when updating gamma.
#'
#' @export
BvsGibbs <- function(covariate, Y, nMc, v1, v0, w, gamma0, sigma0, a, b){
  n <- dim(covariate)[1]
  p <- dim(covariate)[2]
  X <- covariate

  ## Pre-allocation
  save.beta <- matrix(NA, nMc, p)
  save.gamma <- matrix(NA, nMc, p)
  save.odds <- matrix(NA, nMc, p)
  save.sigma2 <- rep(NA, nMc)

  gamma <- gamma0
  sigma2 <- sigma0
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y


  for (i in 1:nMc){
    for (j in 1:p) {
      ## Udpate gamma
      col.index <- setdiff(1:p, j)
      num.M <- X[, col.index]%*%beta[col.index]
      num.V <- 1 / sigma2 * diag(n) - X[, j]%*%t(X[, j]) / (sigma2 * sigma2 * (1 / sigma2 * t(X[, j])%*%X[, j] + 1 / v1[p])[1])
      num.V <- solve(num.V)
      den.M <- num.M
      den.V <- sigma2 * diag(n)
      log.odds <- -1 / 2 * log(det(2 * pi* num.V)) - 1 / 2 * t(Y - num.M) %*% solve(num.V) %*% (Y - num.M) +
        1 / 2 * log(det(2 * pi* den.V)) + 1 / 2 * t(Y - den.M) %*% solve(den.V) %*% (Y - den.M)
      gamma[j] <- ifelse(log.odds > 600, 1, rbinom(1, size = 1, prob = exp(log.odds) / (exp(log.odds) + 1)))

      ## Update beta
      if (gamma[j] != 0 ){
        invK <- solve(t(X[, j])%*%X[, j] / sigma2 + 1 / v1[p])
        M <- invK%*%(t(X[, j])%*%(Y - num.M) / sigma2)
        beta[j] <- rmvnorm(1, M, invK)
      } else {
        beta[j] <- 0
      }
    }
    save.beta[i, ] <- beta
    save.gamma[i,] <- gamma

    ## Update sigma2.
    a.tilde <- a + n / 2
    b.tilde <- sum((Y - X %*% beta) ** 2) / 2 + b
    sigma2 <- rgamma(1, a.tilde, b.tilde)
    save.sigma2[i] <- sigma2
  }

  return(list(sample.beta = save.beta, sample.gamma = save.gamma,
              sample.sigma2 = save.sigma2, logodds = save.odds))
}
