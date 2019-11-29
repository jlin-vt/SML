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


#' Gibbs sampler for Bayesian variables selection requires with continuous gaussian.
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

  save.beta <- matrix(NA, nMc, p)
  save.gamma <- matrix(NA, nMc, p)
  save.odds <- matrix(NA, nMc, p)
  save.sigma2 <- rep(NA, nMc)

  gamma <- gamma0
  sigma2 <- sigma0

  for (i in 1:nMc){
    ## Update beta
    invSIG <- diag(c(gamma * (1 / v1 ** 2) + (1 - gamma) * (1 / v0 ** 2)))
    K <- t(X) %*% X / sigma2 + invSIG
    invK <- solve(K)
    M <- (t(X) %*% Y) / sigma2
    beta <- t(rmvnorm(n = 1, mean = invK%*%M, sigma = invK))
    save.beta[i, ] <- beta

    ## Update gamma from (conditional) posterior odds
    gamma <- rep(NA, p)
    log.odds <- -log(v1) - beta[1:p] ** 2 / (2 * v1 ** 2) + log(v0) +
      beta[1:p] ** 2 / (2 * v0 ** 2) + log(w / (1 - w))

    save.odds[i, ] <- log.odds
    for (j in 1:p){
      if (log.odds[j] > 0)  gamma[j] <- 1
      else {
        gamma[j] <- rbinom(1, size = 1, prob = exp(log.odds[j]) / (exp(log.odds[j]) + 1))
      }
    }
    save.gamma[i,] <- gamma

    ## Update sigma2
    a.tilde <- a + n / 2
    b.tilde <- sum((Y - X %*% beta) ** 2) / 2 + b
    sigma2 <- rgamma(1, a.tilde, b.tilde)
    save.sigma2[i] <- sigma2
  }
  return(list(sample.beta = save.beta, sample.gamma = save.gamma,
              sample.sigma2 = save.sigma2, logodds = save.odds))
}
