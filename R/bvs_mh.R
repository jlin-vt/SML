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


#' Metropolis–Hastings algorithm for Bayesian variables selection with stochastic search.
#' For sampling scheme, we use switch and swap proposal (Brown, Vannucci, Fearn (1998), JRSS B).
#'
#' @param covariate n by p covariate matrix (does not include the column of ones).
#' @param Y continuous response, n by 1.
#' @param nMc number of MCMC iterations.
#' @param v1 p by 1 vector, the large s.d. of beta when gamma_j=1.
#' @param v0 p by 1 vector, the small s.d. of beta when gamma_j=0.
#' @param w the probablity that gamma_j=1 for each variable.
#' @param a prior parameters for sigma^2. (Inverse Gamma parameters).
#' @param b prior parameters for sigma^2. (Inverse Gamma parameters).
#' @param gamma0:initial value of tao, contains zeros and ones.
#'
#' @return sample.gamma all posterior samples of gamma
#'
#' @export
BvsMH <- function(covariate,Y,nMc,v1,v0, w, gamma0, a, b, phi){
  n <- dim(covariate)[1]
  p <- dim(covariate)[2]
  X <- cbind(1, covariate)

  save.gamma <- matrix(NA, nMc, p)
  save.ratio <- rep(NA, nMc)
  gamma <- gamma0
  phi <- 0.5
  accept.count <- 0

  for (i in 1:nMc){
    # Propose a new gamma
    gamma.n <- gamma
    if (sum(gamma) == p)	 gamma.n[sample(1:p, 1)] <- 0
    if (sum(gamma) == 0) 	 gamma.n[sample(1:p, 1)] <- 1
    if ((sum(gamma) < p) & (sum(gamma) > 0)){
      if (runif(1) > phi) {
        chid <- sample(1:p, 1)
        gamma.n[chid] <- 1 - gamma[chid]
      } else{
        l0 <- which(gamma == 0)
        l1 <- which(gamma == 1)
        gamma.n[l0[sample(1:length(l0), 1)]] <- 1
        gamma.n[l1[sample(1:length(l1), 1)]] <- 0
      }
    }

    Sig.n <- diag(c(10 ** 6, gamma.n * v1 ** 2 + (1 - gamma.n) * v0 ** 2)) # assume that R=I.
    invSig.n <- diag(c(1 / 10 ** 6, gamma.n * (1 / v1 ** 2) + (1 - gamma.n) * (1 / v0 ** 2)))
    A.new <- b + (sum(Y ** 2) - t(Y) %*% X %*% solve(t(X) %*% X + invSig.n) %*% t(X) %*% Y) / 2

    Sig.old <- diag(c(10 ** 6, gamma * v1 ** 2 + (1 - gamma) * v0 ** 2))
    invSig.old <- diag(c(1 / 10 ** 6, gamma * (1 / v1 ** 2) + (1 - gamma)*(1 / v0 ** 2)))
    A.old <- b + (sum(Y ** 2) - t(Y) %*% X %*% solve(t(X) %*% X + invSig.old) %*% t(X) %*% Y) / 2

    log.R <- -log(det(Sig.n %*% (t(X) %*% X) + diag(p + 1))) / 2 - (n / 2 + a) * log(A.new) +
      log(det(Sig.old %*% t(X) %*% X + diag(p + 1))) / 2 + (n / 2 + a) * log(A.old) +
      (sum(gamma.n) - sum(gamma)) * log(w / (1 - w))

    save.ratio[i] <- log.R
    if (log(runif(1)) < log.R) {
      gamma <- gamma.n
      accept.count <- accept.count + 1}

    save.gamma[i,] <- gamma

  }

  return(list(sample.gamma = save.gamma,
              post.ratio = save.ratio,
              accept.count = accept.count)
         )
}
