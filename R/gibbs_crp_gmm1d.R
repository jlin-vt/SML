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


#' Gibbs sampling for the Gaussian Dirichlet process mixture model.
#'
#' @param x data (vector of length N).
#' @param alpha DP concentration parameter.
#' @param tau0 normal-gamma prior shape.
#' @param beta0 normal-gamma prior rate (inverse scale).
#' @param kappa0 normal-gamma prior precision scaling parameter.
#' @param nIter number of Gibbs iterations.
#
#' @return C N x nIter matrix of cluster assignments.
#'
#' @export
DpmmGibbs <- function(x, alpha, tau0, beta0, mu0, kappa0, nIter) {
  N <- length(x)
  p <- rep(1, 1, N)
  C <- matrix(data = NA, nrow = N, ncol = nIter); c <-rep(1, 1, N)
  m <- rep(0, 1, N); m[1] <- N
  logpost <- rep(1, 1, nIter); logprior <- rep(1, 1, N); loglik <- rep(1, 1, N)
  ix <- 1:N

  for (i in 1:nIter)
  {
    print(paste("Iteration", i))

    for (n in 1:N)
    {
      ## All customers except n
      cn <- c[1:N!=n]

      ## Count cluster assignments
      for (k in 1:N)
      {
        m[k] <- sum(cn == k)
      }

      if (all(m > 0))
      {
        ## Active dishes
        K.active <- ix[m > 0]
      }
      else
      {
        ## Active dishes + 1 new dish
        K.active <- c(ix[m > 0], min(ix[m == 0]))
      }

      for (k in K.active)
      {
        if (m[k] > 0)
        {
          ## Prior for old dish
          logprior[k] <- log(m[k])
        }
        else
        {
          ## Prior for new dish
          logprior[k] <- log(alpha)
        }
        ## Calculate log likelihood
        loglik[k] <- DpmmLoglik(x[n], x[c==k], tau0, beta0, mu0, kappa0)
      }

      ## Posterior
      post <- NormalizeLogprob(logprior[K.active] + loglik[K.active])

      ## Update cluster assignment
      c[n] <- sample(K.active, 1, rep = TRUE, prob = post)
    }

    C[,i] <- c
  }

  return (C)
}

#' Generalized t-Distribution log PDF
T_logpdf <- function(x, mu, v, df){
  g <- lgamma(0.5 * df + 0.5) - lgamma(0.5 * df) - log(sqrt(df * pi * v))
  return(g - 0.5 * (df + 1) * log(1 + (1 / df) * ((x - mu) ** 2) / v))
}

#' Normalize a log probability vector without numerical underflow/overflow
#'
#' @export
NormalizeLogprob <- function(log.prob){
  ## Compute the log of the normalizing constant
  g <- log(sum(exp(log.prob - max(log.prob)))) + max(log.prob)

  ## Find probabilities by exponentiating normalized log probabilities
  return(exp(log.prob - g))
}

#' Calculate the log likelihood for the Gaussian DP mixture model
#' Mean and variance parameters marginalized under normal-gamma prior
#' This corresponds to a generalized Student's t-distribution
#' NB: some constant terms are ignored
#'
#' @export
DpmmLoglik <- function(xn, x, tau0, beta0, mu0, kappa0){
  N <- length(x)
  kappaN <- kappa0 + N
  tauN <- tau0 + N / 2

  if (N > 0){
    ## If there previous observations, use the current posterior
    xm <- mean(x)
    betaN <- beta0 + 0.5 * sum((x - xm) ** 2) + (kappa0 * N *(xm - mu0) ** 2) / (2 * kappaN)
    muN <- (kappa0 * mu0 + N * xm) / kappaN
  } else {
    ## If there are no previous observations, revert to the prior
    betaN <- beta0
    muN <- mu0
  }

  return(T_logpdf(xn, muN, betaN * (kappaN + 1) / (tauN * kappaN), 2 * kappaN))
}
