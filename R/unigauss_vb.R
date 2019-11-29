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


#' Univariate Variational Bayesian for 1-d Gaussian.
#'
#' @param data data.
#' @param prior a list of priors.
#' @param post a list of posteriors.
#' @param Options a list of options.
#' @param TrueNormalGammaPdf true joint distribution.
#' @param vbPost estimated joint normal-gmma density.
#'
#' @return a list of objects:
#'    Lbound Lower bound of likelihood function.
#'    aN hyperparameter of lambda.
#'    bN hyperparameter of lambda.
#'    muN hyperparameter of mu.
#'    kappaN hyperparameter of mu.

UniGaussVb = function (data, prior, post, Options, TrueNormalGammaPdf, vbPost){
  ## Estimation
  m <- mean(data)
  s <- sum(data)
  sSq <- sum(data*data)
  xbar <- m
  sigma2Hat <- sSq/N - xbar ** 2

  ## Load initial prior
  a0 <- prior$a0
  b0 <- prior$b0
  mu0 <- prior$mu0
  kappa0 <- prior$kappa0

  ## Load posterior
  aN <<- post$aN
  bN <<- post$bN
  muN <<- post$muN
  kappaN <<- post$kappaN

  ## Options for the algorithm
  maxIter <- Options$maxIter
  iter <- Options$iter
  Lq <- Options$Lq
  converged <- Options$converged
  tol <- Options$tol
  Lbound <- Options$Lbound

  ## Main algorithm
  while ((!converged) & (iter < maxIter)) {
    LqOld <- Lq

    ## Update q(mu)
    elambda <- aN / bN
    muN <- (kappa0 * mu0 + N * m) / (kappa0 + N)
    kappaN <- (kappa0 + N) * elambda

    cat("iter ", iter, "\n")
    cat("aN ", aN, "\n")

    if (iter == 1){
      dens <- outer(mu, lambda, TrueNormalGammaPdf)
      contour(mu, lambda, dens,
              xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)
      dens <- outer(mu, lambda, vbPost)
      contour(mu, lambda, dens,
              xlab = 'mu', ylab = 'lambda', col = 'red', add = T, drawlabels = F)
    }

    ## Update q(lambda)
    emu <- muN
    emuSquare <- 1 / kappaN + muN ** 2
    aN <- a0 + (N + 1) / 2
    bN <- b0 + 1/2 * ((sSq + kappa0 * mu0 ** 2) - 2 * emu *(s + kappa0 * mu0) + emuSquare * (kappa0 + N))
    if (iter == 1){
      dens <- outer(mu, lambda, TrueNormalGammaPdf)
      contour(mu, lambda, dens,
              xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)
      dens <- outer(mu, lambda, vbPost)
      contour(mu, lambda, dens,
              xlab = 'mu', ylab = 'lambda', col = 'red', add = T, drawlabels = F)
    }

    ## Lower bound
    Lq <- 1/2 * log(1 / kappaN) + log(gamma(aN)) - aN * log(bN)
    Lbound[iter] <- Lq
    if (iter > 1){
      if (Lq - LqOld < -tol){
        cat('Lq did not increase, iter =', iter, "\n")
      } else if (abs(Lq - LqOld) < tol){
        converged = TRUE
      }
    }

    iter <- iter + 1

    ## Plot approximates
    cat('Total # of iterations:', iter, "\n")
    dens <- outer(mu, lambda, TrueNormalGammaPdf)
    contour(mu, lambda, dens,
            xlab = 'mu', ylab = 'lambda', col = 'green', drawlabels = F)
    dens <- outer(mu, lambda, vbPost)
    contour(mu, lambda, dens,
            xlab = 'mu', ylab = 'lambda', col = 'red', add = T, drawlabels = F)

    return(list(Lbound = Lbound, aN = aN, bN = bN, kappaN = kappaN, muN = muN))
  }
}
