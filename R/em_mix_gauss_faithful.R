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


#' Demo of EM algorithm for univariate Gaussian Mixture Model.
#' using the Old Faithful dataset.
#'
#' @param params: initial parameters (means, covariances, cluster probability).
#' @param X: data.
#' @param clusters: number of clusters desired.
#' @param tol: tolerance.
#' @param maxits: maximum iterations.
#' @param showits: whether to show iterations.
#
#' @return a list of objects:
#'    probs: probability.
#'    mu: mean.
#'    var: variance.
#'    resp: responsibilities.
#'    cluster: cluster assignment.
#'
#' @seealso \url{https://github.com/m-clark/Miscellaneous-R-Code/tree/master/ModelFitting/EM%20Examples}
#'
#' @export
GaussmixEM <- function(params, X, clusters = 2, tol = .00001, maxits = 100, showits = T){
  ## Starting points
  N <- nrow(X)
  nams <- names(params)
  mu <- params$mu
  var <- params$var
  probs <- params$probs

  ## Initials
  ri <- matrix(0, ncol = clusters, nrow = N)
  it <- 0
  converged <- FALSE

  ## Show iterations
  if (showits)
    cat(paste("Iterations of EM:", "\n"))

  while ((!converged) & (it < maxits)) {
    probsOld <- probs
    muOld <- mu
    varOld <- var
    riOld <- ri

    ## E step: compute responsibilities
    for (k in 1:clusters){
      ri[, k] <- probs[k] * dnorm(X, mu[k], sd = sqrt(var[k]), log = F)
    }
    ri <- ri / rowSums(ri)

    ## M step: compute the weighted average cluster membership size
    rk <- colSums(ri)
    probs <- rk / N
    mu <- (t(X) %*% ri) / rk
    var <- (t(X ** 2) %*% ri) / rk - mu ** 2

    parmlistold <- rbind(probsOld, muOld, varOld)
    parmlistcurrent <- rbind(probs, mu, var)
    it <- it + 1

    ## If showits true, & it =1 or divisible by 5 print message
    if (showits & it == 1 | it%%5 == 0){
      cat(it, "\n")
    }
    converged <- max(abs(parmlistold - parmlistcurrent)) < tol
  }

  ## Create cluster membership
  clust <- which(round(ri) == 1, arr.ind = T)

  ## Order accoring to row rather than cluster
  clust <- clust[order(clust[, 1]), 2]
  out <- list(probs = probs, mu = mu, var = var, resp = ri, cluster = clust)
}
