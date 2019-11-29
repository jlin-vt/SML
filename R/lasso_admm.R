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


#' Return proximal function of none smooth function h.
#'
#' @export
Prox.L1 <- function(x, lambda) sign(x) * pmax(abs(x) - lambda, 0)

#' ADMM alogrithm for solving convex optimization with l1-norm.
#'
#' @param rho.
#' @param solATA.
#' @param ATb.
#' @param p.
#' @param maxiter.
#'
#' @return rrs residuals.
#'
#' @seealso \url{https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture19.pdf}
#'
#' @export
ADMM <- function(rho, solATA, ATb, p, maxiter){

  ## Initial
  x <- y <- z <- rep(0, p)
  rrs <- rep(0, maxiter)

  for(i in 1:maxiter){
    ## Store temp y
    y.old <- y

    ## Update x
    x <- solATA %*% (ATb + rho* (y - z))

    ## Update y
    y <- Prox.L1(x + z, lambda / rho)

    ## Update z
    z <- z + x - y

    ## Residual
    rr <- norm(x - y, "2")
    ss <- rho*norm(y - y.old, "2")

    cat(rr, ss, fill = TRUE)

    rrs[i] <- rr
  }

  return(rrs)
}

#' ADMM alogrithm for solving convex optimization with l1-norm.
#'
#' @param rho.
#' @param solATA.
#' @param ATb.
#' @param p.
#' @param maxiter.
#
#' @return rrs residuals.
#'
#' @seealso \url{https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture19.pdf}
#'
#' @export
ADMM.fast <- function(rho, solATA, ATb, p, maxiter){

  ## Initial
  x <- y <- z <- rep(0, p)
  hat.y <- hat.z <- rep(0, p)
  rrs.fast <- rep(0, maxiter)
  alpha <- c(1, rep(0, maxiter - 1))

  for(i in 1:maxiter){
    ## Store temp z
    z.old <- z
    y.old <- y

    ## Update x
    x <- solATA %*% (ATb + rho * (hat.y - hat.z))

    ## Update y
    y <- Prox.L1(x + hat.z, lambda / rho)

    ## Update z
    z <- hat.z + x - y

    ## Residual
    rr <- norm(x - y, "2")
    ss <- rho*norm(y - y.old, "2")

    cat(rr, ss, fill = TRUE)

    rrs.fast[i] <- rr

    # Accelerated stage
    alpha[i + 1] <- (1 + sqrt(1 + 4*alpha[i]))/2
    hat.y <- y + (alpha[i] - 1)/alpha[i + 1] * (y - y.old)
    hat.z <- z + (alpha[i] - 1)/alpha[i + 1] * (z - z.old)
  }

  return(rrs.fast)
}
