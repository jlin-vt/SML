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


#' Variational Bayes EM algorithm for Gaussian Mixture Model.
#'
#' @param D the dimension.
#' @param N the number of Data points.
#' @param x training data of size DxN.
#' @param mix gmm model initialize with netlab's GmmVBEM function.
#' @param PriorPar structure containing priors.
#' @param options options for maxIter, threshold etc. etc.
#
#' @return a list of objects:
#'    L likelihood function.
#'    m centers of clusters.
#'
#' @seealso \url{http://www.cs.ubc.ca/~murphyk/Software/VBEMGMM/index.html}
#'
#' @export
GmmVbem <- function(x, mix, PriorPar, options){
  ## Initialize variables
  D <- dim(x)[1];
  N <- dim(x)[2];
  K <- mix$ncentres;
  eps <- 2.2204e-16
  likIncr <- options$threshold + eps;
  L <- rep(likIncr, options$maxIter);
  logLambdaTilde <- matrix(0, 1, K);
  E <- matrix(0, nrow = N, ncol = K);
  trSW <- matrix(0, 1, K);
  xbarWxbar <- matrix(0, 1, K);
  mWm <- matrix(0, 1, K);
  trW0invW <- matrix(0, 1, K);

  ## Prior over the mixing coefficients - CB 10.39
  alpha0 <- PriorPar$alpha;
  ## Prior over gaussian means -  CB 10. 40
  m0 <- PriorPar$mu;
  ## Prior over the gaussian variance - CB 10.40
  beta0 <- PriorPar$beta;
  ## Wishart prior variables - CB 10.40
  W0 <- PriorPar$W;
  W0inv <- solve(W0);
  v0 <- PriorPar$v;

  ## Use 'responsibilities' from initialization to set sufficient statistics -
  ## CB 10.51-10.53.
  Nk <- N * t(mix$priors);
  xbar <- t(mix$centres);
  S <- mix$covars;

  ## Use above sufficient statistics for M step update equations -
  ## CB 10.58, 10.60-10.63
  alpha <- alpha0 + Nk;
  beta <- beta0 + Nk;
  v <- v0 + Nk;
  m <- ((beta0 * m0) %*% matrix(1, 1, K) + (matrix(1, D, 1) %*% t(Nk)) * xbar) / (matrix(1, D, 1) %*% t(beta));
  W <- array(0,c(D, D, K));

  for (k in 1:K) {
    mult1 <- beta0 * Nk[k] / (beta0 + Nk[k]);
    diff3 <- xbar[, k] - m0;
    W[ , , k] <- solve(W0inv + Nk[k] * S[ , ,k]  + mult1 * diff3 %*% t(diff3));
  }

  ## Main loop of algorithm
  for (iter in 1:options$maxIter){
    ## Calculate r - CB 10.64 - 10.66
    psiAlphaHat <- psigamma(sum(alpha), 0);
    logPiTilde <- psigamma(alpha, 0) - psiAlphaHat;
    const <- D * log(2);

    for (k in 1:K){
      t1 <- psigamma(0.5 * kronecker(matrix(1, D, 1), v[k]+1)	 - 0.5 * matrix(1:D), 0);
      logLambdaTilde[k] <- sum(t1) + const + log(det(W[ , , k]));

      for (n in 1:N){
        ## Calculate E
        diff <- x[, n] - m[, k];
        E[n, k] <- D / beta[k] + v[k]  *  t(diff)  %*%  W[ , , k]  %*%  diff;
      }
    }

    ## Calculate rho - CB 10.45 - 10.46
    logRho <- kronecker(matrix(1, N, 1), t(logPiTilde) + 0.5 * logLambdaTilde) - 0.5 * E;
    logSumRho <- apply(logRho, 1, function(x) log(sum(exp(x))));
    logr <- logRho - kronecker(matrix(1, 1, K), logSumRho);
    r <- exp(logr);

    ## Compute N(k) - CB 10.51
    Nk <- exp(apply(logr, 2, function(x) log(sum(exp(x)))))
    ## Add a non-zero term for the components with zero responsibilities
    Nk <- Nk + 1e-10;
    ## Compute xbar(k), S(k) - CB 10.52 - 10.53
    for (k in 1:K){
      xbar[ ,k] <- rowSums(kronecker(matrix(1, D, 1), t(r[, k])) * x) / Nk[k];
      diff1 <- x - kronecker(matrix(1, 1, N), xbar[, k]);
      diff2 <- kronecker(matrix(1, D, 1), t(r[, k])) * diff1;
      S[, , k] <- (diff2 %*% t(diff1)) / Nk[k];
    }

    ## Compute Lower bound (refer to Bishop for these terms) - CB 10.71 - 10.77
    ## C(alpha0)
    logCalpha0 <- lgamma(K * alpha0) - K * lgamma(alpha0);
    ## B(lambda0)
    logB0 <- (v0 / 2) * log(det(W0inv)) - (v0 * D / 2) * log(2) -
      (D * (D - 1) / 4) * log(pi) - sum(lgamma(0.5 * (v0 + 1 -1:D)));
    ## Log(C(alpha))
    logCalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha));
    ## Various other parameters for different terms
    H <- 0;
    for (k in 1:K){
      logBk <- -(v[k] / 2) * log(det(W[ , , k])) - (v[k] * D / 2) * log(2)
      - (D * (D-1) / 4) * log(pi) - sum(lgamma(0.5 * (v[k] + 1 - 1:D)));
      H <- H -logBk - 0.5 * (v[k] -D-1) * logLambdaTilde[k] + 0.5 * v[k] * D;
      ## For Lt1 - third term
      trSW[k] <- sum(diag(v[k] * S[ , , k] %*% W[ , , k]));
      diff <- xbar[, k] - m[, k] ;
      xbarWxbar[k] <- t(diff) %*% W[ , , k] %*% diff;
      ## For Lt4 - Fourth term
      diff <- m[, k] - m0;
      mWm[k] <- t(diff) %*% W[ , , k] %*% diff;
      trW0invW[k] <- sum(diag(W0inv %*% W[ , , k]));
    }

    ## Bishop's Lower Bound
    Lt1 <- 0.5 * sum(Nk * (t(logLambdaTilde) - D / beta - t(trSW) - v * t(xbarWxbar) - D * log(2 * pi)));
    Lt2 <- sum(Nk * logPiTilde)
    Lt3 <- logCalpha0 + (alpha0 -1) * sum(logPiTilde);
    Lt41 <- 0.5 * sum(D * log(beta0 / (2 * pi)) + t(logLambdaTilde) - D * beta0 / beta - beta0 * v * t(mWm));
    Lt42 <- K * logB0 + 0.5 * (v0 - D - 1) * rowSums(logLambdaTilde) - 0.5 * sum(v * t(trW0invW));
    Lt4 <- Lt41 + Lt42;
    Lt5 <- sum(sum(r * logr));
    Lt6 <- sum((alpha - 1) * logPiTilde) + logCalpha;
    Lt7 <- 0.5 * sum(t(logLambdaTilde) + D * log(beta / (2 * pi))) - 0.5 * D * K - H;

    ## Bishop's Lower Bound
    L[iter] <- Lt1 + Lt2 + Lt3 + Lt4 - Lt5 - Lt6 - Lt7; # where is L?

    ## Warning if lower bound decreses
    if (iter > 2 && L[iter] < L[iter - 1]) {
      cat("Lower bound decreased by =", L[iter] - L[iter-1], "\n" )
    }

    ## Begin M step
    ## Compute new parameters - CB 10.58 - 10.63
    alpha <- alpha0 + Nk;
    beta <- beta0 + Nk;
    v <- v0 + Nk;
    m <- (kronecker(matrix(1, 1, K), beta0 * m0) + kronecker(matrix(1, D, 1),  t(Nk)) * xbar) / kronecker(matrix(1, D, 1), t(beta));
    for (k in 1:K) {
      mult1 <- beta0 * Nk[k] / (beta0 + Nk[k]);
      diff3 <- xbar[, k] - m0;
      W[, , k]  <- solve(W0inv + Nk[k] * S[, , k]  + mult1 * diff3 %*% t(diff3));
    }

    ## Interactive plots
    if (options$displayIter){
      cat("iter =", iter, "\n" )
    }

    if (options$displayFig){
      dev.off()
      plot(x[1, ], x[2, ], xlab = "eruptions", ylab = "waiting")
      for (i in 1:K)  ellipse(m[, i], solve(W[, , i]) / (v[i] - D - 1), col = i)
      Sys.sleep(0.1)
    }

    if (iter>1){
      likIncr <- abs((L[iter] - L[iter - 1]) / L[iter - 1]);
    }

    if (likIncr < options$threshold) break
  }

  return(list(L = L, center = m))
}
