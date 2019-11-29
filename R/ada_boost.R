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


#' Demo of Adaboost using Sonar dataaset.
#'
#' The base learner is based on Recursive Partitioning and Regression Trees.
#'
#' @param train training dataset.
#' @param test test dataset.
#' @param index.x index of X.
#' @param index.y index of Y.
#' @param n.iter the number of iterations.
#'
#' @return a list of objects:
#'    train.error errors from training dataset.
#'    test.error errors from test dataset.
#'
#' @seealso \url{http://web.stanford.edu/class/stats202/}
#'
#' @export
AdaBoost <-  function(train, test, index.x, index.y, n.iter){
  ## Initial
  n.train <- nrow(train)
  n.test <- nrow(test)
  y <- train[, index.y]
  x <- train[, index.x]
  y.test <- test[, index.y]
  x.test <- test[, index.x]
  train.error <- rep(0, n.iter)
  test.error <- rep(0, n.iter)
  y.hat <- rep(0, n.train)
  y.hat.test <- rep(0, n.test)
  i <- 1

  ## Main boosting algorithm
  while(i <= n.iter){
    w <- exp(-y * y.hat)
    w <- w / sum(w)
    fit <- rpart(y~., x, w, method = "class")
    g <- -1 + 2 * (predict(fit, x)[, 2] > .5)
    g.test <- -1 + 2 * (predict(fit, x.test)[, 2] > .5)
    e <- sum(w * (y * g < 0))
    alpha <- .5 * log((1 - e) / e)
    alpha <- 0.1 * alpha
    y.hat <- y.hat + alpha * g
    y.hat.test <- y.hat.test + alpha * g.test
    train.error[i] <- sum(1 * y.hat * y < 0) / n.train
    test.error[i] <- sum(1 * y.hat.test * y.test < 0) / n.test
    i <- i + 1
  }

  return(list(train.error = train.error,
              test.error = test.error)
         )
}
