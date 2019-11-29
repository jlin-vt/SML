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


#' Return Euclidean distance matrix for the data.
#'
#' @param points1 data point.
#' @param points2 data point.
#
#' @return distanceMatrix D by D matrix.
#'
#' @seealso \url{http://stackoverflow.com/questions/31571236/my-own-k-means-algorithm-in-r}
#'
#' @export
Euclid <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow = dim(points1)[1], ncol = dim(points2)[1])
  for(i in 1:nrow(points2)) {
    distanceMatrix[, i] <- sqrt(rowSums(t(t(points1) - points2[i, ]) ** 2))
  }

  return(distanceMatrix)
}

#' Demo of k-Means algorithm using the Old Faithful dataset.
#'
#' @param x data.
#' @param centers centers of clusters.
#' @param dist.fun distance function l1/l2.
#' @param nItter number of iterations.
#
#' @return a list of objects:
#'    clusters history of clusters
#'    centers history of centers.
#'
#' @seealso \url{http://www.cs.ubc.ca/~murphyk/Software/VBEMGMM/index.html}
#'
#' @export
Kmeans <- function(x, centers, dist.fun, nItter) {
  cluster.history <- vector(nItter, mode = "list")
  center.history <- vector(nItter, mode = "list")

  for(i in 1:nItter) {
    distsToCenters <- dist.fun(x, centers)
    clusters <- apply(distsToCenters, 1, which.min)
    centers <- apply(x, 2, tapply, clusters, mean)

    # Saving history
    cluster.history[[i]] <- clusters
    center.history[[i]] <- centers
  }
  list(clusters = cluster.history[[nItter]], centers = center.history[[nItter]])
}
