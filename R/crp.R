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


#' Demo of chinese restaurant process.
#'
#' @param num.customers the number of customers (observaions).
#' @param alpha concentration parameter.
#'
#' @return Plot of clusetring results.
#'
#' @seealso \url{https://github.com/tbroderick/bnp_tutorial/tree/2016mlss}
#'
#' @export
CRP <- function(num.customers, alpha) {
  table <- 1
  next.table <- 2

  ## Main algorithm
  for (i in 1:num.customers) {
    if (runif(1, 0, 1) < alpha/(alpha + i)) {
      ## Add a new ball color
      table <- c(table, next.table)
      next.table <- next.table + 1
    } else {
      ## Pick out a ball from the urn, and add back a ball of the same color
      select.table <- table[sample(1:length(table), 1)]
      table <- c(table, select.table)
    }

    ## Interactive plot cluster assignments in order of appearance
    K <- unique(table)
    plot(seq(1, i + 1), table, xlab = "Sample index", ylab = "Cluster by order of appearance",
         xlim = c(0, max(10, i)), ylim = c(0, max(10, length(table))),
         pch = 19, main = bquote(rho ~ "~Dirichlet"  # ~'('~.(a)~',...,'~.(a)~')'
                                 ~ ", K = " ~ .(K)))
    ## Press 'x' to stop
    line <- readline()
    if (line == "x")
      return("done")
  }

}
