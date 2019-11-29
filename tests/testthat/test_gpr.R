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

testthat::context("Unit tests for gpr.R")

test_that("input must have same length", {
  set.seed(99)
  n = 10
  x = runif(n)
  output = CalcSigma(x, x)
  expect_equal(nrow(output), n)
  expect_equal(nrow(output), n)
})