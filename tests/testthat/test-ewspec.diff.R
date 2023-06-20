test_that("ewspec.diff executes", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.diff(x)
  expect_equal(class(x.s), "list")
})

test_that("ewspec.trend executes with supplied matrix", {
  skip_on_cran()
  x <- stats::rnorm(256)
  A <- Cmat.calc(J=8)
  x.s <- ewspec.trend(x, supply.mat = TRUE, mat = solve(A))
  expect_equal(class(x.s), "list")
})

test_that("ewspec.diff warns for boundary handling on non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(200)
  expect_warning(ewspec.diff(x),
                 "Data length is not power of two. Boundary correction has been applied.")
})

test_that("lag is numeric", {
  expect_error(ewspec.diff(rnorm(64), lag = "1"),
               "The lag parameter should be a positive integer.")
})

test_that("lag is positive", {
  expect_error(ewspec.diff(rnorm(64), lag = -2),
               "The lag parameter should be a positive integer.")
})

