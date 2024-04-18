test_that("ewspec.diff executes", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.diff(x)
  expect_equal(class(x.s), "list")
})

test_that("ewspec.diff executes with diff.number = 2", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.diff(x, diff.number = 2)
  expect_equal(class(x.s), "list")
})

test_that("ewspec.diff executes with supplied matrix", {
  skip_on_cran()
  x <- stats::rnorm(256)
  A <- Cmat.calc(J = 8)
  A1 <- Atau.mat.calc(8)
  x.s <- ewspec.diff(x, supply.inv.mat = TRUE, inv.mat = solve(2 * A - 2 * A1))
  expect_equal(class(x.s), "list")
})

test_that("Supplied matrix is square", {
  skip_on_cran()
  x <- stats::rnorm(256)
  expect_error(
    ewspec.diff(x, supply.inv.mat = TRUE, inv.mat = matrix(0, 7, 8)),
    "Supplied inverse matrix must be square"
  )
})

test_that("Supplied matrix is large enough", {
  skip_on_cran()
  x <- stats::rnorm(256)
  A <- Cmat.calc(J = 4)
  expect_error(
    ewspec.diff(x, supply.inv.mat = TRUE, inv.mat = solve(A), max.scale = 6),
    "Dimension of supplied inverse matrix must be larger than max.scale"
  )
})


test_that("ewspec.diff warns for boundary handling on non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(200)
  expect_warning(
    ewspec.diff(x),
    "Data length is not power of two. Boundary correction has been applied."
  )
})

test_that("ewspec.diff warns when diff.number not recognised", {
  skip_on_cran()
  x <- stats::rnorm(256)
  expect_warning(
    ewspec.diff(x, diff.number = 3),
    "Function only implements 1st or 2nd differences. Using 1st difference."
  )
})

test_that("ewspec.diff warns when diff.number = 2 and lag != 1", {
  skip_on_cran()
  x <- stats::rnorm(256)
  expect_warning(
    ewspec.diff(x, diff.number = 2, lag = 2),
    "When diff.number = 1, only lag = 1 is supported. Resetting lag to be 1."
  )
})

test_that("lag is numeric", {
  expect_error(
    ewspec.diff(rnorm(64), lag = "1"),
    "The lag parameter should be a positive integer."
  )
})

test_that("lag is positive", {
  expect_error(
    ewspec.diff(rnorm(64), lag = -2),
    "The lag parameter should be a positive integer."
  )
})

test_that("ewspec.diff executes with median smoothing", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.diff(x, smooth.type = "median")
  expect_equal(class(x.s), "list")
})

test_that("ewspec.diff executes with epan smoothing", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.diff(x, smooth.type = "epan")
  expect_equal(class(x.s), "list")
})

test_that("ewspec.diff recognises smoothing", {
  skip_on_cran()
  x <- stats::rnorm(256)
  expect_error(
    ewspec.diff(x, smooth.type = "bartlett"),
    "Smoothing type must be one of 'mean', 'median', or 'epan'."
  )
})
