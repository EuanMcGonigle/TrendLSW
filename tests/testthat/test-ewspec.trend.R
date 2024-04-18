test_that("ewspec.trend executes", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.trend(x)
  expect_equal(class(x.s), "list")
})

test_that("ewspec.trend executes using different an and gen wavelets", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.trend(x, an.filter.number = 4, gen.filter.number = 5)
  expect_equal(class(x.s), "list")
})

test_that("ewspec.trend executes with supplied matrix", {
  skip_on_cran()
  x <- stats::rnorm(256)
  A <- Cmat.calc(J = 8)
  x.s <- ewspec.trend(x, supply.inv.mat = TRUE, inv.mat = solve(A))
  expect_equal(class(x.s), "list")
})

test_that("Supplied matrix is square", {
  skip_on_cran()
  x <- stats::rnorm(256)
  expect_error(
    ewspec.trend(x, supply.inv.mat = TRUE, inv.mat = matrix(0, 7, 8)),
    "Supplied inverse matrix must be square"
  )
})

test_that("Supplied matrix is large enough", {
  skip_on_cran()
  x <- stats::rnorm(256)
  A <- Cmat.calc(J = 4)
  expect_error(
    ewspec.trend(x, supply.inv.mat = TRUE, inv.mat = solve(A), max.scale = 6),
    "Dimension of supplied inverse matrix must be larger than max.scale"
  )
})

test_that("ewspec.trend warns for boundary handling on non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(200)
  expect_warning(
    ewspec.trend(x),
    "Data length is not power of two. Boundary correction has been applied."
  )
})


test_that("error for NA data", {
  x <- rep(NA, 64)
  expect_error(
    ewspec.trend(x),
    "Data contains mising values."
  )
})

test_that("data is numeric", {
  x <- c("1", "2")
  expect_error(
    ewspec.trend(x),
    "Data is not numeric"
  )
})

test_that("binwidth is integer", {
  expect_error(
    ewspec.trend(rnorm(64), binwidth = 9.7),
    "binwidth parameter must be an integer."
  )
})

test_that("max.scale is integer", {
  expect_error(
    ewspec.trend(rnorm(64), max.scale = 3.7),
    "max.scale parameter must be an integer."
  )
})

test_that("max.scale isn't too large", {
  expect_warning(
    ewspec.trend(rnorm(64), max.scale = 9),
    "max.scale parameter is outside valid range. Resetting to default value."
  )
})

test_that("ewspec.trend executes with median smoothing", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.trend(x, smooth.type = "median")
  expect_equal(class(x.s), "list")
})

test_that("ewspec.trend executes with epan smoothing", {
  skip_on_cran()
  x <- stats::rnorm(256)
  x.s <- ewspec.trend(x, smooth.type = "epan")
  expect_equal(class(x.s), "list")
})

test_that("ewspec.trend recognises smoothing", {
  skip_on_cran()
  x <- stats::rnorm(256)
  expect_error(
    ewspec.trend(x, smooth.type = "bartlett"),
    "Smoothing type must be one of 'mean', 'median', or 'epan'."
  )
})
