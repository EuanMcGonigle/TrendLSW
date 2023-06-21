test_that("wav.diff.trend.est executes", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s)
  expect_equal(class(x.t), "numeric")
})

test_that("wav.diff.trend.est executes with boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, boundary.handle = TRUE)
  expect_equal(class(x.t), "numeric")
})

test_that("wav.diff.trend.est executes with non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(230) + seq(from = 0, to = 4, length = 230)
  x.s <- suppressWarnings(ewspec.diff(x))
  expect_warning(
    wav.diff.trend.est(x, x.s),
    "Data length is not power of two. Boundary correction has been applied."
  )
})

test_that("wav.diff.trend.est executes with hard thresholding", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, thresh.type = "hard")
  expect_equal(class(x.t), "numeric")
})

test_that("wav.diff.trend.est executes with normal = FALSE", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, normal = FALSE)
  expect_equal(class(x.t), "numeric")
})

test_that("wav.diff.trend.est executes with normal = FALSE and boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, normal = FALSE, boundary.handle = TRUE)
  expect_equal(class(x.t), "numeric")
})
