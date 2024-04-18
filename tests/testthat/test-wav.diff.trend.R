test_that("wav.diff.trend.est executes", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s)
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
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

test_that("wav.diff.trend.est executes with soft thresholding", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, thresh.type = "soft")
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with decimated transform", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, transform.type = "dec")
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with decimated transform and boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, transform.type = "dec", boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with soft thresholding and boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, thresh.type = "soft", boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
})


test_that("wav.diff.trend.est executes with hard thresholding", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, thresh.type = "hard")
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with normal = FALSE", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, normal = FALSE)
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with normal = FALSE and boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  x.t <- wav.diff.trend.est(x, x.s, normal = FALSE, boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
})

test_that("wav.diff.trend.est executes with T.CI = TRUE", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  expect_equal(
    class(wav.diff.trend.est(x, x.s, T.CI = TRUE)), "list"
  )
})

test_that("wav.diff.trend.est executes with T.CI = TRUE and confint.type = 'normal'", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.s <- ewspec.diff(x)
  expect_equal(
    class(wav.diff.trend.est(x, x.s, T.CI = TRUE, confint.type = "normal")), "list"
  )
})

test_that("wav.diff.trend.est executes with non-dyadic data and T.CI = TRUE", {
  skip_on_cran()
  x <- stats::rnorm(230) + seq(from = 0, to = 4, length = 230)
  x.s <- suppressWarnings(ewspec.diff(x))
  expect_warning(
    wav.diff.trend.est(x, x.s, T.CI = TRUE),
    "Data length is not power of two. Boundary correction has been applied."
  )
})

test_that("reps recognised as numeric", {
  expect_error(
    wav.diff.trend.est(stats::rnorm(64), reps = "2"),
    "Number of bootstrap replications should be a single positive integer."
  )
})

test_that("reps recognised as integer", {
  expect_error(
    wav.diff.trend.est(stats::rnorm(64), reps = 0.6),
    "Number of bootstrap replications should be a single positive integer."
  )
})
