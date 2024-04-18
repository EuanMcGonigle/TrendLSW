test_that("wav.trend.est executes", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, boundary.handle = FALSE)
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with decimated tranform", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, transform.type = "dec")
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with confidence interval", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, T.CI = TRUE, transform.type = "dec")
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with confidence interval and non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(230) + seq(from = 0, to = 4, length = 230)
  x.t <- suppressWarnings(wav.trend.est(x, T.CI = TRUE, transform.type = "dec"))
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(230) + seq(from = 0, to = 4, length = 230)
  x.t <- suppressWarnings(wav.trend.est(x))
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with confidence interval and boundary handling", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, T.CI = TRUE, transform.type = "dec", boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
})

test_that("T.CI is logical", {
  expect_error(
    wav.trend.est(stats::rnorm(64), T.CI = "true", transform.type = "dec"),
    "Parameter T.CI must be logical variable"
  )
})


test_that("sig.lvl recognised", {
  expect_error(
    wav.trend.est(stats::rnorm(64), sig.lvl = 2),
    "Error: sig.lvl must be a number between 0 and 1."
  )
})

test_that("transform.type recognised", {
  expect_error(
    wav.trend.est(stats::rnorm(64), transform.type = "packet"),
    "Parameter T.transform must be either 'dec' or 'nondec'"
  )
})

test_that("max.scale is integer", {
  expect_error(
    wav.trend.est(stats::rnorm(64), max.scale = 3.7),
    "max.scale parameter must be an integer."
  )
})

test_that("max.scale isn't too large", {
  expect_warning(
    wav.trend.est(stats::rnorm(64), max.scale = 9),
    "max.scale parameter is outside valid range. Resetting to default value."
  )
})

test_that("error for NA data", {
  x <- rep(NA, 64)
  expect_error(
    wav.trend.est(x),
    "Data contains mising values."
  )
})


test_that("data is numeric", {
  x <- c("1", "2")
  expect_error(
    wav.trend.est(x),
    "Data is not numeric"
  )
})
