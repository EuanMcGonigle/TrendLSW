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

test_that("wav.trend.est executes with nondecimated tranform", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, transform.type = "nondec")
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with confidence interval", {
  skip_on_cran()
  x <- stats::rnorm(256) + seq(from = 0, to = 4, length = 256)
  x.t <- wav.trend.est(x, calc.confint = TRUE)
  expect_equal(class(x.t), "list")
})

test_that("wav.trend.est executes with confidence interval and non-dyadic data", {
  skip_on_cran()
  x <- stats::rnorm(230) + seq(from = 0, to = 4, length = 230)
  x.t <- suppressWarnings(wav.trend.est(x, calc.confint = TRUE))
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
  x.t <- wav.trend.est(x, calc.confint = TRUE, boundary.handle = TRUE)
  expect_equal(class(x.t), "list")
})

test_that("calc.confint is logical", {
  expect_error(
    wav.trend.est(stats::rnorm(64), calc.confint = "true"),
    "Parameter calc.confint must be logical variable"
  )
})

test_that("calc.confint can't be used with nondec transform", {
  expect_error(
    wav.trend.est(stats::rnorm(64), calc.confint = TRUE, transform.type = "nondec"),
    "Only transform.type = 'dec' is supported for calculating confidence
            intervals"
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
    "Parameter transform.type must be either 'dec' or 'nondec'"
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
