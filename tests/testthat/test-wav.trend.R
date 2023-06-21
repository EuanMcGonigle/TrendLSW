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

test_that("transform.type recognised", {
  expect_error(
    wav.trend.est(stats::rnorm(64), transform.type = "packet"),
    "Parameter transform.type must be either 'dec' or 'nondec'"
  )
})
