test_that("lacf.calc executes with spectrum supplied", {
  skip_on_cran()
  x <- stats::arima.sim(model = list(ar = 0.5), n = 512)
  x.s <- ewspec.trend(x)
  expect_equal(class(lacf.calc(x, spec.est = x.s)), "lacf")
})

test_that("lacf.calc executes with spectrum not supplied", {
  skip_on_cran()
  x <- stats::arima.sim(model = list(ar = 0.5), n = 512)
  expect_equal(class(lacf.calc(x)), "lacf")
})

test_that("lacf.calc warns with large lag.max", {
  skip_on_cran()
  x <- stats::arima.sim(model = list(ar = 0.5), n = 512)
  expect_warning(lacf.calc(x, filter.number = 1, family = "DaubExPhase",
                           lag.max = 2^11),
                 "lag.max too high. Have reset it to  511 . Higher lags are zero")
})

