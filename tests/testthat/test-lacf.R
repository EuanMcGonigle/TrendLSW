test_that("TLSWlacf executes", {
  skip_on_cran()
  x <- stats::arima.sim(model = list(ar = 0.5), n = 512)
  x.TLSW <- TLSW(x)
  expect_equal(class(TLSWlacf(x.TLSW)), "lacf")
})

test_that("TLSWlacf warns with large lag.max", {
  skip_on_cran()
  x <- stats::arima.sim(model = list(ar = 0.5), n = 512)
  x.TLSW <- TLSW(x, S.filter.number = 1)
  expect_warning(
    TLSWlacf(x.TLSW,
      lag.max = 2^11
    ),
    "lag.max too high. Have reset it to  511 . Higher lags are zero"
  )
})

test_that("TLSWlacf rejects negative lag.max", {
  skip_on_cran()
  x <- stats::rnorm(64)
  x.TLSW <- TLSW(x)
  expect_error(
    TLSWlacf(x.TLSW, lag.max = -4),
    "Argument lag.max should be a nonegative integer."
  )
})


test_that("TLSWlacf rejects none TLSW object", {
  skip_on_cran()
  x <- stats::rnorm(64)
  expect_error(
    TLSWlacf(x),
    "Argument x.TLSW should be an object of class TLSW."
  )
})
