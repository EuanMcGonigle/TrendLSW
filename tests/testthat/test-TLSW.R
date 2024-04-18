test_that("TLSW executes", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW(x)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW executes without spec est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW(x, do.spec.est = FALSE)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW executes with differenced spec est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW(x, S.do.diff = TRUE, do.trend.est = FALSE)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW executes with nonlinear trend est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW(x, T.est.type = "nonlinear")
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW recognises T.thresh.type", {
  skip_on_cran()
  x <- stats::rnorm(128)
  expect_error(
    TLSW(x, T.thresh.type = "bayes", T.est.type = "nonlinear"),
    "The parameter T.thresh.type must be either 'hard' or 'soft'."
  )
})

test_that("TLSW executes with supplied inv mat", {
  skip_on_cran()
  x <- stats::rnorm(128)
  inv.mat <- solve(Cmat.calc(J = log2(128)))
  x.lsw <- TLSW(x, S.inv.mat = inv.mat)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW executes with nonlinear wavelet trend estimator using non-differenced spec est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW(x, T.CI = TRUE, T.est.type = "nonlinear")
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW gives same output: spec", {
  skip_on_cran()
  spec <- matrix(0, nrow = 9, ncol = 512)
  spec[1, ] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
  trend <- seq(from = 0, to = 5, length = 512)
  set.seed(1)
  x <- TLSWsim(trend = trend, spec = spec)
  x.TLSW <- TLSW(x)
  expect_snapshot_value(x.TLSW$spec.est$S$D, style = "deparse")
})

test_that("TLSW gives same output: trend", {
  skip_on_cran()
  spec <- matrix(0, nrow = 9, ncol = 512)
  spec[1, ] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
  trend <- seq(from = 0, to = 5, length = 512)
  set.seed(1)
  x <- TLSWsim(trend = trend, spec = spec)
  x.TLSW <- TLSW(x)
  expect_snapshot_value(x.TLSW$trend.est, style = "deparse")
})
