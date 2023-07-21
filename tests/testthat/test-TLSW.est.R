test_that("TLSW.est executes", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW.est(x)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW.est executes without spec est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW.est(x, do.spec.est = FALSE)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW.est executes with differenced spec est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW.est(x, WP.do.diff = TRUE, do.trend.est = FALSE)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW.est executes with nonlinear trend est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW.est(x, T.est.type = "nonlinear")
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW.est executes with supplied inv mat", {
  skip_on_cran()
  x <- stats::rnorm(128)
  inv.mat <- solve(Cmat.calc(J = log2(128)))
  x.lsw <- TLSW.est(x, WP.inv.mat = inv.mat)
  expect_equal(class(x.lsw), "TLSW")
})

test_that("TLSW.est executes with nonlinear wavelet trend estimator using non-differenced spec est", {
  skip_on_cran()
  x <- stats::rnorm(128)
  x.lsw <- TLSW.est(x, T.confint = TRUE, T.est.type = "nonlinear")
  expect_equal(class(x.lsw), "TLSW")
})
