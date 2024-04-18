spec1 <- wavethresh::cns(256)
spec1 <- wavethresh::putD(spec1, level = 6, rep(1, 256))

spec2 <- matrix(0, 8, 2^8)
spec2[1, ] <- rep(1, 256)

spec3 <- vector(mode = "list", length = 8)
spec3[[1]] <- function(u) {
  1
}

test_that("TLSWsim executes with wd spec", {
  x <- TLSWsim(spec = spec1)
  expect_equal(class(x), "numeric")
})

test_that("TLSWsim executes with matrix spec", {
  x <- TLSWsim(spec = spec2)
  expect_equal(class(x), "numeric")
})

test_that("TLSWsim executes with list spec", {
  x <- TLSWsim(spec = spec3)
  expect_equal(class(x), "numeric")
})

test_that("TLSWsim executes with numeric trend", {
  x <- TLSWsim(trend = rep(0, 256), spec = spec1)
  expect_equal(class(x), "numeric")
})

test_that("TLSWsim executes with function trend", {
  x <- TLSWsim(trend = function(u) {
    0
  }, spec = spec1)
  expect_equal(class(x), "numeric")
})

test_that("spectrum is positive", {
  spec1$D[1] <- -10
  expect_error(
    TLSWsim(spec = spec1),
    "All spectral elements must be non-negative."
  )
})

test_that("innov.func argument is a function", {
  expect_error(
    TLSWsim(spec = spec1, innov.func = 1),
    "Argument'innov.func' should be a function."
  )
})

test_that("innov.func argument is an rnorm type function", {
  expect_error(
    TLSWsim(spec = spec1, innov.func = function(u) {
      1
    }),
    "Invalid 'innov.func' argument: should be in the rnorm family of functions."
  )
})


test_that("spec matrix dimensions match", {
  expect_error(
    TLSWsim(spec = spec2[1:4, ]),
    "Dimensions of spec matrix incorrect. The integer part of log2\\(number of columns\\) should equal the number of rows."
  )
})

test_that("trend and spec dimensions match", {
  expect_error(
    TLSWsim(trend = rep(0, 10), spec = spec1),
    "Length of trend does not match dimensions of spec."
  )
})

test_that("TLSWsim gives same output", {
  skip_on_cran()
  spec <- matrix(0, nrow = 9, ncol = 512)
  spec[1, ] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
  trend <- seq(from = 0, to = 5, length = 512)
  set.seed(1)
  x <- TLSWsim(trend = trend, spec = spec)
  expect_snapshot_value(round(x, digits = 3), style = "deparse")
})

test_that("TLSWsim works on non-dyadic data", {
  skip_on_cran()
  expect_equal(200, length(TLSWsim(rep(0,200), spec2[1:7,1:200])))
})

test_that("TLSWsim rejects non-dyadic spec and function trend", {
  skip_on_cran()
  expect_error(TLSWsim(function(z){0}, spec2[1:7,1:200]),
               "If spec has a non-dyadic number of columns, then trend must be a numeric vector.")
})

