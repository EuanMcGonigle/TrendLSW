test_that("summary TLSW gives same output",{
  skip_on_cran()
  spec <- matrix(0, nrow = 9, ncol = 512)
  spec[1,] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
  trend <- seq(from = 0, to = 5, length = 512)
  set.seed(1)
  x <- TLSWsim(trend = trend, spec = spec)
  x.TLSW <- TLSW(x)
  expect_snapshot(summary(x.TLSW))
})

test_that("print TLSW gives same output",{
  skip_on_cran()
  spec <- matrix(0, nrow = 9, ncol = 512)
  spec[1,] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
  trend <- seq(from = 0, to = 5, length = 512)
  set.seed(1)
  x <- TLSWsim(trend = trend, spec = spec)
  x.TLSW <- TLSW(x)
  expect_snapshot_value(print(x.TLSW))
})

test_that("plot TLSW gives same output",{
  skip_on_cran()
  spec <- matrix(0, nrow = 9, ncol = 512)
  spec[1,] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
  trend <- seq(from = 0, to = 5, length = 512)
  set.seed(1)
  x <- TLSWsim(trend = trend, spec = spec)
  x.TLSW <- TLSW(x)
  expect_snapshot_file(pdf(file='plotTLSW.pdf');plot(x.TLSW);dev.off)
})
