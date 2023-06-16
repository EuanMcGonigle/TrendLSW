set.seed(1)
x1 <- rnorm(400) + seq(from = 0, to = 2, length = 400)
x2 <- rnorm(512) + seq(from = 0, to = 2, length = 512)

test_that("boundary handling type is recognised", {
  expect_error(
    get.boundary.timeseries(x, type = "symmetric"),
    "Error: type of boundary handling must be either 'TLSW' or
            'LSW.diff'."
  )
})

test_that("get.boundary.timeseries executes on dyadic data with param 'type' = 'TLSW'", {
  skip_on_cran()
  expect_equal(class(get.boundary.timeseries(x1, type = "TLSW")), "numeric")
})

test_that("get.boundary.timeseries executes on non-dyadic data with param 'type' = 'TLSW'", {
  skip_on_cran()
  expect_equal(class(get.boundary.timeseries(x2, type = "TLSW")), "numeric")
})

test_that("get.boundary.timeseries executes on dyadic data with param 'type' = 'LSW.diff'", {
  skip_on_cran()
  expect_equal(class(get.boundary.timeseries(x1, type = "LSW.diff")), "numeric")
})

test_that("get.boundary.timeseries executes on non-dyadic data with param 'type' = 'LSW.diff'", {
  skip_on_cran()
  expect_equal(class(get.boundary.timeseries(x2, type = "LSW.diff")), "numeric")
})

test_that("get.boundary.timeseries executes on odd-lengthed data", {
  skip_on_cran()
  set.seed(2)
  x.odd <- rnorm(373)
  expect_equal(class(get.boundary.timeseries(x.odd)), "numeric")
})
