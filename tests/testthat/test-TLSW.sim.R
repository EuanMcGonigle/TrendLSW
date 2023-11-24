spec <- wavethresh::cns(256)
spec <- wavethresh::putD(spec, level = 6, rep(1, 256))

dists <- c("norm", "pois", "exp", "chisq", "t")

for (i in seq_len(length(dists))) {
  test_that(paste0("TLSWsim with distribution #", i, " executes"), {
    skip_on_cran()
    expect_equal(class(TLSWsim(spec = spec, distribution = dists[i])), "numeric")
  })
}

test_that("distribution is recognised", {
  expect_error(
    TLSWsim(spec = spec, distribution = "cauchy"),
    "Error: distribution must be one of 'norm', 'pois',
               'exp, 'chisq', or 't'."
  )
})

test_that("spectrum is positive", {
  spec$D[1] <- -10
  expect_error(
    TLSWsim(spec = spec),
    "All spectral elements must be non-negative."
  )
})


test_that("df is positive", {
  expect_error(
    TLSWsim(spec = spec, df = -2),
    "parameter df must be positive."
  )
})

test_that("rate is positive", {
  expect_error(
    TLSWsim(spec = spec, rate = -2),
    "The rate parameter must be positive."
  )
})

test_that("df is bigger than 2 for t distribution", {
  expect_error(
    TLSWsim(spec = spec, df = 1, distribution = "t"),
    "for t distribution, parameter df must be greater than 2."
  )
})
