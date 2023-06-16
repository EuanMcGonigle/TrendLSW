test_that("Atau.mat.calc executes", {
  skip_on_cran()
  J <- 5
  A.J <- Atau.mat.calc(J = J)
  expect_equal(class(A.J), c("matrix", "array"))
  expect_equal(ncol(A.J), J)
  expect_equal(nrow(A.J), J)
})

test_that("positive J required", {
  expect_error(
    Atau.mat.calc(J = -5),
    "The parameter J should be a positive integer."
  )
})

test_that("numeric J required", {
  expect_error(
    Atau.mat.calc(J = "5"),
    "The parameter J should be a positive integer."
  )
})

test_that("integer J required", {
  expect_error(
    Atau.mat.calc(J = "5.5"),
    "The parameter J should be a positive integer."
  )
})

test_that("positive lag required", {
  expect_error(
    Atau.mat.calc(lag = 0),
    "The lag parameter should be a positive integer."
  )
})

test_that("numeric lag required", {
  expect_error(
    Atau.mat.calc(lag = "1"),
    "The lag parameter should be a positive integer."
  )
})

test_that("integer lag required", {
  expect_error(
    Atau.mat.calc(lag = 4.5),
    "The lag parameter should be a positive integer."
  )
})
