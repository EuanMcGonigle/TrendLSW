test_that("Cmat.calc executes with larger gen wavelet", {
  skip_on_cran()
  J <- 5
  C <- Cmat.calc(J = J, gen.filter.number = 2)
  expect_equal(class(C), c("matrix", "array"))
  expect_equal(ncol(C), J)
  expect_equal(nrow(C), J)
})

test_that("Cmat.calc executes with larger an wavelet", {
  skip_on_cran()
  J <- 5
  C <- Cmat.calc(J = J, an.filter.number = 2)
  expect_equal(class(C), c("matrix", "array"))
  expect_equal(ncol(C), J)
  expect_equal(nrow(C), J)
})

test_that("positive J required", {
  expect_error(
    Cmat.calc(J = -5),
    "The parameter J should be a positive integer."
  )
})

test_that("numeric J required", {
  expect_error(
    Cmat.calc(J = "5"),
    "The parameter J should be a positive integer."
  )
})

test_that("integer J required", {
  expect_error(
    Cmat.calc(J = "5.5"),
    "The parameter J should be a positive integer."
  )
})
