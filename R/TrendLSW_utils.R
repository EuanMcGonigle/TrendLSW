#' @title Shift a Vector
#' @description Internal function for shifting a vector
#' @keywords internal
#' @noRd

do.shift <- function(v, places, dir = "right") {
  vnew <- NULL
  d <- substring(dir, 1, 1)
  if (d == "r" & places < 0) {
    d <- "l"
  } else {
    if (d == "l" & places < 0) {
      d <- "r"
    }
  }
  n <- length(v)
  if (n == 1) {
    places <- 0
  }
  p <- abs(places)
  if (p == 0) {
    vnew <- v
  } else {
    if (d == "r") {
      vnew <- c(v[(n - p + 1):n], v[1:(n - p)])
    } else {
      vnew <- c(v[(p + 1):n], v[1:p])
    }
  }
  vnew
}

#' @title Compute Covariance Matrix
#' @description Internal function for computing covariance matrix
#' @keywords internal
#' @noRd


create.covmat <- function(lacf, data.len) {
  lacf <- lacf$lacf

  max.lag <- length(lacf[1, ])

  cov.mat <- matrix(0, nrow = data.len, ncol = data.len)

  for (row in 1:data.len) {
    for (column in row:(min((max.lag - 1 + row), data.len))) {
      cov.mat[row, column] <- lacf[row, (abs(column - row) + 1)]
      cov.mat[column, row] <- cov.mat[row, column]
    }
  }

  cov.mat
}



