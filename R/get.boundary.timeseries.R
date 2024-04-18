#' @title Calculate Boundary Extended Time Series
#' @description Internal function to calculate the boundary extended time series, to be used within
#' the \code{\link{TLSW}} function. Not recommended for general usage.
#' @param x The time series used to calculate the boundary extended version.
#' @param type The type of boundary handling used. Either \code{"TLSW"} (default) for boundary
#' handling as described in McGonigle, E. T., Killick, R., and Nunes, M. (2022a),
#' or \code{"LSW.diff"} for a periodic version of this, used for the differencing-based functions
#' as described in McGonigle, E. T., Killick, R., and Nunes, M. (2022b).
#' @return A vector of 4 times the length of the input vector.
#' @seealso \code{\link{TLSW}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022a). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022b). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' @keywords internal
get.boundary.timeseries <- function(x, type = "TLSW") {
  stopifnot("Error: type of boundary handling must be either 'TLSW' or
            'LSW.diff'." = type == "TLSW" || type == "LSW.diff")


  x.len <- length(x)
  J <- wavethresh::IsPowerOfTwo(x.len)

  # use a proportion of the data to fit a pre-estimate of the trend, in order to extend the time series:

  s.right <- seq(from = 0, to = (x.len - 1) / x.len, length = x.len)[(floor(19 * x.len / 20)):x.len]

  s.left <- seq(from = 0, to = (x.len - 1) / x.len, length = x.len)[1:(floor(x.len / 20))]

  x.right <- x[(floor(19 * x.len / 20)):x.len]

  x.left <- x[1:(floor(x.len / 20))]

  L.right <- stats::lm(x.right ~ s.right)

  L.left <- stats::lm(x.left ~ s.left)

  bh.right <- L.right$coefficients[1] + L.right$coefficients[2]

  bh.left <- L.left$coefficients[1] - L.left$coefficients[2] / x.len

  if (type == "LSW.diff") {
    bh.series1 <- c(x - bh.right + bh.left, x, x + bh.right - bh.left)
  } else {
    bh.series1 <- c(-rev(x) + 2 * bh.left, x, -rev(x) + 2 * bh.right)
  }

  if (is.na(J) == TRUE) {
    l <- 2^floor(log2(length(bh.series1)))
    k <- floor((3 * x.len - l) / 2)
    if (x.len %% 2 == 0) {
      bh.series2 <- bh.series1[(k + 1):(3 * x.len - k)]
    } else {
      bh.series2 <- bh.series1[(k + 1):(3 * x.len - k - 1)]
    }

    return(bh.series2)
  } else {
    l <- length(bh.series1)

    bh.series2 <- c((x)[(x.len / 2 + 1):x.len] - abs(2 * bh.right - 2 * bh.left), bh.series1, (x)[1:(x.len / 2)] + abs(2 * bh.right - 2 * bh.left))

    return(bh.series2)
  }
}
