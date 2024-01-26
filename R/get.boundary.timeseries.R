#' @title Calculate Boundary Extended Time Series
#' @description Internal function to calculate the boundary extended time series, to be used within
#' the \code{\link{TLSW}} function. Not recommended for general usage.
#' @param x The time series used to calculate the boundary extended version.
#' @param type The type of boundary handling used. Either "TLSW" for boundary
#' handling as described in McGonigle, E. T., Killick, R., and Nunes, M. (2022a),
#' or "LSW.diff" for a periodic version of this, used for the differencing-based functions
#' as described in McGonigle, E. T., Killick, R., and Nunes, M. (2022b).
#' @return A vector of 4 times the length of the input vector.
#' @seealso \code{\link{TLSW}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022a). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' #' McGonigle, E. T., Killick, R., and Nunes, M. (2022b). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' @keywords internal
get.boundary.timeseries <- function(x, type = c("TLSW", "LSW.diff")[1]) {
  stopifnot("Error: type of boundary handling must be either 'TLSW' or
            'LSW.diff'." = type == "TLSW" || type == "LSW.diff")


  x.len <- length(x)
  J <- wavethresh::IsPowerOfTwo(x.len)

  s <- seq(from = 0, to = (x.len - 1) / x.len, length = x.len)

  L <- stats::lm(x ~ stats::poly(s, 3, raw = TRUE))

  bh.right <- stats::predict(L, newdata = data.frame(s = 1))

  bh.left <- stats::predict(L, newdata = data.frame(s = -1 / x.len))

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
