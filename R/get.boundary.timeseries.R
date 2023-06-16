#' @title Calculates Boundary Extended Time Series
#' @description A function to calculate the boundary extended time series, to be used within
#' the \code{\link{ewspec.trend}}, \code{\link{ewspec.diff}},
#' \code{\link{wav.trend.est}}, and \code{\link{wav.diff.trend.est}} functions.
#' Not recommended for general usage.
#' @param data The time series used to calculate the boundary extended version.
#' @param type The type of boundary handling used. Either "TLSW" for boundary
#' handling as described in McGonigle, E. T., Killick, R., and Nunes, M. (2022)
#' Trend locally stationary wavelet processes, or "LSW.diff" for a periodic
#' version of this, used for the differencing-based functions.
#' @return A vector of 4 times the length of the input vector.
#' @seealso \code{\link{ewspec.trend}}, \code{\link{ewspec.diff}},
#' \code{\link{wav.trend.est}}, \code{\link{wav.diff.trend.est}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @examples
#' trend <- seq(from = 0, to = 2, length = 400)
#' set.seed(1)
#'
#' x <- rnorm(400) + trend
#'
#' x.b <- get.boundary.timeseries(x)
#'
#' plot.ts(x.b)
#' @export
get.boundary.timeseries <- function(data, type = c("TLSW", "LSW.diff")[1]) {
  stopifnot("Error: type of boundary handling must be either 'TLSW' or
            'LSW.diff'." = type == "TLSW" || type == "LSW.diff")


  data.len <- length(data)
  J <- wavethresh::IsPowerOfTwo(data.len)

  s <- seq(from = 0, to = (data.len - 1) / data.len, length = data.len)

  L <- stats::lm(data ~ poly(s, 3, raw = TRUE))

  bh.right <- stats::predict(L, newdata = data.frame(s = 1))

  bh.left <- stats::predict(L, newdata = data.frame(s = -1 / data.len))

  if (type == "LSW.diff") {
    bh.series1 <- c(data - bh.right + bh.left, data, data + bh.right - bh.left)
  } else {
    bh.series1 <- c(-rev(data) + 2 * bh.left, data, -rev(data) + 2 * bh.right)
  }

  if (is.na(J) == TRUE) {
    l <- 2^floor(log2(length(bh.series1)))
    k <- floor((3 * data.len - l) / 2)
    if (data.len %% 2 == 0) {
      bh.series2 <- bh.series1[(k + 1):(3 * data.len - k)]
    } else {
      bh.series2 <- bh.series1[(k + 1):(3 * data.len - k - 1)]
    }

    return(bh.series2)
  } else {
    l <- length(bh.series1)

    bh.series2 <- c((data)[(data.len / 2 + 1):data.len] - abs(2 * bh.right - 2 * bh.left), bh.series1, (data)[1:(data.len / 2)] + abs(2 * bh.right - 2 * bh.left))

    return(bh.series2)
  }
}
