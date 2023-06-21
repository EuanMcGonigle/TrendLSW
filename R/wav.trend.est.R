#' @title Linear Wavelet Thresholding Trend Estimation of Time Series
#' @description Computes the linear wavelet thresholding trend estimate for a
#' time series that may be second-order nonstationary. The function calculates
#' the wavelet transform of the time series, sets to zero the non-boundary
#' @param data The time series you want to estimate the trend function of.
#' @param filter.number Selects the index of the wavelet used in the estimation
#' procedure. For Daubechies compactly supported wavelets the filter number is
#' the number of vanishing moments.
#' @param family Selects the wavelet family to use. Recommended to only use the
#' Daubechies compactly supported wavelets DaubExPhase and DaubLeAsymm.
#' @param max.scale Selects the coarsest scale of the wavelet transform to
#' analyse to. Should be a value from \eqn{1} (finest) to \eqn{J-1} (coarsest),
#' where \eqn{T=2^J} is the length of the time series.
#' @param transform.type The type of wavelet transform used. By default, it is "dec"
#' which is the standard discrete wavelet transform. Can also be "nondec",
#' which uses a non-decimated wavelet transform, but a confidence interval
#' cannot be calculated in this case.
#' @param boundary.handle Logical variable. If \code{TRUE}, the time
#' series is boundary corrected, to get a less variable trend estimate at the
#' boundaries of the times series. If \code{FALSE}, no boundary correction is applied.
#' @param calc.confint Logical variable, only to be used if \code{transform.type = TRUE}.
#' If \code{TRUE}, a \code{(1-sig.lvl)} pointwise confidence interval is
#' computed for the trend estimate.
#' @param sig.lvl Used only if \code{calc.confint = TRUE}; a numeric value
#' (\code{0 <= sig.lvl <= 1}) with which a \code{(1-sig.lvl)} pointwise
#' confidence interval for the trend estimate is generated.
#' @param lag.max Used only if \code{calc.confint = TRUE}; a positive integer
#' specifying the maximum lag to which the local autocovariance function is
#' estimated.
#' @param ...  Further arguments to be passed to the \code{\link{ewspec.trend}} call,
#' only to be used if \code{calc.confint = TRUE}.
#' @return A \code{list} object containing the following fields:
#' \item{data}{Input data}
#' \item{filter.number, family}{Inpute wavelet parameters}
#' \item{trend.est}{A vector of length \code{length(data)} containing the trend estimate}
#' \item{calc.confint}{Input parameter}
#' \item{lower.conf}{Returned if \code{calc.confint = TRUE}. The lower limit of the pointwise confidence interval}
#' \item{upper.conf}{Returned if \code{calc.confint = TRUE}. The upper limit of the pointwise confidence interval}
#' \item{sig.lvl}{Returned if \code{calc.confint = TRUE}. The significance level of the pointwise confidence interval}
#' @seealso \code{\link{wav.diff.trend.est}}, \code{\link{ewspec.trend}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @examples
#' # computes trend estimator of simulated linear trend time series
#'
#' set.seed(1)
#'
#' noise <- rnorm(512)
#' trend <- seq(from = 0, to = 5, length = 512)
#' x <- trend + noise
#'
#' trend.est <- wav.trend.est(x, filter.number = 4, family = "DaubLeAsymm", boundary.handle = TRUE)
#'
#' plot.ts(x, lty = 1, col = 8)
#' lines(trend, col = 2, lwd = 2)
#' lines(trend.est$trend.est, col = 4, lwd = 2, lty = 2)
#' @export
wav.trend.est <- function(data, filter.number = 4, family = "DaubLeAsymm",
                          max.scale = floor(log2(length(data)) * 0.7),
                          transform.type = c("dec", "nondec")[1],
                          boundary.handle = FALSE, calc.confint = FALSE, sig.lvl = 0.05,
                          lag.max = floor(10 * (log10(length(data)))), ...) {
  # this function carries out wavelet thresholding of a time series to obtain a
  # trend estimate. All non-boundary wavelet coefficients up to a specified scale
  # are set to zero.

  data.check <- ewspec.checks(
    data = data, max.scale = max.scale, lag = 1,
    binwidth = 1, boundary.handle = boundary.handle
  )

  data.len <- data.check$data.len
  max.scale <- data.check$max.scale
  boundary.handle <- data.check$boundary.handle
  J <- data.check$J
  dyadic <- data.check$dyadic

  trend.est.check(transform.type = transform.type, calc.confint = calc.confint)

  orig.data <- data
  if (boundary.handle == TRUE) {
    data <- get.boundary.timeseries(data)
  }
  data.len <- length(data)
  J <- wavethresh::IsPowerOfTwo(data.len)


  # below code determines the boundary coefficients for a given wavelet

  boundary.test <- c(rep(0, data.len - 1), 1)

  if (transform.type == "dec") {
    y.wd <- wavethresh::wd(boundary.test, family = family, filter.number = filter.number)
  } else if (transform.type == "nondec") {
    y.wd <- wavethresh::wd(boundary.test, family = family, filter.number = filter.number, type = "station")
  }

  boundary.coefs <- list()

  for (j in (J - 1):(J - max.scale)) {
    temp <- wavethresh::accessD(y.wd, level = j)

    boundary.coefs[[j]] <- (which(temp != 0))
  }

  # perform DWT on series
  if (transform.type == "dec") {
    data.wd <- wavethresh::wd(data, filter.number = filter.number, family = family)
  } else if (transform.type == "nondec") {
    data.wd <- wavethresh::wd(data, filter.number = filter.number, family = family, type = "station")
  }

  data.thresh <- data.wd

  # set to zero the non-boundary coefficients

  for (j in (J - 1):(J - max.scale)) {
    temp <- wavethresh::accessD(data.wd, level = j)

    temp[-boundary.coefs[[j]]] <- 0

    data.thresh <- wavethresh::putD(data.thresh, temp, level = j)
  }

  # perform inverse transform on thresholded coefficients
  if (transform.type == "dec") {
    data_wr <- wavethresh::wr(data.thresh)
  } else if (transform.type == "nondec") {
    data_wr <- wavethresh::AvBasis(wavethresh::convert(data.thresh))
  }

  # subset the longer estimate to get the true estimate


  if (calc.confint == FALSE) {
    if (boundary.handle == TRUE) {
      if (dyadic == TRUE) {
        lower <- 2^(J - 2) + 2^(J - 3) + 1
        upper <- 2^(J - 1) + 2^(J - 3)
      } else {
        lower <- floor((data.len - length(orig.data)) / 2)
        upper <- lower + length(orig.data) - 1
      }
      data_wr <- data_wr[lower:upper]
    }
    return(list(data = orig.data, filter.number = filter.number, family = family, trend.est = data_wr, calc.confint = calc.confint))
  } else {
    spec.est <- ewspec.trend(data, max.scale = max.scale, ..., boundary.handle = FALSE, AutoReflect = FALSE)

    lacf.est <- lacf.calc(data,
      filter.number = spec.est$S$filter$filter.number, family = spec.est$S$filter$family,
      lag.max = lag.max, spec.est = spec.est
    )

    trend.confint <- trend.estCI(
      trend.est = data_wr, lacf.est = lacf.est, filter.number = filter.number,
      family = family, sig.lvl = sig.lvl
    )
    lower.conf <- trend.confint$lower.conf
    upper.conf <- trend.confint$upper.conf

    if (boundary.handle == TRUE) {
      if (dyadic == TRUE) {
        lower <- 2^(J - 2) + 2^(J - 3) + 1
        upper <- 2^(J - 1) + 2^(J - 3)
      } else {
        lower <- floor((data.len - length(orig.data)) / 2)
        upper <- lower + length(orig.data) - 1
      }
      data_wr <- data_wr[lower:upper]
      lower.conf <- lower.conf[lower:upper]
      upper.conf <- upper.conf[lower:upper]
    }

    return(list(
      data = orig.data, filter.number = filter.number, family = family, trend.est = data_wr, calc.confint = calc.confint,
      lower.conf = lower.conf, upper.conf = upper.conf,
      sig.lvl = sig.lvl
    ))
  }
}
