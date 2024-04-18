#' @title Linear Wavelet Thresholding Trend Estimation of Time Series
#' @description Internal function to compute the linear wavelet thresholding trend estimate for a
#' time series that may be second-order nonstationary. The function calculates
#' the wavelet transform of the time series, sets to zero the non-boundary coefficients,
#' then inverts the transform to obtain the estimate.
#' This function is not intended for general use by regular users of the package.
#' @param x The time series you want to estimate the trend function of.
#' @param filter.number Selects the index of the wavelet used in the estimation
#' procedure. For Daubechies compactly supported wavelets the filter number is
#' the number of vanishing moments.
#' @param family Selects the wavelet family to use. Recommended to only use the
#' Daubechies compactly supported wavelets DaubExPhase and DaubLeAsymm.
#' @param max.scale Selects the coarsest scale of the wavelet transform to
#' analyse to. Should be a value from \eqn{1} (finest) to \eqn{J-1} (coarsest),
#' where \eqn{n=2^J} is the length of the time series.
#' @param transform.type The type of wavelet transform used. Can be \code{"dec"}
#' which is the standard discrete wavelet transform or \code{"nondec"} (default),
#' a non-decimated wavelet transform, but a confidence interval
#' cannot be calculated in this case.
#' @param boundary.handle Logical variable. If \code{TRUE}, the time
#' series is boundary corrected, to get a less variable trend estimate at the
#' boundaries of the times series. If \code{FALSE}, no boundary correction is applied.
#' @param T.CI Logical variable, only to be used if \code{transform.type = TRUE}.
#' If \code{TRUE}, a \code{(1-sig.lvl)} pointwise confidence interval is
#' computed for the trend estimate.
#' @param sig.lvl Used only if \code{T.CI = TRUE}; a numeric value
#' (\code{0 <= sig.lvl <= 1}) with which a \code{(1-sig.lvl)} pointwise
#' confidence interval for the trend estimate is generated.
#' @param lag.max Used only if \code{T.CI = TRUE}; a positive integer
#' specifying the maximum lag to which the local autocovariance function is
#' estimated.
#' @param confint.type Used only if \code{T.CI = TRUE}; the type of confidence
#' interval computed. Can be \code{"percentile"}, in which case empirical percentiles are used, or
#' \code{"normal"} (default), in which case the normal approximation is used.
#' @param reps Used only if \code{T.CI = TRUE} and \code{transform.type = "dec"} ; the number
#' of bootstrap replications used to compute the confidence interval.
#' @param spec.est Used only if \code{T.CI = TRUE}; the spectrum estimate of the time series,
#' used to calculate the confidence interval for the trend estimate.
#' @param ...  Further arguments to be passed to the \code{\link{ewspec.trend}} call.
#' @return A \code{list} object containing the following fields:
#' \item{x}{Input data}
#' \item{filter.number, family}{Input wavelet parameters}
#' \item{transform.type, max.scale, boundary.handle, T.CI}{Input parameters}
#' \item{T}{A vector of length \code{length(x)} containing the trend estimate}
#' \item{lower.CI}{Returned if \code{T.CI = TRUE}. The lower limit of the pointwise confidence interval}
#' \item{upper.CI}{Returned if \code{T.CI = TRUE}. The upper limit of the pointwise confidence interval}
#' \item{sig.lvl}{Returned if \code{T.CI = TRUE}. The significance level of the pointwise confidence interval}
#' @seealso \code{\link{TLSW}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @keywords internal
wav.trend.est <- function(x, filter.number = 4, family = "DaubLeAsymm",
                          max.scale = floor(log2(length(x)) * 0.7),
                          transform.type = "nondec",
                          boundary.handle = FALSE, T.CI = FALSE, sig.lvl = 0.05,
                          lag.max = floor(10 * (log10(length(x)))),
                          confint.type = "normal",
                          reps = 199, spec.est = NULL, ...) {
  # this function carries out wavelet thresholding of a time series to obtain a
  # trend estimate. All non-boundary wavelet coefficients up to a specified scale
  # are set to zero.

  x.check <- trend.est.checks(
    x = x, max.scale = max.scale, boundary.handle = boundary.handle,
    transform.type = transform.type, T.CI = T.CI,
    reps = 2, sig.lvl = sig.lvl, est.type = "linear"
  )

  x.len <- x.check$x.len
  max.scale <- x.check$max.scale
  boundary.handle <- x.check$boundary.handle
  J <- x.check$J
  dyadic <- x.check$dyadic

  orig.x <- x
  if (boundary.handle == TRUE) {
    x <- get.boundary.timeseries(x)
  }
  x.len <- length(x)
  J <- wavethresh::IsPowerOfTwo(x.len)


  # below code determines the boundary coefficients for a given wavelet

  boundary.test <- c(rep(0, x.len - 1), 1)

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
    x.wd <- wavethresh::wd(x, filter.number = filter.number, family = family)
  } else if (transform.type == "nondec") {
    x.wd <- wavethresh::wd(x, filter.number = filter.number, family = family, type = "station")
  }

  x.thresh <- x.wd

  # set to zero the non-boundary coefficients

  for (j in (J - 1):(J - max.scale)) {
    temp <- wavethresh::accessD(x.wd, level = j)

    temp[-boundary.coefs[[j]]] <- 0

    x.thresh <- wavethresh::putD(x.thresh, temp, level = j)
  }

  # perform inverse transform on thresholded coefficients
  if (transform.type == "dec") {
    x_wr <- wavethresh::wr(x.thresh)
  } else if (transform.type == "nondec") {
    x_wr <- wavethresh::AvBasis(wavethresh::convert(x.thresh))
  }

  # subset the longer estimate to get the true estimate


  if (T.CI == FALSE) {
    if (boundary.handle == TRUE) {
      if (dyadic == TRUE) {
        lower <- 2^(J - 2) + 2^(J - 3) + 1
        upper <- 2^(J - 1) + 2^(J - 3)
      } else {
        lower <- floor((x.len - length(orig.x)) / 2)
        upper <- lower + length(orig.x) - 1
      }
      x_wr <- x_wr[lower:upper]
    }
    return(list(
      x = orig.x, T = x_wr, filter.number = filter.number, family = family,
      transform.type = transform.type, max.scale = max.scale,
      boundary.handle = boundary.handle, T.CI = T.CI
    ))
  } else {
    if (transform.type == "dec") {
      if (!is.null(spec.est)) {
        if (boundary.handle == TRUE) {
          spec.est2 <- ewspec.trend(
            x = x, an.filter.number = spec.est$S$filter$filter.number,
            an.family = spec.est$S$filter$family,
            binwidth = spec.est$binwidth,
            max.scale = spec.est$max.scale, S.smooth = spec.est$S.smooth,
            smooth.type = spec.est$smooth.type,
            AutoReflect = FALSE, boundary.handle = FALSE
          )
        } else {
          spec.est2 <- spec.est
        }
      } else {
        spec.est2 <- ewspec.trend(x, max.scale = max.scale, ..., AutoReflect = FALSE)
      }


      lacf.est <- TLSW.TLSWlacf(x,
        filter.number = spec.est2$S$filter$filter.number, family = spec.est2$S$filter$family,
        lag.max = lag.max, spec.est = spec.est2
      )

      trend.CI <- trend.estCI(
        trend.est = x_wr, lacf.est = lacf.est, filter.number = filter.number,
        family = family, sig.lvl = sig.lvl
      )
      lower.CI <- trend.CI$lower.CI
      upper.CI <- trend.CI$upper.CI

      if (boundary.handle == TRUE) {
        if (dyadic == TRUE) {
          lower <- 2^(J - 2) + 2^(J - 3) + 1
          upper <- 2^(J - 1) + 2^(J - 3)
        } else {
          lower <- floor((x.len - length(orig.x)) / 2)
          upper <- lower + length(orig.x) - 1
        }
        x_wr <- x_wr[lower:upper]
        lower.CI <- lower.CI[lower:upper]
        upper.CI <- upper.CI[lower:upper]
      }
    } else {
      if (boundary.handle == TRUE) {
        if (dyadic == TRUE) {
          lower <- 2^(J - 2) + 2^(J - 3) + 1
          upper <- 2^(J - 1) + 2^(J - 3)
        } else {
          lower <- floor((x.len - length(orig.x)) / 2)
          upper <- lower + length(orig.x) - 1
        }
        x_wr <- x_wr[lower:upper]
      }

      trend.CI <- trend.est.CI.bootstrap(
        x = orig.x,
        trend.est = x_wr, spec.est = spec.est, filter.number = filter.number,
        family = family, max.scale = max.scale, boundary.handle = boundary.handle,
        reps = reps, sig.lvl = sig.lvl, confint.type = confint.type,
        diff = FALSE
      )
      lower.CI <- trend.CI[1, ]
      upper.CI <- trend.CI[2, ]
    }

    return(list(
      x = orig.x, T = x_wr, lower.CI = lower.CI, upper.CI = upper.CI,
      sig.lvl = sig.lvl, filter.number = filter.number, family = family, transform.type = transform.type,
      max.scale = max.scale, boundary.handle = boundary.handle, T.CI = T.CI
    ))
  }
}
