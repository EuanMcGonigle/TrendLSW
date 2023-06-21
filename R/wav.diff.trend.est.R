#' @title Wavelet Thresholding Trend Estimation of Time Series
#' @description Computes the wavelet thresholding trend estimate for a time
#' series that may be second-order nonstationary. The function calculates the
#' wavelet transform of the time series, thresholds the coefficients based on
#' an estimate of their variance, and inverts to give the trend estimate.
#' @details Estimates the trend function of a locally stationary time series, by
#' incorporating the evolutionary wavelet spectrum estimate in a wavelet
#' thresholding procedure. To use this function, first compute the spectral
#' estimate of the time series, using the function ewspec.diff.
#'
#' The function works as follows:
#'
#' 1. The wavelet transform of the time series is calculated.
#'
#' 2. The wavelet coefficients are individually thresholded using the universal
#' threshold \eqn{\hat{\sigma}\sqrt(2 log T)}, where \eqn{\hat{sigma}^2} is an estimate of their variance. The variance
#' estimate is calculated using the spectral estimate, supplied by the user in
#' the \code{spec} argument.
#'
#' 3. The inverse wavelet transform is applied to obtain the final estimate.
#' @param data The time series you want to estimate the trend function of.
#' @param spec.est You must supply the estimate of the evolutionary wavelet
#' spectrum of the time series. This is the output of the \code{ewspec.diff}
#' function.
#' @param filter.number Selects the index of the wavelet used in the estimation
#' procedure. For Daubechies compactly supported wavelets the filter number is
#' the number of vanishing moments.
#' @param thresh.type The type of thresholding function used. Currently only
#' "soft" and "hard" are available. Recommended to use "soft".
#' @param normal If TRUE, uses a threshold assuming the data are normally
#' distributed. If FALSE, uses a larger threshold to reflect non-normality.
#' @param family Selects the wavelet family to use. Recommended to only use the
#' Daubechies compactly supported wavelets DaubExPhase and DaubLeAsymm.
#' @param max.scale Selects the number of scales of the wavelet transform to
#' apply thresholding to. Should be a value from 1 (finest) to J-1 (coarsest),
#' where T=2^J is the length of the time series. Recommended to use 2J/3
#' scales.
#' @param boundary.handle Logical variable, decides if boundary handling should
#' be applied to the time series before estimation.
#' @param calc.confint Logical variable. If \code{TRUE}, a bootstrapped \code{(1-sig.lvl)}
#' pointwise confidence interval is computed for the trend estimate.
#' @param sig.lvl Used only if \code{calc.confint = TRUE}; a numeric value
#' (\code{0 <= sig.lvl <= 1}) with which a \code{(1-sig.lvl)} pointwise
#' confidence interval for the trend estimate is generated.
#' @param reps Used only if \code{calc.confint = TRUE}; the number of bootstrap
#' replications used to calcualte the confidence interval.
#' @param ...  Further arguments to be passed to the \code{\link{ewspec.diff}}
#' call, only to be used if \code{calc.confint = TRUE}.
#' @return A vector of length \code{length(data)} containing the trend estimate.
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link{ewspec.diff}}, \code{\link{wav.trend.est}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' @examples
#' spec <- wavethresh::cns(1024, filter.number = 4)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 1024))^2)
#'
#' set.seed(120)
#'
#' noise <- wavethresh::LSWsim(spec)
#' sine_trend <- -2 * sin(seq(from = 0, to = 2 * pi, length = 1024)) -
#'   1.5 * cos(seq(from = 0, to = pi, length = 1024))
#'
#' x <- sine_trend + noise
#'
#' spec.est <- ewspec.diff(data = x, family = "DaubExPhase", filter.number = 4, max.scale = 7)
#'
#' trend.est <- wav.diff.trend.est(data = x, spec = spec.est)
#'
#' plot.ts(x, lty = 1, col = 8)
#' lines(sine_trend, col = 2, lwd = 2)
#' lines(trend.est, col = 4, lwd = 2, lty = 2)
#' @export
wav.diff.trend.est <- function(data, spec.est, filter.number = 4, thresh.type = "soft",
                               normal = TRUE, family = "DaubLeAsymm",
                               max.scale = floor(0.7 * log2(length(data))),
                               boundary.handle = FALSE, calc.confint = FALSE,
                               reps = 199, sig.lvl = 0.05, ...) {
  data.check <- ewspec.checks(
    data = data, max.scale = max.scale, lag = 1,
    binwidth = 1, boundary.handle = boundary.handle
  )

  data.len <- data.check$data.len
  max.scale <- data.check$max.scale
  boundary.handle <- data.check$boundary.handle
  J <- data.check$J
  dyadic <- data.check$dyadic

  spec <- spec.est$S
  # by default, do T.I. denoising:

  if (boundary.handle == TRUE) {
    orig.data <- data
    data <- get.boundary.timeseries(data)
  }

  data.len <- length(data)
  J <- wavethresh::IsPowerOfTwo(data.len)


  data.wd <- wavethresh::wd(data, filter.number = filter.number, family = family, type = "station")

  # calculate C to use in variance estimate of wavelet coefficients

  C <- Cmat.calc(
    J = max.scale, an.filter.number = filter.number,
    gen.filter.number = spec$filter$filter.number,
    an.family = family, gen.family = spec$filter$family
  )

  spec.mat <- matrix(0, nrow = max.scale, ncol = length(wavethresh::accessD(spec, level = 0)))

  for (j in 1:max.scale) {
    spec.mat[j, ] <- wavethresh::accessD(spec, level = wavethresh::nlevelsWT(spec) - j)
  }

  # calculate variance estimate matrix

  var.mat <- C %*% spec.mat

  var.mat <- replace.neg.values(var.mat, max.scale)

  data.wd.thresh <- data.wd


  if (boundary.handle == TRUE) {
    # below code determines the boundary coefficients for a given wavelet

    boundary.test <- c(rep(0, data.len - 1), 1)

    y.wd <- wavethresh::wd(boundary.test, family = family, filter.number = filter.number, type = "station")

    # create EWS estimate matrix

    bc.var.mat <- matrix(0, nrow = max.scale, ncol = data.len)

    lower <- floor((data.len - length(orig.data)) / 2)
    upper <- lower + length(orig.data) - 1

    bc.var.mat[, lower:upper] <- var.mat[, 1:length(orig.data)]


    for (j in 1:max.scale) {
      dj <- wavethresh::accessD(data.wd, level = J - j)

      if (normal == TRUE) {
        thresh <- sqrt(2 * bc.var.mat[j, ] * log(data.len))
      } else {
        thresh <- sqrt(bc.var.mat[j, ]) * log(data.len)
      }

      temp1 <- dj

      temp1[abs(dj) < thresh] <- 0

      if (thresh.type == "soft") {
        temp2 <- dj[abs(dj) >= thresh]
        temp1[abs(dj) >= thresh] <- sign(temp2) * (abs(temp2) - thresh[abs(dj) >= thresh])
      }

      data.wd.thresh <- wavethresh::putD(data.wd.thresh, level = J - j, temp1)
    }

    # invert the thresholded object to get trend estimate:

    trend.est <- wavethresh::AvBasis(wavethresh::convert(data.wd.thresh))

    if (boundary.handle == TRUE) {
      if (dyadic == TRUE) {
        lower <- 2^(J - 2) + 2^(J - 3) + 1
        upper <- 2^(J - 1) + 2^(J - 3)
      } else {
        lower <- floor((data.len - length(orig.data)) / 2)
        upper <- lower + length(orig.data) - 1
      }
      trend.est <- trend.est[lower:upper]
    }
  } else {
    # threshold the wavelet coefficients using the variance estimate and user
    # inputted rules

    for (j in 1:max.scale) {
      dj <- wavethresh::accessD(data.wd, level = J - j)

      if (normal == TRUE) {
        thresh <- sqrt(2 * var.mat[j, ] * log(data.len))
      } else {
        thresh <- sqrt(var.mat[j, ]) * log(data.len)
      }

      temp1 <- dj

      temp1[abs(dj) < thresh] <- 0

      if (thresh.type == "soft") {
        temp2 <- dj[abs(dj) >= thresh]
        temp1[abs(dj) >= thresh] <- sign(temp2) * (abs(temp2) - thresh[abs(dj) >= thresh])
      }

      data.wd.thresh <- wavethresh::putD(data.wd.thresh, level = J - j, temp1)
    }

    # invert the thresholded object to get trend estimate:

    trend.est <- wavethresh::AvBasis(wavethresh::convert(data.wd.thresh))
  }

  if (calc.confint == TRUE) {
    trend.confint <- trend.estCI.diff(
      data = orig.data, trend.est = trend.est, spec.est = spec.est,
      filter.number = filter.number, thresh.type = thresh.type,
      normal = normal, boundary.handle = boundary.handle,
      family = family, max.scale = max.scale,
      reps = reps, sig.lvl = sig.lvl, ...
    )
    return(list(trend.est = trend.est, conf.int = trend.confint))
  } else {
    return(trend.est)
  }
}
