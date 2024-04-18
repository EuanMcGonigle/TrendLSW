#' @title Estimation of Evolutionary Wavelet Spectrum of Non-Zero Mean Time
#' Series via Differencing
#' @description Internal function to estimate the evolutionary wavelet spectrum (EWS) of
#' a time series that may include a trend component. The estimate is computed
#' by taking the non-decimated wavelet transform of the first differenced time
#' series data, squaring it; smoothing using a running mean and then correcting
#' for bias using the appropriate correction matrix. Inherits the smoothing
#' functionality from the \code{ewspec3} function in the R package \code{locits}.
#' This function is not intended for general use by regular users of the package.
#' @details Computes an estimate of the evolutionary wavelet spectrum of a
#' time series that displays nonstationary mean and autocovariance. The
#' estimation procedure is as follows:
#'
#' 1. The time series is differenced to remove the trend.
#'
#' 2. The squared modulus of the non-decimated wavelet transform is computed,
#' known as the raw wavelet periodogram. This is returned by the function.
#'
#' 3. The raw wavelet periodogram is smoothed using a running mean smoother.
#'
#' 4. The smoothed periodogram is bias corrected using the inverse of the bias
#' matrix.
#'
#' The final estimate, stored in the S component, can be plotted using the plot
#' function, please see the example below.
#' @param x The time series you wish to analyse.
#' @param lag An integer indicating which lag to use for differencing.
#' @param filter.number The index number for the wavelet used to analyse the
#' time series. For the "DaubExPhase" family, the filter number can be between
#' 1 to 10. For the "DaubLeAsymm" family, the filter number can be between 4 to
#' 10.
#' @param family The family of the wavelet used. It is recommended to use
#' either the Daubechies Extremal Phase family, or the Daubechies Least
#' Asymmetric family, corresponding to the "DaubExPhase" or the "DaubLeAsymm"
#' options.
#' @param binwidth The bin width of the running mean smoother used to smooth
#' the raw wavelet periodogram.
#' @param diff.number The number of differences used to remove the trend of the
#' series. A first difference is recommended as default.
#' @param max.scale The coarsest level to which the time series is analysed to.
#' Should be a positive integer less than \eqn{J}, where \eqn{T=2^J} is the length of the
#' time series. The default setting is \eqn{0.7J}, to control for bias.
#' @param S.smooth Argument that dictates if smoothing is performed on the raw
#' wavelet periodogram.
#' @param smooth.type String indicating which type of smoothing to use on wavelet periodogram.
#' Can be \code{"mean"}, \code{"median"}, or \code{"epan"}. Default is \code{"epan"}.
#' @param boundary.handle Logical variable, if TRUE, then boundary handling
#' will be applied when computing the periodogram. Recommended to set as FALSE,
#' will be set as TRUE automatically if non-dyadic data is used.
#' @param AutoReflect As in \code{wavethresh}. Decides whether or not the time series
#' is reflected when computing the wavelet transform. Strongly recommended to
#' leave as FALSE.
#' @param supply.inv.mat Not intended for general use. If TRUE, user must supply the
#' appropriate correction matrix.
#' @param inv.mat Not intended for general use. If supply.mat is TRUE, user must supply the
#' appropriate correction matrix used to correct the raw wavelet periodogram.
#' Equal to \eqn{(2A-2A_1)^{-1}} for first differences.
#' @return A list object, containing the following fields:
#'  \item{S}{The evolutionary wavelet spectral estimate of the input data. This object is of
#' class wd and so can be plotted and printed in the usual way using wavethresh
#' functionality. }
#' \item{WavPer}{ The raw wavelet periodogram of the input
#' data. The EWS estimate (above) is the smoothed corrected version of the
#' wavelet periodogram.}
#' \item{SmoothWavPer}{ The smoothed, un-corrected raw
#' wavelet periodogram of the input data. }
#' \item{max.scale, boundary.handle, S.smooth, smooth.type, binwidth, lag, diff.number}{Input parameters}
#' @seealso \code{\link{TLSW}}, \code{\link[locits]{ewspec3}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' wavethresh::plot.wd(spec.est$S)
#' @keywords internal
ewspec.diff <- function(x, lag = 1, filter.number = 4, family = "DaubExPhase",
                        binwidth = floor(2 * sqrt(length(x))), diff.number = 1,
                        max.scale = floor(log2(length(x)) * 0.7), S.smooth = TRUE,
                        smooth.type = "epan",
                        boundary.handle = FALSE, AutoReflect = FALSE,
                        supply.inv.mat = FALSE, inv.mat = NULL) {
  # function that computes the spectral estimate of a time series that has a trend.

  x.check <- ewspec.checks(
    x = x, max.scale = max.scale, lag = lag,
    binwidth = binwidth, boundary.handle = boundary.handle,
    S.smooth = S.smooth, smooth.type = smooth.type
  )

  x.len <- x.check$x.len
  max.scale <- x.check$max.scale
  boundary.handle <- x.check$boundary.handle
  J <- x.check$J
  dyadic <- x.check$dyadic

  if (boundary.handle == TRUE) {
    x <- get.boundary.timeseries(x, type = "LSW.diff")
  }

  # difference x to remove trend/seasonality and calculate the appropriate correction matrix
  # and its inverse:
  if (diff.number != 1 && diff.number != 2) {
    diff.number <- 1
    warning("Function only implements 1st or 2nd differences. Using 1st difference.")
  }

  if (diff.number == 1) {
    diff.x <- c(diff(x, lag))
    diff.x <- c(diff.x, rep(0, lag))
  } else if (diff.number == 2) {
    diff.x <- c(diff(diff(x)), 0, 0)
    if (lag != 1) {
      lag <- 1
      warning("When diff.number = 1, only lag = 1 is supported. Resetting lag to be 1.")
    }
  }

  if (supply.inv.mat == FALSE) {
    A <- wavethresh::ipndacw(J = -max.scale, filter.number = filter.number, family = family)
    A1 <- Atau.mat.calc(J = max.scale, filter.number = filter.number, family = family, lag = lag)

    if (diff.number == 1) {
      inv.mat <- solve(2 * A - 2 * A1)
    } else if (diff.number == 2) {
      A2 <- Atau.mat.calc(J = max.scale, filter.number = filter.number, family = family, lag = 2)
      inv.mat <- solve(6 * A - 8 * A1 + 2 * A2)
    }
  } else {
    inv.mat <- supply.mat.check(inv.mat = inv.mat, max.scale = max.scale)
  }

  # calculate raw wavelet periodogram which we need to correct:
  if (smooth.type == "median" || smooth.type == "epan") {
    x.wd <- locits::ewspec3(diff.x,
      filter.number = filter.number, family = family,
      binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = FALSE
    )
    x.wd <- WP.manual.smooth(
      x.wd = x.wd, smooth.type = smooth.type,
      max.scale = max.scale, binwidth = binwidth
    )
  } else {
    x.wd <- locits::ewspec3(diff.x,
      filter.number = filter.number, family = family,
      binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = S.smooth
    )
  }

  if (boundary.handle == TRUE) {
    x.wd <- smooth.wav.per.calc(
      x.wd = x.wd, J = J, x.len = x.len,
      filter.number = filter.number, family = family,
      dyadic = dyadic, max.scale = max.scale
    )
  }


  l <- S.calc(
    x.wd = x.wd, max.scale = max.scale, J = J, inv.mat = inv.mat,
    filter.number = filter.number, family = family
  )

  l$max.scale <- max.scale
  l$boundary.handle <- boundary.handle
  l$S.smooth <- S.smooth
  if (S.smooth == TRUE) {
    l$smooth.type <- smooth.type
    l$binwidth <- binwidth
  }
  l$lag <- lag
  l$diff.number <- diff.number

  return(l)
}
