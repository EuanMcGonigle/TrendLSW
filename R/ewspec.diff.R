#' @title Estimation of Evolutionary Wavelet Spectrum of Non-Zero Mean Time
#' Series via Differencing.
#' @description Estimates the evolutionary wavelet spectrum (EWS) of
#' a time series that may include a trend component. The estimate is computed
#' by taking the non-decimated wavelet transform of the first differenced time
#' series data, squaring it; smoothing using a running mean and then correcting
#' for bias using the appropriate correction matrix. Inherits the smoothing
#' functionality from the \code{ewspec3} function in the R package \code{locits}.
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
#' @param data The time series you wish to analyse.
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
#' @param WP.smooth Argument that dictates if smoothing is performed on the raw
#' wavelet periodogram.
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
#' \itemize{
#'  \item{S}{The evolutionary wavelet spectral estimate of the input data. This object is of
#' class wd and so can be plotted and printed in the usual way using wavethresh
#' functionality. }
#' \item{WavPer}{ The raw wavelet periodogram of the input
#' data. The EWS estimate (above) is the smoothed corrected version of the
#' wavelet periodogram.}
#' \item{SmoothWavPer}{ The smoothed, un-corrected raw
#' wavelet periodogram of the input data. }
#' }
#' @seealso \code{\link{ewspec}}, \code{\link{ewspec3}},
#' \code{\link{ewspec.trend}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' @examples
#' spec <- wavethresh::cns(1024)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 1024))^2)
#'
#' set.seed(2352)
#'
#' noise <- wavethresh::LSWsim(spec)
#' trend <- c(seq(from = 0, to = 4, length = 400), seq(from = 4, to = 0, length = 624))
#'
#' x <- trend + noise
#'
#' spec.est <- ewspec.diff(x, family = "DaubExPhase", filter.number = 1, max.scale = 7)
#'
#' quick.spec.plot(spec.est$S)
#' @export
ewspec.diff <- function(data, lag = 1, filter.number = 1, family = "DaubExPhase",
                        binwidth = floor(2 * sqrt(length(data))), diff.number = 1,
                        max.scale = floor(log2(length(data)) * 0.7), WP.smooth = TRUE,
                        boundary.handle = FALSE, AutoReflect = FALSE,
                        supply.inv.mat = FALSE, inv.mat = NULL) {
  # function that computes the spectral estimate of a time series that has a trend.

  data.check <- ewspec.checks(
    data = data, max.scale = max.scale, lag = lag,
    binwidth = binwidth, boundary.handle = boundary.handle
  )

  data.len <- data.check$data.len
  max.scale <- data.check$max.scale
  boundary.handle <- data.check$boundary.handle
  J <- data.check$J
  dyadic <- data.check$dyadic

  if (boundary.handle == TRUE) {
    data <- get.boundary.timeseries(data, type = "LSW.diff")
  }

  # difference data to remove trend/seasonality and calculate the appropriate correction matrix
  # and its inverse:
  if (diff.number != 1 && diff.number != 2) {
    diff.number <- 1
    warning("Function only implements 1st or 2nd differences. Using 1st difference.")
  }

  if (diff.number == 1) {
    diff.data <- c(diff(data, lag))
    diff.data <- c(diff.data, rep(0, lag))
  } else if (diff.number == 2) {
    diff.data <- c(diff(diff(data)), 0, 0)
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

  data.wd <- locits::ewspec3(diff.data,
    filter.number = filter.number, family = family,
    binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = WP.smooth
  )

  if (boundary.handle == TRUE) {
    data.wd <- smooth.wav.per.calc(
      data.wd = data.wd, J = J, data.len = data.len,
      filter.number = filter.number, family = family,
      dyadic = dyadic, max.scale = max.scale
    )
  }


  l <- S.calc(
    data.wd = data.wd, max.scale = max.scale, J = J, inv.mat = inv.mat,
    filter.number = filter.number, family = family
  )

  return(l)
}
