#' @title Estimation of Evolutionary Wavelet Spectrum for Non-Zero Mean Time Series
#' @description Internal function to compute the evolutionary wavelet spectrum (EWS) estimate from
#' a time series that may include a trend component. The estimate is computed
#' by taking the non-decimated wavelet transform of the time series data,
#' squaring it, smoothing using a running mean, and then correction for bias
#' using the appropriate correction matrix.
#' This function is not intended for general use by regular users of the package.
#' @details Estimates the evolutionary wavelet spectrum of a
#' time series that displays a smooth mean and nonstationary autocovariance.
#' The estimation procedure is as follows:
#'
#' 1. The squared modulus of the non-decimated wavelet transform is computed,
#' known as the raw wavelet periodogram. This is returned by the function.
#'
#' 2. The raw wavelet periodogram is smoothed using a running mean smoother.
#'
#' 3. The smoothed periodogram is bias corrected using the inverse of the bias
#' matrix. The correction is applied across the finest max.scale scales. If the
#' analysing wavelet and generating wavelet are different, this is given by the
#' inverse of the \eqn{C} matrix defined in McGonigle et al. (2022). If they are the
#' same, this is the inverse of the \eqn{A} matrix, defined in Nason et al. (2000).
#' If you are unsure on the filter and wavelet choices, it is recommended to
#' use the same wavelet for generating and analysing purposes.
#'
#' The final estimate, stored in the S component, can be plotted using the plot
#' function, please see the example below.
#' @param x The time series you wish to analyse.
#' @param an.filter.number The index number for the wavelet used to analyse the
#' time series. For the "DaubExPhase" family, the filter number can be between
#' 1 to 10. For the "DaubLeAsymm" family, the filter number can be between 4 to
#' 10. Similarly for gen.filter.number.
#' @param an.family The family of the analysing wavelet. It is recommended to
#' use either the Daubechies Extremal Phase family, or the Daubechies Least
#' Asymmetric family, corresponding to the "DaubExPhase" or the "DaubLeAsymm"
#' options. Similarly for gen.family.
#' @param gen.filter.number The index number for the wavelet that generates the
#' stochastic component of the time series.
#' @param gen.family The family of the generating wavelet.
#' @param binwidth The bin width of the running mean smoother used to smooth
#' the raw wavelet periodogram.
#' @param max.scale The coarsest level to which the time series is analysed to.
#' Should be a positive integer less than J, where T=2^J is the length of the
#' time series. The default setting is 0.7J, to control for bias from the trend
#' and boundary effects.
#' @param S.smooth Argument that dictates if smoothing is performed on the raw
#' wavelet periodogram.
#' @param AutoReflect As in wavethresh. Decides whether or not the time series
#' is reflected when computing the wavelet transform. Helps estimation at the
#' boundaries.
#' @param supply.inv.mat Not intended for general use. If TRUE, user must supply the
#' appropriate correction matrix
#' @param inv.mat If supply.mat is TRUE, user must supply the appropriate
#' correction matrix used to correct the raw wavelet periodogram. Equal to \eqn{C^{-1}}.
#' @param boundary.handle Logical variable, if TRUE, the time series is
#' boundary corrected, to get a more accurate spectrum estimate at the
#' boundaries of the times series. If FALSE, no boundary correction is applied.
#' Recommended to use TRUE.
#' @param smooth.type String indicating which type of smoothing to use on wavelet periodogram.
#' Can be \code{"mean"}, \code{"median"}, or \code{"epan"}. Default is \code{"epan"}.
#' @return A list object, containing the following fields:
#' \item{S}{The evolutionary wavelet spectral estimate of the input data. This object is of
#' class wd and so can be plotted and printed in the usual way using wavethresh
#' functionality. }
#' \item{WavPer}{ The raw wavelet periodogram of the input
#' data. The EWS estimate (above) is the smoothed corrected version of the raw
#' wavelet periodogram.}
#' \item{SmoothWavPer}{ The smoothed, un-corrected raw
#' wavelet periodogram of the input data. }
#' \item{max.scale, boundary.handle, S.smooth, smooth.type, binwidth}{Input parameters}
#' @seealso \code{\link{TLSW}}, \code{\link[wavethresh]{wd}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and
#' adaptive estimation of the evolutionary wavelet spectrum. \emph{Journal of
#' the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{62(2)}, 271--292.
#' @keywords internal
ewspec.trend <- function(x, an.filter.number = 4, an.family = "DaubExPhase",
                         gen.filter.number = an.filter.number, gen.family = an.family,
                         binwidth = floor(2 * sqrt(length(x))),
                         max.scale = floor(log2(length(x)) * 0.7), S.smooth = TRUE,
                         smooth.type = "epan",
                         AutoReflect = TRUE, supply.inv.mat = FALSE, inv.mat = NULL,
                         boundary.handle = TRUE) {
  # function that computes the spectral estimate of a time series that has a smooth trend
  # that can be zeroed out by the wavelet coefficients.

  # user chooses a maximum scale of the wavelet transform to analyse, and
  # binwidth of the running mean smoother.

  x.check <- ewspec.checks(
    x = x, max.scale = max.scale, lag = 1,
    binwidth = binwidth, boundary.handle = boundary.handle,
    S.smooth = S.smooth, smooth.type = smooth.type
  )

  x.len <- x.check$x.len
  max.scale <- x.check$max.scale
  boundary.handle <- x.check$boundary.handle
  J <- x.check$J
  dyadic <- x.check$dyadic

  # calculate the appropriate correction matrix and its inverse:

  if (supply.inv.mat == FALSE) {
    if (an.filter.number == gen.filter.number && an.family == gen.family) {
      A <- wavethresh::ipndacw(J = -max.scale, filter.number = an.filter.number, family = an.family)
      inv.mat <- solve(A)
    } else {
      C <- Cmat.calc(
        J = max.scale, gen.filter.number = gen.filter.number,
        an.filter.number = an.filter.number, gen.family = gen.family,
        an.family = an.family
      )
      inv.mat <- solve(C)
    }
  } else {
    inv.mat <- supply.mat.check(inv.mat = inv.mat, max.scale = max.scale)
  }

  if (boundary.handle == TRUE) {
    x <- get.boundary.timeseries(x, type = "TLSW")
  }
  # calculate raw wavelet periodogram which we need to correct:

  if (smooth.type == "median" || smooth.type == "epan") {
    x.wd <- locits::ewspec3(x,
      filter.number = an.filter.number, family = an.family,
      binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = FALSE
    )
    x.wd <- WP.manual.smooth(
      x.wd = x.wd, smooth.type = smooth.type,
      max.scale = max.scale, binwidth = binwidth
    )
  } else {
    x.wd <- locits::ewspec3(x,
      filter.number = an.filter.number, family = an.family,
      binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = S.smooth
    )
  }

  # access smoothed, uncorrected wavelet periodogram:

  if (boundary.handle == TRUE) {
    x.wd <- smooth.wav.per.calc(
      x.wd = x.wd, J = J, x.len = x.len,
      filter.number = an.filter.number, family = an.family,
      dyadic = dyadic, max.scale = max.scale
    )
  }

  l <- S.calc(
    x.wd = x.wd, max.scale = max.scale, J = J, inv.mat = inv.mat,
    filter.number = gen.filter.number, family = gen.family
  )

  l$max.scale <- max.scale
  l$boundary.handle <- boundary.handle
  l$S.smooth <- S.smooth
  if (S.smooth == TRUE) {
    l$smooth.type <- smooth.type
    l$binwidth <- binwidth
  }
  return(l)
}
