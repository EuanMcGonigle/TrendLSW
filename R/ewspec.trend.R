#' @title Estimation of Evolutionary Wavelet Spectrum for Non-Zero Mean Time Series
#' @description Computes the evolutionary wavelet spectrum (EWS) estimate from
#' a time series that may include a trend component. The estimate is computed
#' by taking the non-decimated wavelet transform of the time series data,
#' squaring it, smoothing using a running mean, and then correction for bias
#' using the appropriate correction matrix.
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
#' @param data The time series you wish to analyse.
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
#' @param WP.smooth Argument that dictates if smoothing is performed on the raw
#' wavelet periodogram.
#' @param AutoReflect As in wavethresh. Decides whether or not the time series
#' is reflected when computing the wavelet transform. Helps estimation at the
#' boundaries.
#' @param supply.mat Not intended for general use. If TRUE, user must supply the
#' appropriate correction matrix
#' @param mat If supply.mat is TRUE, user must supply the appropriate
#' correction matrix used to correct the raw wavelet periodogram. Equal to \eqn{C^{-1}}.
#' @param boundary.handle Logical variable, if TRUE, the time series is
#' boundary corrected, to get a more accurate spectrum estimate at the
#' boundaries of the times series. If FALSE, no boundary correction is applied.
#' Recommended to use TRUE.
#' @return A list object, containing the following fields:
#' \itemize{\item{S}{The evolutionary wavelet spectral estimate of the input data. This object is of
#' class wd and so can be plotted and printed in the usual way using wavethresh
#' functionality. }
#' \item{WavPer}{ The raw wavelet periodogram of the input
#' data. The EWS estimate (above) is the smoothed corrected version of the raw
#' wavelet periodogram. }
#' \item{SmoothWavPer}{ The smoothed, un-corrected raw
#' wavelet periodogram of the input data. }
#' }
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and
#' adaptive estimation of the evolutionary wavelet spectrum. \emph{Journal of
#' the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{62(2)}, 271--292.
#' @examples
#' # simulates an example time series and estimates its evolutionary wavelet spectrum
#'
#' spec <- wavethresh::cns(512)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2)
#'
#' noise <- wavethresh::LSWsim(spec)
#' trend <- seq(from = 0, to = 5, length = 512)
#'
#' x <- trend + noise
#'
#' spec.est <- ewspec.trend(x,
#'   an.filter.number = 4, an.family = "DaubExPhase",
#'   gen.filter.number = 1, gen.family = "DaubExPhase"
#' )
#'
#' quick.spec.plot(spec.est$S)
#' @export
ewspec.trend <- function(data, an.filter.number = 10, an.family = "DaubLeAsymm",
                         gen.filter.number = an.filter.number, gen.family = an.family,
                         binwidth = floor(2 * sqrt(length(data))),
                         max.scale = floor(log2(length(data)) * 0.7), WP.smooth = TRUE,
                         AutoReflect = TRUE, supply.mat = FALSE, mat = NULL,
                         boundary.handle = TRUE) {
  # function that computes the spectral estimate of a time series that has a smooth trend
  # that can be zeroed out by the wavelet coefficients.

  # user chooses a maximum scale of the wavelet transform to analyse, and
  # binwidth of the running mean smoother.

  data.check <- ewspec.checks(data = data, max.scale = max.scale, lag = 1,
                                     binwidth = binwidth, boundary.handle = boundary.handle)

  data.len <- data.check$data.len
  max.scale <- data.check$max.scale
  boundary.handle <- data.check$boundary.handle
  J <- data.check$J
  dyadic <- data.check$dyadic

  # calculate the appropriate correction matrix and its inverse:

  if (supply.mat == FALSE) {
    if (an.filter.number == gen.filter.number && an.family == gen.family) {
      A <- wavethresh::ipndacw(J = -max.scale, filter.number = an.filter.number, family = an.family)
      inv <- solve(A)
    } else {
      C <- Cmat.calc(
        J = max.scale, gen.filter.number = gen.filter.number,
        an.filter.number = an.filter.number, gen.family = gen.family,
        an.family = an.family
      )
      inv <- solve(C)
    }
  } else {
    stopifnot("Supplied inverse matrix must be square" = nrow(mat) == ncol(mat))
    stopifnot("Dimension of supplied inverse matrix must be larger than max.scale"
              = nrow(mat) >= max.scale)
    inv <- mat[1:max.scale,1:max.scale]
  }

  if (boundary.handle == TRUE) {
    data <- get.boundary.timeseries(data, type = "TLSW")
  }
  # calculate raw wavelet periodogram which we need to correct:

  data.wd <- locits::ewspec3(data,
    filter.number = an.filter.number, family = an.family,
    binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = WP.smooth
  )

  calc.final.spec <- function(spec, dyadic, data.len) {
    if (dyadic == TRUE) {
      final_spec <- wavethresh::cns(2^(spec$nlevels - 2),
        filter.number = spec$filter$filter.number,
        family = spec$filter$family
      )

      lower <- 2^(spec$nlevels - 2) + 2^(spec$nlevels - 3) + 1
      upper <- 2^(spec$nlevels - 1) + 2^(spec$nlevels - 3)


      for (j in 0:(spec$nlevels - 3)) {
        bh_d <- wavethresh::accessD(spec, level = j + 2)[lower:upper]

        final_spec <- wavethresh::putD(final_spec, level = j, bh_d)
      }

      return(final_spec)
    } else {
      est.spec.J <- spec$nlevels

      final.spec.J <- floor(log2(data.len)) + 1

      final_spec <- wavethresh::cns(2^final.spec.J,
        filter.number = spec$filter$filter.number,
        family = spec$filter$family
      )

      lower <- floor((2^est.spec.J - data.len) / 2)
      upper <- lower + data.len - 1


      for (j in 0:(final.spec.J - 1)) {
        bh_d <- c(wavethresh::accessD(spec, level = j + (est.spec.J - final.spec.J))[lower:upper], rep(0, 2^final.spec.J + lower - upper - 1))

        final_spec <- wavethresh::putD(final_spec, level = j, bh_d)
      }

      return(final_spec)
    }
  }


  # access smoothed, uncorrected wavelet periodogram:

  if (boundary.handle == TRUE) {
    temp <- locits::ewspec3(rep(0, 2^J), filter.number = an.filter.number, family = an.family)

    temp$SmoothWavPer <- calc.final.spec(data.wd$SmoothWavPer, dyadic = dyadic, data.len = data.len)
    temp$WavPer <- calc.final.spec(data.wd$WavPer, dyadic = dyadic, data.len = data.len)

    if (max.scale < J) {
      for (j in 0:(J - 1 - max.scale)) {
        temp$WavPer <- wavethresh::putD(temp$WavPer, level = j, rep(0, 2^J))
        temp$SmoothWavPer <- wavethresh::putD(temp$SmoothWavPer, level = j, rep(0, 2^J))
      }
    }


    data.wd <- temp
  }

  uncor.spec <- data.wd$SmoothWavPer

  # perform correction: mutiply by inverse matrix, non-estimated scales are set to zero.

  uncor.spec.mat <- matrix(0, nrow = max.scale, ncol = 2^J)

  for (j in 1:max.scale) {
    uncor.spec.mat[j, ] <- wavethresh::accessD(uncor.spec, level = J - j)
  }

  # perform correction step:

  cor.spec.mat <- inv %*% uncor.spec.mat

  # now fill in wd object with final spectrum estimate.

  S <- wavethresh::cns(2^J, filter.number = gen.filter.number, family = gen.family)

  for (j in 1:max.scale) {
    S <- wavethresh::putD(S, level = J - j, cor.spec.mat[j, ])
  }

  # return final estimate, along with smoothed and unsmoothed periodogram:

  l <- list(S = S, WavPer = data.wd$WavPer, SmoothWavPer = uncor.spec)

  return(l)
}
