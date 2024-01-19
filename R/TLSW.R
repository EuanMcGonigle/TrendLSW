#' @title Estimate Trend and Spectrum of Trend Locally Stationary Wavelet Process
#' @description Using wavelet-based methods, estimates the trend and evolutionary
#' wavelet spectrum (EWS) of a nonstationary time series. For trend estimation, either
#' a linear or nonlinear wavelet thresholding estimator can be computed. Spectral
#'  estimation can be performed with or without differencing the times series first.
#'
#' @param x The time series you wish to analyse.
#' @param do.spec.est Logical variable, indicating whether spectral estimation is to be performed on the time series.
#' @param do.trend.est Logical variable, indicating whether trend estimation is to be performed on the time series.
#' @param S.filter.number The index number for the wavelet used for spectrum estimation.
#' @param S.family The family of the wavelet used for spectrum estimation.
#' @param S.smooth A logical variable to indicate whether smoothing is performed on the wavelet periodogram.
#' @param S.smooth.type String indicating which type of smoothing to use on wavelet periodogram.
#' Can be one of
#' \itemize{
#' \item{\code{"mean"}: running mean smoother.}
#' \item{ \code{"median"} running median smoother.}
#' \item{\code{"epan"} Epanechnikov kernel smoother.}
#' }
#' @param S.binwidth The bin width of the smoother used to smooth
#' the raw wavelet periodogram.
#' @param S.max.scale The coarsest wavelet scale used to estimate the spectrum.
#' Should be a positive integer less than \eqn{J}, where \eqn{T=2^J} is the length of the
#' time series.
#' @param S.boundary.handle Logical variable, if TRUE, the time series is
#' boundary corrected, to get a more accurate spectrum estimate at the
#' boundaries of the times series. If FALSE, no boundary correction is applied.
#' Recommended to use TRUE.
#' @param S.inv.mat The user can pre-calculate and supply the appropriate
#' correction matrix used to correct the raw wavelet periodogram. If left blank,
#' then the correction matrix is calculated when performing spectral estimation.
#' @param S.do.diff Logical variable, indicating if the time series is to be
#' differenced before spectral estimation is performed.
#' @param S.lag The lag of differencing used, only applicable if \code{S.do.dif = TRUE}.
#' @param S.diff.number The number of differencing operations performed,
#' only applicable if \code{S.do.dif = TRUE}. A first difference is recommended as default.
#' @param T.est.type String indicating type of wavelet thresholding used. Can be "linear", which means
#' that all non-boundary wavelet coefficients are set to zero, or "nonlinear", where
#' each wavelet coefficient is thresholded using a time-varying, noise-dependent threshold.
#' @param T.filter.number The index number for the wavelet used for trend estimation.
#' @param T.family The family of the wavelet used for trend estimation.
#' @param T.transform String giving the type of wavelet transform used for trend estimation.
#' Can be "dec", in which case a standard (decimated) wavelet transform is used, or "nondec",
#' in which case a nondecimated transform is used.
#' @param T.boundary.handle Logical variable, if \code{TRUE}, the time series is
#' boundary corrected when estimating the trend.
#' @param T.max.scale Integer variable, selects the number of scales of the wavelet transform to
#' apply thresholding to for trend estimation.
#' @param T.confint Logical variable. If \code{TRUE}, a \code{(1-T.sig.lvl)} pointwise confidence interval is
#' computed for the trend estimate. For \code{T.est.type = "linear"}, this is
#' computed using the asymptotic distribution of the trend estimator,
#' whilst for \code{T.est.type = "nonlinear"}, it is computed via bootstrapping.
#' @param T.sig.lvl Used only if \code{T.confint = TRUE}; a numeric value
#' (\code{0 <= T.sig.lvl <= 1}) with which a \code{(1-T.sig.lvl)} pointwise
#' confidence interval for the trend estimate is generated.
#' @param T.lacf.max.lag Used only if \code{T.est.type = "linear"} and  \code{calc.confint = TRUE};
#' the maximum lag of the autocovariance to compute needed for calculating the asymptotic confidence interval.
#' @param T.reps Used only if \code{T.est.type = "nonlinear"} and  \code{calc.confint = TRUE}; the number of bootstrap
#' replications used to calculate the confidence interval.
#' @param T.thresh.type String variable, used only if \code{T.est.type = "nonlinear"}; the type of
#' thresholding function used in the trend estimation. Currently only
#' "soft" and "hard" are available.
#' @param T.thresh.normal Logical variable, used only if \code{T.est.type = "nonlinear"};
#' if \code{TRUE}, uses a threshold assuming the data are normally
#' distributed. If \code{FALSE}, uses a larger threshold to reflect non-normality.
#' @param gen.filter.number The index number for the wavelet that generates the
#' stochastic component of the time series. For the "DaubExPhase" family, the filter number can be between
#' 1 to 10. For the "DaubLeAsymm" family, the filter number can be between 4 to 10.
#' Recommended to leave as the default, set to the same as \code{S.filter.number}.
#' @param gen.family The family of the generating wavelet. It is recommended to
#' use either the Daubechies Extremal Phase family, or the Daubechies Least
#' Asymmetric family, corresponding to the "DaubExPhase" or the "DaubLeAsymm"
#' options. Recommended to leave as the default, set to the same as \code{S.family}.
#' @return An object of class \code{"TLSW"}, a list that contains the following components:
#'    \item{x}{Input data}
#'    \item{do.spec.est}{Input parameter, logical variable specifying if spectral estimation was performed.}
#'    \item{spec.est}{A list object, returned if \code{do.spec.est = TRUE}. Contains relevant input parameters
#'    and the following fields related to the spectrum estimate:
#' \itemize{
#'  \item{S}{: The evolutionary wavelet spectral estimate of the input data. This object is of
#' class wd and so can be plotted and printed in the usual way using wavethresh
#' functionality. }
#' \item{WavPer}{: The raw wavelet periodogram of the input
#' data. The EWS estimate (above) is the smoothed corrected version of the
#' wavelet periodogram.}
#' \item{SmoothWavPer}{: The smoothed, un-corrected raw
#' wavelet periodogram of the input data. }
#' }
#' }
#'    \item{do.trend.est}{Input parameter, logical variable specifying if trend estimation was performed.}
#'    \item{trend.est}{A list object, returned if \code{do.trend.est = TRUE}. Contains relevant input parameters
#'    and the following fields related to the trend estimate:
#'    \itemize{
#' \item{trend.est}{: A vector of length \code{length(x)} containing the trend estimate.}
#' \item{lower.confint}{: Returned if \code{calc.confint = TRUE}. The lower limit of the pointwise confidence interval.}
#' \item{upper.confint}{: Returned if \code{calc.confint = TRUE}. The upper limit of the pointwise confidence interval.}
#' \item{T.est.type}{: Type of trend estimator computed, either 'linear' or 'nonlinear'.}
#'    }}
#' @references #' McGonigle, E. T., Killick, R., and Nunes, M. (2022a). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022b). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' @seealso \code{\link{plot.TLSW}}, \code{\link{summary.TLSW}}, \code{\link{print.TLSW}}
#' @examples
#' # simulates an example time series and estimates its trend and evolutionary wavelet spectrum
#'
#' spec <- matrix(0, nrow = 10, ncol = 2^10)
#'
#' spec[1,] = seq(from = 1, to = 10, length = 1024)
#'
#' trend <- sin(pi * (seq(from = 0, to = 4, length = 1024)))
#'
#' set.seed(1)
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' plot.ts(x)
#'
#' x.TLSW <- TLSW(x)
#'
#' summary(x.TLSW)
#'
#' @export
TLSW <- function(x, do.spec.est = TRUE, do.trend.est = TRUE,  S.filter.number = 4,
                     S.family = "DaubExPhase", S.smooth = TRUE,
                     S.smooth.type = c("mean", "median", "epan")[1],
                     S.binwidth = floor(2 * sqrt(length(x))),
                     S.max.scale = floor(log2(length(x)) * 0.7),
                     S.boundary.handle = TRUE, S.inv.mat = NULL,
                     S.do.diff = FALSE, S.lag = 1, S.diff.number = 1,
                     T.est.type = c("linear", "nonlinear")[1], T.filter.number = 4,
                     T.family = "DaubExPhase", T.transform = c("dec", "nondec")[1],
                     T.boundary.handle = TRUE, T.max.scale = floor(log2(length(x)) * 0.7),
                     T.confint = FALSE, T.sig.lvl = 0.05, T.lacf.max.lag = floor(10 * (log10(length(x)))),
                     T.reps = 199, T.thresh.type = c("hard", "soft")[1], T.thresh.normal = TRUE,
                     gen.filter.number = S.filter.number, gen.family = S.family) {
  stopifnot("Both the do.spec.est and do.trend.est parameters have been set to FALSE,
            at least one should be TRUE." = do.spec.est == TRUE || do.trend.est == TRUE)
  if (is.null(S.inv.mat)) {
    supply.inv.mat <- FALSE
  } else {
    supply.inv.mat <- TRUE
  }

  if (do.trend.est == TRUE && do.spec.est == FALSE && (T.est.type == "nonlinear" || T.confint == TRUE)) {
    do.spec.est <- TRUE
    warning("Spectral estimate is needed for trend estimation. Setting do.spec.est = TRUE.")
  }

  if (do.spec.est == TRUE) {
    if (S.do.diff == TRUE) {
      x.spec <- ewspec.diff(
        x = x, lag = S.lag, filter.number = S.filter.number,
        family = S.family, binwidth = S.binwidth, diff.number = S.diff.number,
        max.scale = S.max.scale, S.smooth = S.smooth,
        boundary.handle = S.boundary.handle, AutoReflect = FALSE,
        supply.inv.mat = supply.inv.mat, inv.mat = S.inv.mat
      )
    } else {
      x.spec <- ewspec.trend(
        x = x, gen.filter.number = gen.filter.number,
        gen.family = gen.family, an.filter.number = S.filter.number,
        an.family = S.family, binwidth = S.binwidth,
        max.scale = S.max.scale, S.smooth = S.smooth,
        boundary.handle = S.boundary.handle,
        supply.inv.mat = supply.inv.mat, inv.mat = S.inv.mat,
        smooth.type = S.smooth.type
      )
    }
  }

  if (do.trend.est == TRUE) {
    if (T.est.type == "linear") {
      x.trend <- wav.trend.est(
        x = x, filter.number = T.filter.number,
        family = T.family, max.scale = T.max.scale,
        transform.type = T.transform,
        boundary.handle = T.boundary.handle, calc.confint = T.confint,
        sig.lvl = T.sig.lvl, lag.max = floor(10 * (log10(length(x))))
      )
    } else {
      x.trend <- wav.diff.trend.est(
        x = x, spec.est = x.spec, filter.number = T.filter.number,
        family = T.family, max.scale = T.max.scale,
        transform.type = T.transform,
        thresh.type = T.thresh.type, normal = T.thresh.normal,
        boundary.handle = T.boundary.handle, calc.confint = T.confint,
        reps = T.reps, sig.lvl = T.sig.lvl
      )
    }
    x.trend$T.est.type <- T.est.type
    x.trend <- x.trend[names(x.trend) != "x"]
  }

  if (do.spec.est == TRUE && do.trend.est == TRUE) {
    out <- list(x = x, do.spec.est = do.spec.est, spec.est = x.spec, do.trend.est = do.trend.est, trend.est = x.trend)
  } else if (do.trend.est == FALSE) {
    out <- list(x = x, do.spec.est = do.spec.est, spec.est = x.spec, do.trend.est = do.trend.est)
  } else {
    out <- list(x = x, do.spec.est = do.spec.est, do.trend.est = do.trend.est, trend.est = x.trend)
  }

  class(out) <- "TLSW"
  return(out)
}
