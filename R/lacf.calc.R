#' @title Compute Localised Autocovariance Estimate from Spectrum Estimate
#' @description Computes the local autocovariance and autocorrelation estimates, given an
#' input of a spectrum estimate. Provides the same functionality as the
#' function \code{lacf} from the \code{locits} package, but user provides the spectrum
#' estimate in the argument.
#' @param x The time series you wish to analyse.
#' @param filter.number Wavelet filter number that generated the time series.
#' @param family Wavelet family that generated the time series.
#' @param spec.est Estimated spectrum from which the lacf estimate will be
#' calculated, the component \code{spec.est} in the output of the \code{TLSW} functions.
#' @param lag.max The maximum lag of acf required. If NULL then the same
#' default as in the regular acf function is used.
#' @param ... Further arguments to be passed to perform spectrum estimation, only to be used if \code{spec.est} is not supplied.
#' @return An object of class \code{lacf} which contains the autocovariance.
#' @seealso \code{\link{lacf}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' Nason, G. P. (2013). A test for second-order stationarity and approximate
#' confidence intervals for localized autocovariances for locally stationary
#' time series. \emph{Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology)}, \bold{75(5)}, 879--904.
#' @examples
#'
#' ## ---- computes estimate of local autocovariance function
#'
#' ## ---- example where LSW process is generated using the Haar wavelet
#'
#' spec <- wavethresh::cns(512)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2)
#'
#' noise <- wavethresh::LSWsim(spec)
#' trend <- seq(from = 0, to = 10, length = 512)
#'
#' x <- trend + noise
#'
#' ## ---- first estimate the spectrum using Daubechies EP4 wavelet:
#'
#' x.TLSW <- TLSW(x)
#'
#' spec.est <- x.TLSW$spec.est
#'
#' #---- estimate the lacf:
#'
#' lacf.est <- lacf.calc(x = x, spec.est = spec.est)
#'
#' plot.ts(lacf.est$lacf[, 1])
#' @export
lacf.calc <- function(x, filter.number = 4, family = "DaubExPhase",
                      spec.est = NULL, lag.max = NULL, ...) {
  stopifnot("Paramter lag.max should be a nonegative integer." = lag.max >= 0)

  if (is.null(spec.est)) {
    spec.est <- ewspec.trend(
      x = x, an.filter.number = filter.number, an.family = family,
      gen.filter.number = filter.number, gen.family = family, ...
    )
  }

  dsname <- deparse(substitute(x))

  S <- spec.est$S
  SmoothWP <- spec.est$SmoothWavPer

  J <- S$nlevels
  Smat <- matrix(S$D, nrow = 2^J, ncol = J)[1:length(x), ]
  Psi <- wavethresh::PsiJmat(-J, filter.number = filter.number, family = family)
  nc <- ncol(Psi)
  L <- (nc - 1) / 2
  dimnames(Psi) <- list(NULL, c(-L:0, 1:L))
  if (is.null(lag.max)) {
    lag.max <- floor(10 * (log10(length(x))))
  }
  if (L + 1 + lag.max > ncol(Psi)) {
    warning(paste(
      "lag.max too high. Have reset it to ",
      ncol(Psi) - L - 1, ". Higher lags are zero"
    ))
    lag.max <- ncol(Psi) - L - 1
  }
  the.lacf <- Smat %*% Psi[, (L + 1):(L + 1 + lag.max)]
  the.lacor <- sweep(the.lacf, 1, the.lacf[, 1], FUN = "/")
  l <- list(
    lacf = the.lacf, lacr = the.lacor, name = dsname,
    date = date(), SmoothWP = SmoothWP, S = S, J = J
  )
  class(l) <- "lacf"
  return(l)
}
