#' @title Compute Localised Autocovariance Estimate of a TLSW Object
#' @description Computes the local autocovariance and autocorrelation estimates, given an
#' input of an object of class TLSW containing the estimated spectrum. Provides the same functionality as the
#' function \code{lacf} from the \code{locits} package, but user provides an object of
#' class \code{TLSW} as the main argument.
#' @param x.TLSW a \code{TLSW} object.
#' @param lag.max The maximum lag of acf required. If NULL then the same
#' default as in the regular acf function is used.
#' @return An object of class \code{lacf} which contains the following components:
#'  \itemize{
#'  \item \code{lacf}: a matrix containing the estimate of the local autocovariance. Columns represent lags
#'  (beginning at lag 0), and rows represent time points.
#' \item \code{lacr}: a matrix containing the estimate of the local autocorrelation. Columns represent lags
#' (beginning at lag 0), and rows represent time points.
#' \item \code{name}: the name of the time series (if applicable).
#' \item \code{date}: the date the function was executed.
#' \item \code{SmoothWP}: The smoothed, un-corrected raw
#' wavelet periodogram of the input data.
#' \item \code{S}: the spectral estimate used to compute the local autocovariance.
#' \item \code{J}: the number of total wavelet scales.
#' }
#' @seealso \code{\link[locits]{lacf}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#'
#' Nason, G. P. (2013). A test for second-order stationarity and approximate
#' confidence intervals for localized autocovariances for locally stationary
#' time series. \emph{Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology)}, \bold{75(5)}, 879--904.
#'
#' Nason, G. P. (2016). locits: Tests of stationarity and localized autocovariance.
#' R package version 1.7.3.
#' @examples
#'
#' ## ---- computes estimate of local autocovariance function
#'
#'
#' spec <- matrix(0, nrow = 9, ncol = 512)
#' spec[2, ] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2
#'
#' trend <- seq(from = 0, to = 10, length = 512)
#'
#' set.seed(123)
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' ## ---- first estimate the spectrum:
#'
#' x.TLSW <- TLSW(x)
#'
#' #---- estimate the lacf:
#'
#' lacf.est <- TLSWlacf(x.TLSW)
#'
#' #---- plot the variance (lag 0 lacf) over time:
#'
#' plot.ts(lacf.est$lacf[, 1], ylab = "Variance")
#' @export
TLSWlacf <- function(x.TLSW, lag.max = NULL) {
  stopifnot("Argument lag.max should be a nonegative integer." = lag.max >= 0)
  stopifnot("Argument x.TLSW should be an object of class TLSW." = isa(x.TLSW, "TLSW"))

  x <- x.TLSW$x
  dsname <- deparse(substitute(x))

  S <- x.TLSW$spec.est$S
  SmoothWP <- x.TLSW$spec.est$SmoothWavPer

  J <- S$nlevels
  Smat <- matrix(S$D, nrow = 2^J, ncol = J)[1:length(x), ]
  Psi <- wavethresh::PsiJmat(-J,
    filter.number = S$filter$filter.number,
    family = S$filter$family
  )
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
