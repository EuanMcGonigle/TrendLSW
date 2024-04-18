#' @title Print an Object of Class \code{TLSW}
#'
#' @description Prints a \code{TLSW} object, alongside summary information.
#' The first part prints details of the class, specifically the names of elements within.
#' Then prints out the summary, which gives information about a \code{TLSW} object. If spectral
#' estimation was performed, then the type of smoothing and binwidth is printed,
#' along with the differencing performed if it is used, the maximum wavelet
#' scale analysed, and whether or not boundary handling was used. If trend
#' estimation is performed, then the type of wavelet thresholding and transform
#' used is printed, as well as the maximum wavelet scale used, whether or not boundary handling was used,
#' and the significance of the confidence interval if it was calculated.

#'
#' @param x A \code{TLSW} object.
#' @param ... Other arguments.
#'
#' @return No return value, called for side effects
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @seealso \code{\link{TLSW}}, \code{\link{summary.TLSW}}
#' @export
#'
#' @examples
#' # simulates an example time series and estimates its trend and evolutionary wavelet spectrum
#'
#' spec <- wavethresh::cns(512)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2)
#'
#' trend <- seq(from = 0, to = 5, length = 512)
#'
#' set.seed(1)
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' x.TLSW <- TLSW(x)
#'
#' print(x.TLSW)
#'
print.TLSW <- function(x, ...) {
  cat("Class 'TLSW' : Trend Locally Stationary Wavelet Object:\n")
  cat("       ~~~~  : List with", length(x), "components with names\n")
  cat("             ", names(x), "\n\n")
  cat("\nsummary(.):\n-----------\n")
  summary.TLSW(x)
}
