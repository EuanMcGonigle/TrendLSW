#' @title Prints an object of class \code{TLSW}
#'
#' @description Prints a \code{TLSW} object, giving summary information.
#'
#' @param x A \code{TLSW} object.
#' @param ... Other arguments.
#'
#' @return None
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @seealso \code{\link{TLSW.est}}, \code{\link{summary.TLSW}}
#' @export
#'
#' @examples
print.TLSW <- function(x, ...) {
  cat("Class 'TLSW' : Trend Locally Stationary Wavelet Object:\n")
  cat("       ~~~~  : List with", length(x), "components with names\n")
  cat("             ", names(x), "\n\n")
  cat("\nsummary(.):\n-----------\n")
  summary.TLSW(x)
}
