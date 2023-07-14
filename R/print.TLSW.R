#' Title
#'
#' @param x a \code{TLSW} object
#' @param ... not in use
#'
#' @export
#'
#' @examples
print.TLSW <- function(x, ...) {
    cat("Class 'TLSW' : Trend Locally Stationary Wavelet Object:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    summary.TLSW(x)
  }
