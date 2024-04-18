#' @title Summary of Output Provided by the \code{TLSW} Function
#' @description Summary method for objects of class \code{TLSW}.
#'
#' @details Prints out information about a \code{TLSW} object. If spectral
#' estimation was performed, then the type of smoothing and binwidth is printed,
#' along with the differencing performed if it is used, the maximum wavelet
#' scale analysed, and whether or not boundary handling was used. If trend
#' estimation is performed, then the type of wavelet thresholding and transform
#' used is printed, as well as the maximum wavelet scale used, whether or not boundary handling was used,
#' and the significance of the confidence interval if it was calculated.
#'
#' @param object A \code{TLSW} object.
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
#' @seealso \code{\link{TLSW}}, \code{\link{print.TLSW}}
#' @export
#'
#' @examples
#' # simulates an example time series and estimates its trend and evolutionary wavelet spectrum
#'
#' spec <- matrix(0, nrow = 10, ncol = 2^10)
#'
#' spec[1, ] <- seq(from = 1, to = 10, length = 1024)
#'
#' trend <- sin(pi * (seq(from = 0, to = 4, length = 1024)))
#'
#' set.seed(1)
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' x.TLSW <- TLSW(x)
#'
#' summary(x.TLSW)
#'
summary.TLSW <- function(object, ...) {
  if (object$do.spec.est == TRUE) {
    cat("Spectral estimation was performed:\n")
    if (object$spec.est$S.smooth == TRUE) {
      cat("-smoothing was performed using ", object$spec.est$smooth.type, " smoothing with binwidth ", object$spec.est$binwidth, ".\n", sep = "")
    } else {
      cat("-no smoothing was performed.\n")
    }
    if (!is.null(object$spec.est$lag)) {
      if (object$spec.est$diff.number == 2) {
        cat("-time series was second differenced before wavelet transform applied.\n")
      } else {
        cat("-time series was first differenced at lag ", object$spec.est$lag, " before wavelet transform applied.\n", sep = "")
      }
    }
    cat("-maximum wavelet scale analysed is scale ", object$spec.est$max.scale, ".\n", sep = "")
    if (object$spec.est$boundary.handle == TRUE) {
      cat("-boundary handling was used.\n")
    } else {
      cat("-no boundary handling was perfomed.\n")
    }
  } else {
    cat("Spectral estimation was not performed.\n")
  }
  cat("-----------\n")

  if (object$do.trend.est == TRUE) {
    cat("Trend estimation was performed:\n")
    if (object$trend.est$T.est.type == "linear") {
      cat("-estimation was performed using a ", object$trend.est$transform.type, "imated wavelet transform with ", object$trend.est$T.est.type, " thresholding.\n", sep = "")
    } else {
      cat("-estimation was performed using a ", object$trend.est$transform.type, "imated wavelet transform with ", object$trend.est$T.est.type, " thresholding using a ", object$trend.est$thresh.type, " threshold.\n", sep = "")
    }

    cat("-maximum wavelet scale analysed is scale ", object$trend.est$max.scale, ".\n", sep = "")
    if (object$trend.est$boundary.handle == TRUE) {
      cat("-boundary handling was used.\n")
    } else {
      cat("-no boundary handling was perfomed.\n")
    }
    if (object$trend.est$T.CI == TRUE) {
      cat("-a pointwise ", 100 * (1 - object$trend.est$sig.lvl), "% confidence interval was calculated.", sep = "")
    }
  } else {
    cat("Trend estimation was not performed.\n")
  }
}
