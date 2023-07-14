#' @title Summary of output provided by the \code{TLSW.est} function
#' @description Summary method for objects of class \code{TLSW}
#'
#'
#' @param object a \code{TLSW} object
#' @param ... not in use
#'
#' @export
#'
#' @examples
summary.TLSW <- function(object, ...) {

  if(object$do.spec.est == TRUE){
    cat("Spectral estimation was performed:\n")
    if(object$spec$WP.smooth == TRUE){
      cat("-smoothing was performed using ", object$spec$smooth.type, " smoothing with binwidth ", object$spec$binwidth, ".\n", sep = "")
    } else{
      cat("-no smoothing was performed.\n")
    }
    if(!is.null(object$spec$lag)){
      if(object$spec$diff.number == 2){
        cat("-time series was second differenced before wavelet transform applied.\n")
      } else{
        cat("-time series was first differenced at lag ", object$spec$lag," before wavelet transform applied.\n", sep = "")
      }

    }
    cat("-maximum wavelet scale analysed is scale ", object$spec$max.scale, ".\n", sep = "")
    if(object$spec$boundary.handle == TRUE){
      cat("-boundary handling was used.\n")
    } else{
      cat("-no boundary handling was perfomed.\n")
    }

  } else {
    cat("Spectral estimation was not performed.\n")
  }
  cat("----------------\n")

  if(object$do.trend.est == TRUE){
    cat("Trend estimation was performed:\n")
    if(object$trend$T.est.type == "linear"){
      cat("-estimation was performed using a ", object$trend$T.transform, "imated wavelet transform with ", object$trend$T.est.type," thresholding.\n", sep = "")
    } else {
      cat("-estimation was performed using a ", object$trend$T.transform, "imated wavelet transform with ", object$trend$T.est.type," thresholding using a ", object$trend$T.thresh.type, " threshold.\n", sep = "")
    }

    cat("-maximum wavelet scale analysed is scale ", object$trend$max.scale, ".\n", sep = "")
    if(object$trend$boundary.handle == TRUE){
      cat("-boundary handling was used.\n")
    } else{
      cat("-no boundary handling was perfomed.\n")
    }

  } else {
    cat("Trend estimation was not performed.\n")
  }


}
