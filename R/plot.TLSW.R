#' @title Plot Trend and/or Spectrum Information in a \code{TLSW} Object
#'
#' @description Plots information contained within a \code{TLSW} object.
#'
#'
#' @param x A \code{TLSW} object
#' @param plot.type A string object indicating what is to be plotted. Can be:
#' \itemize{
#' \item{\code{"trend"}}{: will plot the trend estimate only.}
#' \item{\code{"spec"}}{: will plot the spectral estimate only.}
#' \item{\code{"both"}}{: will plot both the trend and spectral estimate.}
#' }
#' @param trend.plot.args A list object, that includes any choices for the graphical parameters used for plotting the trend estimate.
#' @param spec.plot.args A list object, that includes any choices for the graphical parameters used for plotting the spectral estimate.
#' @param ... Any additional arguments that will be applied to the graphical parameters of both the trend and spectrum plotting.
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @export
#'
#' @importFrom graphics lines par polygon axis segments title
#' @examples
#' # Simulates an example time series and estimates its trend and evolutionary wavelet spectrum.
#' # Then plots both estimates.
#'
#' spec <- wavethresh::cns(512)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2)
#'
#' trend <- seq(from = 0, to = 5, length = 512)
#'
#' set.seed(1)
#'
#' x <- TLSW.sim(trend = trend, spec = spec)
#'
#' x.TLSW <- TLSW.est(x)
#'
#'plot(x.TLSW, trend.plot.args = list(ylab = "Simulated Data"))
plot.TLSW <- function(x, plot.type = c("both", "trend", "spec")[1],
                      trend.plot.args, spec.plot.args, ...){

  if (x$do.trend.est == FALSE) {
    plot.type <- "spec"
  } else if (x$do.spec.est == FALSE) {
    plot.type <- "trend"
  }

  old.par <- par(no.readonly = TRUE)

  both.args <- list(...)

  if(plot.type == "both") {

    max.plot.scale <- wavethresh::nlevelsWT(x$spec.est$S)

    if(x$trend$calc.confint == FALSE){
      y.min <- min(x$x, x$trend$trend.est)
      y.max <- max(x$x, x$trend$trend.est)
    } else {
      y.min <- min(x$x, x$trend$trend.est, x$trend.est$lower.confint)
      y.max <- max(x$x, x$trend$trend.est, x$trend.est$upper.confint)
    }

    par(mfrow = c(1,2), mar = c(4, 4, 2, 1) + 0.1)

    if(missing(trend.plot.args)){
      trend.plot.args <- list(type = "l", xlab = "Time",
                              ylab = expression(X[t]),
                              ylim = c(y.min, y.max),
                              main = "", col = "grey60")
    } else{
      if(!("type" %in% names(trend.plot.args))){
        trend.plot.args$type <- "l"
      }
      if(!("xlab" %in% names(trend.plot.args))){
        trend.plot.args$xlab <- "Time"
      }
      if(!("ylab" %in% names(trend.plot.args))){
        trend.plot.args$ylab <- expression(X[t])
      }
      if(!("ylim" %in% names(trend.plot.args))){
        trend.plot.args$ylim <- c(y.min, y.max)
      }
      if(!("main" %in% names(trend.plot.args))){
        trend.plot.args$main <- ""
      }
      if(!("col" %in% names(trend.plot.args))){
        trend.plot.args$col <- "grey60"
      }
    }
    if(missing(spec.plot.args)){
      spec.plot.args <- list(ylabchars = (1:max.plot.scale),
                              xlab = "Time", ylab = "Scale", main = "",
                              sub = "")
    } else{
      if(!("ylabchars" %in% names(spec.plot.args))){
        spec.plot.args$ylabchars <- (1:max.plot.scale)
      }
      if(!("xlab" %in% names(spec.plot.args))){
        spec.plot.args$xlab <- "Time"
      }
      if(!("ylab" %in% names(spec.plot.args))){
        spec.plot.args$ylab <- "Scale"
      }
      if(!("sub" %in% names(spec.plot.args))){
        spec.plot.args$sub <- ""
      }
      if(!("main" %in% names(spec.plot.args))){
        spec.plot.args$main <- ""
      }
    }

  } else if(plot.type == "trend"){

    if(x$trend$calc.confint == FALSE){
      y.min <- min(x$x, x$trend$trend.est)
      y.max <- max(x$x, x$trend$trend.est)
    } else {
      y.min <- min(x$x, x$trend$trend.est, x$trend.est$lower.confint)
      y.max <- max(x$x, x$trend$trend.est, x$trend.est$upper.confint)
    }

    par(mfrow = c(1,1), mar = c(4, 4, 2, 1) + 0.1)

    if(missing(trend.plot.args)){
      trend.plot.args <- list(type = "l", xlab = "Time",
                              ylab = expression(X[t]),
                              ylim = c(y.min, y.max),
                              main = "Trend Estimate", col = "grey60")
    } else{
      if(!("type" %in% names(trend.plot.args))){
        trend.plot.args$type <- "l"
      }
      if(!("xlab" %in% names(trend.plot.args))){
        trend.plot.args$xlab <- "Time"
      }
      if(!("ylab" %in% names(trend.plot.args))){
        trend.plot.args$ylab <- expression(X[t])
      }
      if(!("ylim" %in% names(trend.plot.args))){
        trend.plot.args$ylim <- c(y.min, y.max)
      }
      if(!("main" %in% names(trend.plot.args))){
        trend.plot.args$main <- "Trend Estimate"
      }
      if(!("col" %in% names(trend.plot.args))){
        trend.plot.args$col <- "grey60"
      }
    }

  } else if(plot.type == "spec"){

    max.plot.scale <- wavethresh::nlevelsWT(x$spec.est$S)

    par(mfrow = c(1,1),mar = c(4, 4, 2, 1) + 0.1)

    if(missing(spec.plot.args)){
      spec.plot.args <- list(ylabchars = (1:max.plot.scale),
                              xlab = "Time", ylab = "Scale", main = "Spectral Estimate",
                              sub = "")
    } else{
        if(!("ylabchars" %in% names(spec.plot.args))){
          spec.plot.args$ylabchars <- (1:max.plot.scale)
        }
        if(!("xlab" %in% names(spec.plot.args))){
          spec.plot.args$xlab <- "Time"
        }
        if(!("ylab" %in% names(spec.plot.args))){
          spec.plot.args$ylab <- "Scale"
        }
        if(!("sub" %in% names(spec.plot.args))){
          spec.plot.args$sub <- ""
        }
        if(!("main" %in% names(spec.plot.args))){
          spec.plot.args$main <- "Spectral Estimate"
        }
    }
  } else{
    stop("Parameter plot.type not recognised. Can be one of 'both', 'trend,', or 'spec'.")
  }

  if (plot.type == "trend" || plot.type == "both") {
    time.index <- 1:length(x$x)

    trend.plot.args.replace <- both.args[names(both.args) %in% names(trend.plot.args)]

    trend.plot.args[names(trend.plot.args.replace)] <- trend.plot.args.replace

    trend.plot.args.add <- both.args[ !names(both.args) %in% names(trend.plot.args) ]

    trend.plot.args <- as.list(c(trend.plot.args, trend.plot.args.add))

    do.call(plot, c(x$x ~ time.index, trend.plot.args))

    lines(x$trend$trend.est, col = 2, lwd = 2)

    if(x$trend$calc.confint == TRUE){


      polygon(c(time.index,rev(time.index)), c(x$trend.est$lower.confint, rev(x$trend.est$upper.confint)),
              col = "#0000FF33", border = NA)
      lines(x$trend.est$lower.confint, col="blue",lty=2)
      lines(x$trend.est$upper.confint, col="blue",lty=2)

    }


  }

  if (plot.type == "spec" || plot.type == "both") {

    spec.plot.args$n <- length(x$x)
    spec.plot.args$x <- x$spec.est$S

    spec.plot.args.replace <- both.args[names(both.args) %in% names(spec.plot.args)]

    spec.plot.args[names(spec.plot.args.replace)] <- spec.plot.args.replace

    spec.plot.args.add <- both.args[ !names(both.args) %in% names(spec.plot.args) ]

    spec.plot.args <- as.list(c(spec.plot.args, spec.plot.args.add))

    do.call(spec.plot, spec.plot.args)

  }

  if(plot.type == "both"){
    title(substitute(paste(bold("Trend and Spectral Estimate"))),
          line = -1, outer = TRUE)
  }

  suppressWarnings(par(old.par))

}
