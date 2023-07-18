#' @title Plot Trend and/or Spectrum Information in a \code{TLSW} Object
#'
#' @description Plots information contained within a \code{TLSWß} object.
#'
#'
#' @param x A \code{TLSW} object
#' @param plot.type A string object indicating what is to be plotted. Can be "trend", in which case
#' the trend estimate (and associated confidence intervals if calculated) are plotted, or "spec",
#' in which case the spectral estimate is plotted, or "both", in which case both estimates are plotted.
#' @param ... Additional plotting arguments used for the spectrum plotting.
#'
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @export
#'
#' @importFrom graphics lines par polygon
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
#' x <- TLSW.sim(trend = trend, spec = spec)
#'
#' x.TLSW <- TLSW.est(x)
#'
#'plot(x.TLSW)
plot.TLSW <- function(x, plot.type = c("both", "trend", "spec")[1], ...) {
  if (x$do.trend.est == FALSE) {
    plot.type <- "spec"
  } else if (x$do.spec.est == FALSE) {
    plot.type <- "trend"
  }

  if( plot.type == "both") {
    par(mfrow = c(1,2))
  }

  if (plot.type == "trend" || plot.type == "both") {

    if(x$trend$calc.confint == FALSE){
      y.min <- min(x$x, x$trend$trend.est)
      y.max <- max(x$x, x$trend$trend.est)
      plot(x$x, type = "l", xlab = "Time", ylab = expression(X[t]), ylim = c(y.min,y.max))
      lines(x$trend$trend.est, col = 2, lwd = 2)

    } else {
      y.min <- min(x$x, x$trend$trend.est, x$trend.est$lower.confint)
      y.max <- max(x$x, x$trend$trend.est, x$trend.est$upper.confint)
      time.index <- 1:length(x$x)
      plot(x$x, type = "l", xlab = "Time", ylab = expression(X[t]), ylim = c(y.min,y.max))
      polygon(c(time.index,rev(time.index)), c(x$trend.est$lower.confint, rev(x$trend.est$upper.confint)),
      col = 'grey85', border = NA)
      lines(x$trend.est$lower.confint, col="blue",lty=2)
      lines(x$trend.est$upper.confint, col="blue",lty=2)
      lines(x$trend$trend.est, col = 2, lwd = 2)
    }


  }

  if (plot.type == "spec" || plot.type == "both") {
    max.plot.scale <- wavethresh::nlevelsWT(x$spec.est$S)

    wavethresh::plot.wd(x$spec.est$S, ylabchars = (1:max.plot.scale),
                        xlab = "Time", ylab = "Scale", main = "", sub = "", ...)
  }





}
