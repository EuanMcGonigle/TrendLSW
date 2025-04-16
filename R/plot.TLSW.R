#' @title Plot Trend and/or Spectrum Information in a \code{TLSW} Object
#'
#' @description Plots information contained within a \code{TLSW} object.
#' Depending on the \code{plot.type} option this will produce a plot of the data
#' with trend estimate overlayed, a plot of the spectral estimate, or both (default).
#' If the \code{TLSW} object does not contain trend or spectral estimates and these are requested
#' a warning will be given.
#'
#'
#' @param x A \code{TLSW} object
#' @param plot.type A string object indicating what is to be plotted. Can be:
#' \itemize{
#'    \item \code{"trend"}: will plot the trend estimate only.
#'    \item \code{"spec"}: will plot the spectral estimate only.
#'    \item \code{c("trend", "spec")}: the default value will plot both the trend and spectral estimate.
#' }
#' @param trend.plot.args A list object, that includes any choices for the graphical parameters used for plotting the trend estimate.
#' @param spec.plot.args A list object, that includes any choices for the graphical parameters used for plotting the spectral estimate.
#' @param plot.CI A logical variable. If TRUE, the confidence interval of the trend estimate (if computed) will be included in the plot.
#' @param ... Any additional arguments that will be applied to the graphical parameters of both the trend and spectrum plotting.
#' @details
#' A TLSW object can be plotted using the standard \code{plot} function in R to display the
#' estimated trend function and wavelet spectrum. The estimated trend is visualised using
#' \code{\link[graphics]{plot.default}}. Visualisation of the estimated spectrum is
#' based on \code{\link[wavethresh]{plot.wd}}, for which credit belongs to Guy Nason.
#' Graphical parameters for customising the display of the trend or spectrum plots should be given
#' to the \code{trend.plot.args} and \code{spec.plot.args} arguments respectively.
#' For graphical parameters for the trend plot:
#' \itemize{
#' \item Parameters related to the overall plot should be provided as they usually would be when using the \code{plot} function,
#' in the \code{trend.plot.args} list object. For example, to change the title of the plot to "Plot", use \code{main = "Plot"}.
#' \item Parameters affecting the display of the estimated trend line should begin with the
#' prefix \code{"T."}. For example, to set the colour of the trend line to blue, use
#' \code{T.col = "blue"}.
#' \item Parameters affecting the display of the confidence interval lines should begin with the
#' prefix \code{"CI."}. For example, to set the line width of the confidence interval to 2, use
#' \code{CI.lwd = 2}.
#' \item Parameters affecting the display of the polygon drawn by the confidence interval
#'  should begin with the prefix \code{"poly."}. For example, to set the colour of the
#'  confidence interval region to green, use \code{poly.col = "green"}.
#' }
#' @return No return value, called for side effects
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @seealso \code{\link{TLSW}}, \code{\link{summary.TLSW}}, \code{\link{print.TLSW}}, \code{\link[wavethresh]{plot.wd}}
#' @export
#'
#' @importFrom graphics lines par polygon axis segments title points
#' @examples
#' # Simulates an example time series and estimates its trend and evolutionary wavelet spectrum.
#' # Then plots both estimates.
#'
#' spec <- matrix(0, nrow = 9, ncol = 512)
#'
#' spec[1, ] <- 4 + 4 * sin(seq(from = 0, to = 2 * pi, length = 512))^2
#'
#' trend <- seq(from = 0, to = 10, length = 512) + 2 * sin(seq(from = 0, to = 2 * pi, length = 512))
#'
#' set.seed(1)
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' x.TLSW <- TLSW(x)
#'
#' plot(x.TLSW, trend.plot.args = list(
#'   ylab = "Simulated Data", T.col = 4,
#'   T.lwd = 2, T.lty = 2
#' ))
#'
plot.TLSW <- function(x, plot.type = c("trend", "spec"),
                      trend.plot.args, spec.plot.args, plot.CI,
                      ...) {
  if (any(plot.type == "trend") & (x$do.trend.est == FALSE)) {
    plot.type <- plot.type[-which(plot.type == "trend")]
    warning("trend plot requested but no trend estimated in TLSW object")
  } else if (any(plot.type == "spec") & (x$do.spec.est == FALSE)) {
    plot.type <- plot.type[-which(plot.type == "spec")]
    warning("spec plot requested but no spec estimated in TLSW object")
  }

  if (length(plot.type) == 0) {
    stop("No plot.type specified")
  }

  both.args <- list(...)

  if (any(plot.type == "trend")) {
    if (x$trend$T.CI == FALSE) {
      y.min <- min(x$x, x$trend$trend.est)
      y.max <- max(x$x, x$trend$trend.est)
      if(missing(plot.CI)){
        plot.CI <- FALSE
      }
      if(plot.CI == TRUE){
        warning("Argument plot.CI set to TRUE, but no CI was calculated. Will not display CIs.")
        plot.CI <- FALSE
      }
    } else {
      if(missing(plot.CI)){
        plot.CI <- TRUE
      }
      if(plot.CI == TRUE){
        y.min <- min(x$x, x$trend$trend.est, x$trend.est$lower.CI)
        y.max <- max(x$x, x$trend$trend.est, x$trend.est$upper.CI)
      }else{
        y.min <- min(x$x, x$trend$trend.est)
        y.max <- max(x$x, x$trend$trend.est)
      }
    }

    if (missing(trend.plot.args)) {
      trend.plot.args <- list(
        type = "l", xlab = "Time",
        ylab = expression(X[t]),
        ylim = c(y.min, y.max),
        main = "Trend Estimate", col = "grey60"
      )
      T.args <- list(col = 2, lwd = 2)
      poly.args <- list(col = "#0000FF40", border = NA)
      CI.args <- list(col = "blue", lty = 2)
    } else {
      T.names <- names(trend.plot.args)[grep("T.", names(trend.plot.args))]
      poly.names <- names(trend.plot.args)[grep("poly.", names(trend.plot.args))]
      CI.names <- names(trend.plot.args)[grep("CI.", names(trend.plot.args))]

      name.func <- function(y, num) {
        substr(y, num, nchar(y))
      }

      if (length(T.names) > 0) {
        T.args <- trend.plot.args[T.names]
        names(T.args) <- sapply(names(T.args), name.func, num = 3)
      } else {
        T.args <- list(col = 2, lwd = 2)
      }
      if (length(poly.names) > 0) {
        poly.args <- trend.plot.args[poly.names]
        names(poly.args) <- sapply(names(poly.args), name.func, num = 6)
      } else {
        poly.args <- list(col = "#0000FF40", border = NA)
      }
      if (length(CI.names) > 0) {
        CI.args <- trend.plot.args[CI.names]
        names(CI.args) <- sapply(names(CI.args), name.func, num = 4)
      } else {
        CI.args <- list(col = "blue", lty = 2)
      }

      all.names <- c(T.names, poly.names, CI.names)

      trend.plot.args <- trend.plot.args[!names(trend.plot.args) %in% all.names]

      if (!("type" %in% names(trend.plot.args))) {
        trend.plot.args$type <- "l"
      }
      if (!("xlab" %in% names(trend.plot.args))) {
        trend.plot.args$xlab <- "Time"
      }
      if (!("ylab" %in% names(trend.plot.args))) {
        trend.plot.args$ylab <- expression(X[t])
      }
      if (!("ylim" %in% names(trend.plot.args))) {
        trend.plot.args$ylim <- c(y.min, y.max)
      }
      if (!("main" %in% names(trend.plot.args))) {
        trend.plot.args$main <- "Trend Estimate"
      }
      if (!("col" %in% names(trend.plot.args))) {
        trend.plot.args$col <- "grey60"
      }
    }

    time.index <- 1:length(x$x)

    trend.plot.args.replace <- both.args[names(both.args) %in% names(trend.plot.args)]

    trend.plot.args[names(trend.plot.args.replace)] <- trend.plot.args.replace

    trend.plot.args.add <- both.args[!names(both.args) %in% names(trend.plot.args)]

    trend.plot.args <- as.list(c(trend.plot.args, trend.plot.args.add))

    if (x$trend$T.CI == TRUE) {
      trend.point.type <- trend.plot.args$type
      trend.plot.args$type <- "n"
      do.call(plot, c(x$x ~ time.index, trend.plot.args))

      trend.plot.args$type <- trend.point.type
      trend.plot.args$xlab <- NULL
      trend.plot.args$ylab <- NULL
      trend.plot.args$main <- NULL

      if (plot.CI == TRUE) {
        do.call(polygon, c(
          list(x = c(time.index, rev(time.index)), y = c(x$trend.est$lower.CI, rev(x$trend.est$upper.CI))),
          poly.args
        ))
        do.call(points, c(x$x ~ time.index, trend.plot.args))
        do.call(lines, c(x$trend.est$lower.CI ~ time.index, CI.args))
        do.call(lines, c(x$trend.est$upper.CI ~ time.index, CI.args))
      }else{
        do.call(points, c(x$x ~ time.index, trend.plot.args))
      }
    } else {
      do.call(plot, c(x$x ~ time.index, trend.plot.args))
    }
    do.call(lines, c(x$trend.est$T ~ time.index, T.args))
  }

  if (any(plot.type == "spec")) {
    max.plot.scale <- wavethresh::nlevelsWT(x$spec.est$S)

    if (missing(spec.plot.args)) {
      spec.plot.args <- list(
        ylabchars = (1:max.plot.scale),
        xlab = "Time", ylab = "Scale", main = "Spectral Estimate",
        sub = ""
      )
    } else {
      if (!("ylabchars" %in% names(spec.plot.args))) {
        spec.plot.args$ylabchars <- (1:max.plot.scale)
      }
      if (!("xlab" %in% names(spec.plot.args))) {
        spec.plot.args$xlab <- "Time"
      }
      if (!("ylab" %in% names(spec.plot.args))) {
        spec.plot.args$ylab <- "Scale"
      }
      if (!("sub" %in% names(spec.plot.args))) {
        spec.plot.args$sub <- ""
      }
      if (!("main" %in% names(spec.plot.args))) {
        spec.plot.args$main <- "Spectral Estimate"
      }
    }
    spec.plot.args$n <- length(x$x)
    spec.plot.args$x <- x$spec.est$S

    spec.plot.args.replace <- both.args[names(both.args) %in% names(spec.plot.args)]

    spec.plot.args[names(spec.plot.args.replace)] <- spec.plot.args.replace

    spec.plot.args.add <- both.args[!names(both.args) %in% names(spec.plot.args)]

    spec.plot.args <- as.list(c(spec.plot.args, spec.plot.args.add))

    do.call(spec.plot, spec.plot.args)
  }
}
