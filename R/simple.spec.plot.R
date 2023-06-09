#' @title Quick Plotting of Evolutionary Wavelet Spectrum
#' @description Plots a standard version of the EWS plot.
#' @param spectrum A spectrum object of class wd which you wish to plot
#' @param bylev If TRUE, plots each level of the spectrum on its own individual
#' scaling.
#' @examples
#' spec <- wavethresh::cns(1024)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 1024))^2)
#'
#' simple.spec.plot(spec)
#' @export
simple.spec.plot <- function(spectrum, bylev = FALSE) {
  scales <- wavethresh::nlevelsWT(spectrum)
  if (bylev == FALSE) {
    plot(spectrum, ylabchars = -(1:scales), main = "", sub = "", ylab = "Scale")
  } else {
    plot(spectrum, ylabchars = -(1:scales), main = "", sub = "", ylab = "Scale", scaling = "by.level")
  }
}
