#' @title Simulate Locally Stationary Wavelet Process with Specified
#' Distribution of Random Innovations
#' @description Simulates a locally stationary wavelet process given a spectrum and
#' distribution for the random innovations. Extension of the \code{LSWsim} function
#' from the \code{wavethresh} package.
#' @param spec An object of class wd which contains the spectrum for simulating
#' an LSW process.
#' @param distribution The distribution of the random variables used to
#' simulate the process. Can be "norm", "pois", "exp", "chisq" or
#' "t".
#' @param rate The rate parameter, only used if \code{distribution = "exp"}
#' or \code{distribution = "pois"}
#' @param df The degrees of freedom, only used if \code{distribution = "chisq"}
#' or \code{distribution = "t"}
#' @return A vector simulated from the spectral description given in the spec
#' description. The returned vector will exhibit the spectral characteristics
#' defined by spec.
#' \code{\link{LSWsim}}
#' @examples
#' spec <- wavethresh::cns(1024)
#'
#' spec <- wavethresh::putD(spec, level = 8, seq(from = 2, to = 8, length = 1024))
#'
#' x <- LSWsim.anydist(spec, distribution = "exp")
#'
#' plot.ts(x)
#' @export
LSWsim.anydist <- function(spec, distribution = c("norm", "pois", "exp", "chisq", "t")[1], rate = NULL, df = NULL) {
  if (any(spec$D < 0)) {
    stop("All spectral elements must be non-negative.")
  }

  stopifnot("Error: distribution must be one of 'norm', 'pois',
               'exp, 'chisq', or 't'." = distribution == "norm" || distribution
  == "pois" || distribution == "exp" || distribution == "chisq" ||
    distribution == "t")

  stopifnot("parameter df must be positive." = df > 0)
  stopifnot("The rate parameter must be positive." = rate > 0)

  if (is.null(df)) {
    if (distribution == "chisq") {
      df <- 1 / 2
    } else if (distribution == "t") {
      df <- 5
    }
  }
  if (is.null(rate)) {
    if (distribution == "exp" || distribution == "pois") {
      rate <- 1
    }
  }

  stopifnot("for t distribution, parameter df must be greater than 2." = distribution != "t" || df > 2)


  nlev <- wavethresh::nlevelsWT(spec)
  len <- 2^nlev
  for (i in (nlev - 1):0) {
    v <- wavethresh::accessD(spec, level = i)
    if (distribution == "pois") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rpois(len, rate) - rate)/(sqrt(rate))
    } else if (distribution == "exp") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rexp(len, rate = rate) - 1/rate)*rate
    } else if (distribution == "norm") {
      v <- sqrt(v) * 2^(nlev - i) * stats::rnorm(len, mean = 0, sd = 1)
    } else if (distribution == "chisq") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rchisq(n = len, df = df) - df) / (sqrt(2 * df))
    } else if (distribution == "t") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rt(n = len, df = df) / sqrt(df / (df - 2)))
    }
    spec <- wavethresh::putD(spec, level = i, v = v)
  }
  wavethresh::AvBasis(wavethresh::convert(spec))
}
