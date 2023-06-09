#' @title Simulate locally stationary wavelet processes with random innovations not
#' necessarily Gaussian.
#' @description Simulates a locally stationary wavelet process given a spectrum and
#' distribution for the random innovations. Extension of the LSWsim function
#' from wavethresh.
#' @param spec An object of class wd which contains the spectrum for simulating
#' an LSW process.
#' @param distribution The distribution of the random variables used to
#' simulate the process. Can be "Normal", "Exponential", "Chisquare" or
#' "Poisson".
#' @return A vector simulated from the spectral description given in the spec
#' description. The returned vector will exhibit the spectral characteristics
#' defined by spec.
#' \code{\link{LSWsim}}
#' @examples
#' spec <- wavethresh::cns(1024)
#'
#' spec <- wavethresh::putD(spec, level = 8, seq(from = 2, to = 8, length = 1024))
#'
#' x <- LSWsim.anydist(spec, distribution = "Exponential")
#'
#' plot.ts(x)
#' @export
LSWsim.anydist <- function(spec, distribution = "Normal") {
  # function to simulate LSW processes using innovations from
  # other distribution families

  if (any(spec$D < 0)) {
    stop("All spectral elements must be non-negative")
  }
  nlev <- wavethresh::nlevelsWT(spec)
  len <- 2^nlev
  for (i in (nlev - 1):0) {
    v <- wavethresh::accessD(spec, level = i)
    if (distribution == "Poisson") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rpois(len, 1) - 1)
    } else if (distribution == "Exponential") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rexp(len, rate = 1) - 1)
    } else if (distribution == "Normal") {
      v <- sqrt(v) * 2^(nlev - i) * stats::rnorm(len, mean = 0, sd = 1)
    } else if (distribution == "Chisquare") {
      v <- sqrt(v) * 2^(nlev - i) * (stats::rchisq(n = len, df = 1 / 2) - 1 / 2)
    }
    spec <- wavethresh::putD(spec, level = i, v = v)
  }
  wavethresh::AvBasis(wavethresh::convert(spec))
}
