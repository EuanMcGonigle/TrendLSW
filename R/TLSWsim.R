#' @title Simulate Trend Locally Stationary Wavelet Process
#' @description Simulates a trend locally stationary wavelet process with a given trend function
#' and spectrum. Extends the \code{LSWsim} function from the \code{wavethresh} package.
#'
#' @param trend Either:
#' \itemize{
#' \item{A numeric vector of length \eqn{n} giving the values of the deterministic
#' trend function,}
#' \item{A real-valued function of one argument defined on rescaled time \eqn{[0,1)}.}
#' }
#' When using a numeric vector for \code{trend}, if \eqn{n} is not a power of 2 then \code{spec} must be specified using
#' a numeric matrix of dimensions \eqn{\lfloor \log_2 (n) \rfloor \times n}.
#' @param spec Either:
#' \itemize{
#' \item{A \code{wavethresh} object of class wd which contains the spectrum for simulating
#' an LSW process,}
#' \item{A numeric matrix of dimensions \eqn{J \times n}, where the \eqn{j}-th row corresponds to the spectrum values at scale \eqn{j} and
#' \eqn{\lfloor \log_2 (n) \rfloor = J},}
#' \item{A list of length \eqn{J=\log_2(n)}, where the \eqn{j}-th element of the list is a function of one argument specifying the spectrum
#' function at scale \eqn{j} on rescaled time \eqn{[0,1)}.}
#' }
#' When using a numeric matrix for \code{spec}, if \eqn{n} is not a power of 2 then \code{trend}
#' must be specified using a numeric vector of length \eqn{n}.
#' @param innov.func A function with first argument \code{n} used for simulating the innovations. By default,
#' normal random innovations are sampled using the \code{rnorm} function.
#' @param filter.number The filter number for the wavelet used to simulate the LSW process (default 4)
#' @param family The family of the wavelet used to simulate the LSW process (default \code{DaubExPhase}).
#' @param ... Optional arguments to be passed to the function  \code{innov.func} for
#' sampling the innovation process.
#' @return A \eqn{n}-length vector containing a TLSW process simulated from the trend and spectral description given by the trend and
#' spec arguments.
#' @seealso \code{\link[wavethresh]{LSWsim}}
#' @examples
#'
#' #---- simulate with numeric trend, and spec a wd object as in wavethresh-----
#'
#' spec <- wavethresh::cns(1024)
#'
#' spec <- wavethresh::putD(spec, level = 8, seq(from = 2, to = 8, length = 1024))
#'
#' trend <- sin(pi * (seq(from = 0, to = 4, length = 1024)))
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' plot.ts(x)
#'
#' #---- simulate with numeric trend, and spec a matrix, with non-dyadic n-----
#'
#' spec <- matrix(0, nrow = 9, ncol = 1000)
#'
#' spec[1, ] <- seq(from = 1, to = 10, length = 1000)
#'
#' trend <- sin(pi * (seq(from = 0, to = 4, length = 1000)))
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' plot.ts(x)
#'
#' #---- simulate with functional trend, and spec a list of functions-----
#'
#' spec <- vector(mode = "list", length = 10)
#'
#' spec[[1]] <- function(u) {
#'   1 + 9 * u
#' }
#'
#' trend <- function(u) {
#'   sin(pi * u)
#' }
#'
#' x <- TLSWsim(trend = trend, spec = spec)
#'
#' plot.ts(x)
#'
#' @export
TLSWsim <- function(trend, spec, filter.number = 4, family = "DaubExPhase",
                    innov.func, ...) {
  if (missing(trend)) {
    trend <- function(u) {
      0
    }
  }

  if (missing(innov.func)) {
    innov.func <- stats::rnorm
  }

  if (!is.function(innov.func)) {
    stop("Argument'innov.func' should be a function.")
  }
  if (names(formals(innov.func))[1] != "n") {
    stop("Invalid 'innov.func' argument: should be in the rnorm family of functions.")
  }

  # error check:

  stopifnot("Argument trend must be either a numeric vector or function." = is.function(trend) || is.numeric(trend))
  stopifnot("Argument spec must be either a numeric matrix, list of functions, or wd object." = isa(spec, "list") || isa(spec, "wd") || is.matrix(spec))

  if (isa(spec, "list")) {
    J <- length(spec)
    n <- 2^J

    stopifnot("The length of the argument spec should be at least 2." = J >= 2)

    for (j in 1:J) {
      stopifnot("The elements of the list spec should be a function, or the NULL value." = is.function(spec[[j]]) || is.null(spec[[j]]))
    }

    spec.mat <- matrix(0, nrow = J, ncol = n)
    vals <- (0:(n - 1)) / n

    for (j in 1:J) {
      if (is.null(spec[[j]])) {
        spec.mat[j, ] <- 0
      } else {
        spec.mat[j, ] <- Vectorize(spec[[j]])(vals)
      }
    }

    spec <- mat.to.spec(spec.mat,
      filter.number = filter.number,
      family = family
    )
  } else if (isa(spec, "wd")) {
    J <- spec$nlevels
    n <- 2^J
    family <- spec$filter$family
    filter.number <- spec$filter$filter.number
  } else if (is.matrix(spec)) {
    J <- nrow(spec)
    n <- ncol(spec)

    if (floor(log2(n)) != J) {
      stop("Dimensions of spec matrix incorrect. The integer part of log2(number of columns) should equal the number of rows.")
    }

    if(is.na(wavethresh::IsPowerOfTwo(n))){
      if (!is.numeric(trend)) {
        stop("If spec has a non-dyadic number of columns, then trend must be a numeric vector.")
      }
      extended.spec <- matrix(0, nrow = J + 1, ncol = 2^(J+1))
      for(j in 1:J){
        extended.spec[j, ] <- c(spec[j, ], rep(spec[j, n], 2^(J+1) - n))
      }
      spec <- extended.spec
    }

    spec <- mat.to.spec(spec,
      filter.number = filter.number,
      family = family
    )
  }

  if (is.numeric(trend)) {
    stopifnot("Length of trend does not match dimensions of spec." = n == length(trend))
  }


  if (any(spec$D < 0)) {
    stop("All spectral elements must be non-negative.")
  }

  if (is.function(trend)) {
    vals <- (0:(n - 1)) / n
    T <- Vectorize(trend)(vals)
  } else {
    T <- trend
  }


  if(is.na(wavethresh::IsPowerOfTwo(n))){
    final.n <- 2^(ceiling(log2(n)))
    J <- J + 1
  }else{
    final.n <- n
  }

  for (i in (J - 1):0) {
    v <- wavethresh::accessD(spec, level = i)
    v <- sqrt(v) * 2^(J - i) * innov.func(final.n, ...)
    spec <- wavethresh::putD(spec, level = i, v = v)
  }

  x <- T + wavethresh::AvBasis(wavethresh::convert(spec))[1:n]
  return(x)
}
