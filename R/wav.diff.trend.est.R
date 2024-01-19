#' @title Wavelet Thresholding Trend Estimation of Time Series
#' @description Computes the wavelet thresholding trend estimate for a time
#' series that may be second-order nonstationary. The function calculates the
#' wavelet transform of the time series, thresholds the coefficients based on
#' an estimate of their variance, and inverts to give the trend estimate.
#' @details Estimates the trend function of a locally stationary time series, by
#' incorporating the evolutionary wavelet spectrum estimate in a wavelet
#' thresholding procedure. To use this function, first compute the spectral
#' estimate of the time series, using the function ewspec.diff.
#'
#' The function works as follows:
#'
#' 1. The wavelet transform of the time series is calculated.
#'
#' 2. The wavelet coefficients are individually thresholded using the universal
#' threshold \eqn{\hat{\sigma}\sqrt(2 log T)}, where \eqn{\hat{sigma}^2} is an estimate of their variance. The variance
#' estimate is calculated using the spectral estimate, supplied by the user in
#' the \code{spec} argument.
#'
#' 3. The inverse wavelet transform is applied to obtain the final estimate.
#'
#' @param x The time series you want to estimate the trend function of.
#' @param spec.est You must supply the estimate of the evolutionary wavelet
#' spectrum of the time series. This is the output of the \code{ewspec.diff}
#' function.
#' @param filter.number Selects the index of the wavelet used in the estimation
#' procedure. For Daubechies compactly supported wavelets the filter number is
#' the number of vanishing moments.
#' @param thresh.type The type of thresholding function used. Currently only
#' "soft" and "hard" are available.
#' @param normal If TRUE, uses a threshold assuming the data are normally
#' distributed. If FALSE, uses a larger threshold to reflect non-normality.
#' @param family Selects the wavelet family to use. Recommended to only use the
#' Daubechies compactly supported wavelets DaubExPhase and DaubLeAsymm.
#' @param transform.type String giving the type of wavelet transform used.
#' Can be "dec", in which case a standard (decimated) wavelet transform is used, or "nondec",
#' in which case a nondecimated transform is used.
#' @param max.scale Selects the number of scales of the wavelet transform to
#' apply thresholding to. Should be a value from 1 (finest) to J-1 (coarsest),
#' where T=2^J is the length of the time series. Recommended to use 2J/3
#' scales.
#' @param boundary.handle Logical variable, decides if boundary handling should
#' be applied to the time series before estimation.
#' @param calc.confint Logical variable. If \code{TRUE}, a bootstrapped \code{(1-sig.lvl)}
#' pointwise confidence interval is computed for the trend estimate.
#' @param sig.lvl Used only if \code{calc.confint = TRUE}; a numeric value
#' (\code{0 <= sig.lvl <= 1}) with which a \code{(1-sig.lvl)} pointwise
#' confidence interval for the trend estimate is generated.
#' @param reps Used only if \code{calc.confint = TRUE}; the number of bootstrap
#' replications used to calcualte the confidence interval.
#' @param ...  Further arguments to be passed to the \code{\link{ewspec.diff}}
#' call, only to be used if \code{calc.confint = TRUE}.
#' @return A \code{list} object containing the following fields:
#' \item{x}{Input data}
#' \item{filter.number, family}{Input wavelet parameters}
#' \item{transform.type, max.scale, boundary.handle, thresh.type, normal,  calc.confint}{Input parameters}
#' \item{T}{A vector of length \code{length(x)} containing the trend estimate}
#' \item{lower.confint}{Returned if \code{calc.confint = TRUE}. The lower limit of the pointwise confidence interval}
#' \item{upper.confint}{Returned if \code{calc.confint = TRUE}. The upper limit of the pointwise confidence interval}
#' \item{reps}{Returned if \code{calc.confint = TRUE}. The number of bootstrap replicates used to compute
#'  pointwise confidence interval}
#' \item{sig.lvl}{Returned if \code{calc.confint = TRUE}. The significance level of the pointwise confidence interval}
#' @seealso \code{\link{ewspec.diff}}, \code{\link{wav.trend.est}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' @examples
#' spec <- wavethresh::cns(1024, filter.number = 4)
#' spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 1024))^2)
#'
#' set.seed(120)
#'
#' sine.trend <- -2 * sin(seq(from = 0, to = 2 * pi, length = 1024)) -
#'   1.5 * cos(seq(from = 0, to = pi, length = 1024))
#'
#' x <- TLSWsim(trend = sine.trend, spec = spec)
#'
#' spec.est <- ewspec.diff(x = x, family = "DaubExPhase", filter.number = 4, max.scale = 7)
#'
#' trend.est <- wav.diff.trend.est(x = x, spec = spec.est)
#'
#' plot.ts(x, lty = 1, col = 8)
#' lines(sine.trend, col = 2, lwd = 2)
#' lines(trend.est$T, col = 4, lwd = 2, lty = 2)
#' @keywords internal
#' @noRd
wav.diff.trend.est <- function(x, spec.est, filter.number = 4, family = "DaubExPhase",
                               thresh.type = c("hard","soft")[1], normal = TRUE,
                               transform.type = c("dec", "nondec")[1],
                               max.scale = floor(0.7 * log2(length(x))),
                               boundary.handle = FALSE, calc.confint = FALSE,
                               reps = 199, sig.lvl = 0.05, ...) {
  x.check <- trend.est.checks(
    x = x, max.scale = max.scale, boundary.handle = boundary.handle,
    transform.type = transform.type, calc.confint = calc.confint,
    reps = reps, sig.lvl = sig.lvl, est.type = "nonlinear"
  )

  x.len <- x.check$x.len
  max.scale <- x.check$max.scale
  boundary.handle <- x.check$boundary.handle
  J <- x.check$J
  dyadic <- x.check$dyadic

  spec <- spec.est$S

  orig.x <- x
  if (boundary.handle == TRUE) {
    x <- get.boundary.timeseries(x)
  }

  x.len <- length(x)
  J <- wavethresh::IsPowerOfTwo(x.len)

  if(transform.type == "nondec"){
    x.wd <- wavethresh::wd(x, filter.number = filter.number, family = family, type = "station")
  }else{
    x.wd <- wavethresh::wd(x, filter.number = filter.number, family = family)
  }

  # calculate C to use in variance estimate of wavelet coefficients

  C <- Cmat.calc(
    J = max.scale, an.filter.number = filter.number,
    gen.filter.number = spec$filter$filter.number,
    an.family = family, gen.family = spec$filter$family
  )

  spec.mat <- matrix(0, nrow = max.scale, ncol = length(wavethresh::accessD(spec, level = 0)))

  for (j in 1:max.scale) {
    spec.mat[j, ] <- wavethresh::accessD(spec, level = wavethresh::nlevelsWT(spec) - j)
  }

  # calculate variance estimate matrix

  var.mat <- C %*% spec.mat

  var.mat <- replace.neg.values(var.mat, max.scale)

  x.wd.thresh <- x.wd

  if (boundary.handle == TRUE) {
    # below code determines the boundary coefficients for a given wavelet

    boundary.test <- c(rep(0, x.len - 1), 1)

    y.wd <- wavethresh::wd(boundary.test, family = family, filter.number = filter.number, type = "station")

    # create EWS estimate matrix

    bc.var.mat <- matrix(0, nrow = max.scale, ncol = x.len)

    lower <- floor((x.len - length(orig.x)) / 2)
    upper <- lower + length(orig.x) - 1

    bc.var.mat[, lower:upper] <- var.mat[, 1:length(orig.x)]


    for (j in 1:max.scale) {
      dj <- wavethresh::accessD(x.wd, level = J - j)

      if (normal == TRUE) {
        thresh <- sqrt(2 * bc.var.mat[j, ] * log(x.len))
      } else {
        thresh <- sqrt(bc.var.mat[j, ]) * log(x.len)
      }
      if(transform.type == "dec") {
        thresh <- thresh[(1:(x.len/(2^j)))*2^j-2^j+1]
      }

      temp1 <- dj

      temp1[abs(dj) < thresh] <- 0

      if (thresh.type == "soft") {
        temp2 <- dj[abs(dj) >= thresh]
        temp1[abs(dj) >= thresh] <- sign(temp2) * (abs(temp2) - thresh[abs(dj) >= thresh])
      }

      x.wd.thresh <- wavethresh::putD(x.wd.thresh, level = J - j, temp1)
    }

    # invert the thresholded object to get trend estimate:

    if (transform.type == "dec") {
      trend.est <- wavethresh::wr(x.wd.thresh)
    } else if (transform.type == "nondec") {
      trend.est <- wavethresh::AvBasis(wavethresh::convert(x.wd.thresh))
    }


    if (dyadic == TRUE) {
      lower <- 2^(J - 2) + 2^(J - 3) + 1
      upper <- 2^(J - 1) + 2^(J - 3)
    } else {
      lower <- floor((x.len - length(orig.x)) / 2)
      upper <- lower + length(orig.x) - 1
    }
    trend.est <- trend.est[lower:upper]

  } else {
    # threshold the wavelet coefficients using the variance estimate and user
    # inputted rules

    for (j in 1:max.scale) {
      dj <- wavethresh::accessD(x.wd, level = J - j)

      if (normal == TRUE) {
        thresh <- sqrt(2 * var.mat[j, ] * log(x.len))
      } else {
        thresh <- sqrt(var.mat[j, ]) * log(x.len)
      }
      if(transform.type == "dec") {
        thresh <- thresh[(1:(x.len/(2^j)))*2^j-2^j+1]
      }

      temp1 <- dj

      temp1[abs(dj) < thresh] <- 0

      if (thresh.type == "soft") {
        temp2 <- dj[abs(dj) >= thresh]
        temp1[abs(dj) >= thresh] <- sign(temp2) * (abs(temp2) - thresh[abs(dj) >= thresh])
      }

      x.wd.thresh <- wavethresh::putD(x.wd.thresh, level = J - j, temp1)
    }

    # invert the thresholded object to get trend estimate:

    if (transform.type == "dec") {
      trend.est <- wavethresh::wr(x.wd.thresh)
    } else if (transform.type == "nondec") {
      trend.est <- wavethresh::AvBasis(wavethresh::convert(x.wd.thresh))
    }
  }

  if (calc.confint == TRUE) {
    trend.confint <- trend.estCI.diff(
      x = orig.x, trend.est = trend.est, spec.est = spec.est,
      filter.number = filter.number, thresh.type = thresh.type,
      normal = normal, transform.type = transform.type,
      boundary.handle = boundary.handle, family = family, max.scale = max.scale,
      reps = reps, sig.lvl = sig.lvl, ...
    )
    return(list(
      x = orig.x, T = trend.est, lower.confint = trend.confint[1, ],
      upper.confint = trend.confint[2, ], sig.lvl = sig.lvl, reps = reps,
      filter.number = filter.number, family = family,
      transform.type = transform.type, max.scale = max.scale,
      boundary.handle = boundary.handle, thresh.type = thresh.type,
      normal = normal, calc.confint = calc.confint

    ))
  } else {
    return(list(
      x = orig.x, T = trend.est,  filter.number = filter.number, family = family,
      transform.type = transform.type, max.scale = max.scale,
      boundary.handle = boundary.handle, thresh.type = thresh.type, normal = normal,
      calc.confint = calc.confint
    ))
  }
}
