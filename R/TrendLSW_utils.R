#' @title Compute Covariance Matrix
#' @description Internal function for computing covariance matrix
#' @keywords internal
#' @noRd


create.covmat <- function(lacf, data.len) {
  lacf <- lacf$lacf

  max.lag <- length(lacf[1, ])

  cov.mat <- matrix(0, nrow = data.len, ncol = data.len)

  for (row in 1:data.len) {
    for (column in row:(min((max.lag - 1 + row), data.len))) {
      cov.mat[row, column] <- lacf[floor((row + column) / 2), (abs(column - row) + 1)]
      cov.mat[column, row] <- cov.mat[row, column]
    }
  }

  cov.mat
}


#' @title Calculate Confidence Intervals Of Wavelet-Based Trend Estimate
#' @description Internal function to alculate appropriate confidence intervals for the trend
#' estimate, based on a given trend estimate and covariance matrix of the time
#' series.
#' @param trend.est The trend estimate of the time series.
#' @param lacf.est The estimated local autocovariance of the time series.
#' @param filter.number The index of the wavelet used in the estimation
#' procedure.
#' @param family The wavelet family used in the estimation procedure.
#' @param alpha The significance level of the confidence interval to calculate.
#' Default is 0.95.
#' @return A list object, containing the following objects:
#' \item{trend.est}{The trend estimate of the input data. } \item{trend.var}{
#' The variance estimate of the trend estimate. } \item{lower.conf}{ The lower
#' pointwise confidence interval for the trend estimate. } \item{upper.conf}{
#' The upper pointwise confidence interval for the trend estimate. }
#' @seealso \code{\link{wav.trend.est}}, \code{\link{lacf.calc}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, , 43(6), 895-917.
#' @keywords internal
#' @noRd
trend.estCI <- function(trend.est, lacf.est, filter.number = 4, family = "DaubLeAsymm", alpha = 0.95,
                        max.scale = floor(log2(length(trend.est)) * 0.7)) {
  # function to create confidence interval for the trend estimate

  data.len <- length(trend.est)

  size <- 1 - alpha
  qval <- stats::qnorm(1 - size / 2)

  J <- log2(data.len)

  cov.mat <- create.covmat(lacf.est, data.len)

  # calculate the wavelet transform matrix:

  W <- t(wavethresh::GenW(n = data.len, filter.number = filter.number, family = family))

  boundary_test <- c(rep(0, data.len - 1), 1)

  y_wd <- wavethresh::wd(boundary_test, family = family, filter.number = filter.number)

  boundary_coefs <- list()

  for (j in (J - 1):(J - max.scale)) {
    temp <- wavethresh::accessD(y_wd, level = j)

    boundary_coefs[[j]] <- (which(temp != 0))
  }

  boundary_vec <- rep(1, data.len)

  for (j in (J - 1):(J - max.scale)) {
    boundary_vec[(data.len - 2^(j + 1) + 2):(data.len - 2^(j + 1) + 1 + 2^j)][-boundary_coefs[[j]]] <- 0
  }


  # calculate the linear thresholding operation matrix:

  A <- diag(boundary_vec)

  # calculate the covariance matrix of the trend estimate:

  R <- t(W) %*% A %*% W

  trend.cov <- (R) %*% cov.mat %*% t(R)

  trend.sd <- suppressWarnings(sqrt(diag(trend.cov)))

  v.nas <- which(is.na(trend.sd))

  trend.sd.final <- trend.sd

  if (length(v.nas > 0)) {
    for (i in 1:length(v.nas)) {
      trend.sd.final[v.nas[i]] <- trend.sd[which(!is.nan(trend.sd))][which.min(abs(v.nas[i] - which(!is.nan(trend.sd))))]
    }
  }

  l <- list(trend.est = trend.est, trend.var = trend.sd.final^2, lower.conf = trend.est - trend.sd.final * qval, upper.conf = trend.est + trend.sd.final * qval)

  return(l)
}




#' @title Spectral Estimation Error Checks
#' @description Internal function for error checking for spectral estimation
#' @keywords internal
#' @noRd
ewspec.checks <- function(data, max.scale, binwidth, lag, boundary.handle) {
  if (any(is.na(data))) {
    stop("Data contains mising values.")
  }
  if (!is.numeric(data)) {
    stop("Data is not numeric")
  }

  data.len <- length(data)

  if (max.scale %% 1 != 0) {
    stop("max.scale parameter must be an integer.")
  }
  if (max.scale < 1 || max.scale > floor(log2(data.len))) {
    warning("max.scale parameter is outside valid range. Resetting to default value.")
    max.scale <- floor(log2(data.len) * 0.7)
  }
  if (binwidth %% 1 != 0) {
    stop("binwidth parameter must be an integer.")
  }

  J <- wavethresh::IsPowerOfTwo(data.len)

  if (is.na(J) == TRUE) {
    warning("Data length is not power of two. Boundary correction has been applied.")
    boundary.handle <- TRUE
    dyadic <- FALSE
    J <- floor(log2(data.len)) + 1
  } else {
    dyadic <- TRUE
  }

  if (!is.numeric(lag)) {
    stop("The lag parameter should be a positive integer.")
  }
  if ((length(lag) != 1) || (lag %% 1 != 0) || (lag <= 0)) {
    stop("The lag parameter should be a positive integer.")
  }

  return(list(
    data.len = data.len, max.scale = max.scale, boundary.handle = boundary.handle,
    J = J, dyadic = dyadic
  ))
}


#' @title Calculate Boundary Handled Spectrum
#' @description Internal function for spectral estimation used when there is
#' boundary handling
#' @keywords internal
#' @noRd

calc.final.spec <- function(spec, dyadic, data.len) {
  if (dyadic == TRUE) {
    final_spec <- wavethresh::cns(2^(spec$nlevels - 2),
      filter.number = spec$filter$filter.number,
      family = spec$filter$family
    )

    lower <- 2^(spec$nlevels - 2) + 2^(spec$nlevels - 3) + 1
    upper <- 2^(spec$nlevels - 1) + 2^(spec$nlevels - 3)


    for (j in 0:(spec$nlevels - 3)) {
      bh_d <- wavethresh::accessD(spec, level = j + 2)[lower:upper]

      final_spec <- wavethresh::putD(final_spec, level = j, bh_d)
    }

    return(final_spec)
  } else {
    est.spec.J <- spec$nlevels

    final.spec.J <- floor(log2(data.len)) + 1

    final_spec <- wavethresh::cns(2^final.spec.J,
      filter.number = spec$filter$filter.number,
      family = spec$filter$family
    )

    lower <- floor((2^est.spec.J - data.len) / 2)
    upper <- lower + data.len - 1


    for (j in 0:(final.spec.J - 1)) {
      bh_d <- c(wavethresh::accessD(spec, level = j + (est.spec.J - final.spec.J))[lower:upper], rep(0, 2^final.spec.J + lower - upper - 1))

      final_spec <- wavethresh::putD(final_spec, level = j, bh_d)
    }

    return(final_spec)
  }
}

#' @title Calculate Spectrum Estimate
#' @description Internal function for calculating spectral estimate
#' @keywords internal
#' @noRd

S.calc <- function(data.wd, max.scale, J, inv.mat, filter.number, family) {
  # access smoothed,uncorrected wavelet periodogram:

  uncor.spec <- data.wd$SmoothWavPer

  # perform correction: mutiply by inverse matrix, non-estimated scales are set to zero.

  uncor.spec.mat <- matrix(0, nrow = max.scale, ncol = 2^J)

  for (j in 1:max.scale) {
    uncor.spec.mat[j, ] <- wavethresh::accessD(uncor.spec, level = J - j)
  }

  # perform correction step:

  cor.spec.mat <- inv.mat %*% uncor.spec.mat

  # now fill in wd object with final spectrum estimate.

  S <- wavethresh::cns(2^J, filter.number = filter.number, family = family)

  for (j in 1:max.scale) {
    S <- wavethresh::putD(S, level = J - j, cor.spec.mat[j, ])
  }

  # return final EWS estimate, along with smoothed and unsmoothed periodogram:

  l <- list(S = S, WavPer = data.wd$WavPer, SmoothWavPer = uncor.spec)

  return(l)
}

#' @title Calculate oundary Handled Smoothed Periodogram Estimate
#' @description Internal function for calculating smoothed spectral estimate when
#' boundary handling is used
#' @keywords internal
#' @noRd
smooth.wav.per.calc <- function(data.wd, J, data.len, filter.number, family,
                                dyadic, max.scale) {
  temp <- locits::ewspec3(rep(0, 2^J), filter.number = filter.number, family = family)

  temp$SmoothWavPer <- calc.final.spec(data.wd$SmoothWavPer, dyadic = dyadic, data.len = data.len)
  temp$WavPer <- calc.final.spec(data.wd$WavPer, dyadic = dyadic, data.len = data.len)

  if (max.scale < J) {
    for (j in 0:(J - 1 - max.scale)) {
      temp$WavPer <- wavethresh::putD(temp$WavPer, level = j, rep(0, 2^J))
      temp$SmoothWavPer <- wavethresh::putD(temp$SmoothWavPer, level = j, rep(0, 2^J))
    }
  }

  return(temp)
}

#' @title Matrix Error Checks
#' @description Internal function for check user-supplied matrix
#' @keywords internal
#' @noRd
supply.mat.check <- function(inv.mat, max.scale) {
  stopifnot("Supplied inverse matrix must be square" = nrow(inv.mat) == ncol(inv.mat))
  stopifnot(
    "Dimension of supplied inverse matrix must be larger than max.scale" = nrow(inv.mat) >= max.scale
  )
  inv.mat <- inv.mat[1:max.scale, 1:max.scale]

  return(inv.mat)
}


#' @title Replace Negative Values in Variance Estimate
#' @description Internal function to replace negative values in the variance
#' estimate used in the \code{wav.diff.trend.est} function
#' @keywords internal
#' @noRd
replace.neg.values <- function(var.mat, max.scale) {
  for (j in 1:max.scale) {
    var.row <- var.mat[j, ]

    if (sum(var.row <= 0) > 0) {
      var.row[var.row < 0] <- 0
      var.row0 <- which(var.row == 0)
      var.row.non0 <- var.row[which(var.row != 0)]

      for (i in 1:length(var.row0)) {
        var.mat[j, (var.row0[i])] <- var.row.non0[which.min(abs(var.row0[i] - which(var.row != 0)))]
      }
    }
  }

  var.mat
}


