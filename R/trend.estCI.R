#' @title Calculate Confidence Intervals Of Wavelet-Based Trend Estimate
#' @description Calculates appropriate confidence intervals for the trend
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
#' @export
trend.estCI <- function(trend.est, lacf.est, filter.number = 4, family = "DaubLeAsymm", alpha = 0.95) {
  # function to create confidence interval for the trend estimate

  data.len <- length(trend.est)

  size <- 1 - alpha
  qval <- stats::qnorm(1 - size / 2)

  cov.mat <- create.covmat(lacf.est, data.len)

  # calculate the wavelet transform matrix:

  W <- t(wavethresh::GenW(n = data.len, filter.number = filter.number, family = family))

  scale <- log2(data.len)

  boundary_test <- c(rep(0, data.len - 1), 1)

  y_wd <- wavethresh::wd(boundary_test, family = family, filter.number = filter.number)

  boundary_coefs <- list(1:scale)

  for (i in (scale - 1):1) {
    temp <- wavethresh::accessD(y_wd, level = i)

    boundary_coefs[[i]] <- (which(temp != 0))
  }

  boundary_vec <- c(1, rep(0, data.len - 2), 1)

  for (i in (scale - 1):3) {
    boundary_vec[(data.len - 2^(i + 1) + 2):(data.len - 2^(i + 1) + 1 + 2^i)][boundary_coefs[[i]]] <- 1
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
