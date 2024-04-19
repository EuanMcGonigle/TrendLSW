#' @title Compute Covariance Matrix
#' @description Internal function for computing covariance matrix
#' @keywords internal
#' @noRd


create.covmat <- function(lacf, x.len) {
  lacf <- lacf$lacf

  max.lag <- length(lacf[1, ])

  cov.mat <- matrix(0, nrow = x.len, ncol = x.len)

  for (row in 1:x.len) {
    for (column in row:(min((max.lag - 1 + row), x.len))) {
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
#' The variance estimate of the trend estimate. } \item{lower.CI}{ The lower
#' pointwise confidence interval for the trend estimate. } \item{upper.CI}{
#' The upper pointwise confidence interval for the trend estimate. }
#' @seealso \code{\link{wav.trend.est}}, \code{\link{TLSWlacf}}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, , 43(6), 895-917.
#' @keywords internal
#' @noRd
#'
#'
#'
#'
trend.estCI <- function(trend.est, lacf.est, filter.number = 4, family = "DaubLeAsymm", sig.lvl = 0.05,
                        max.scale = floor(log2(length(trend.est)) * 0.7)) {
  # function to create confidence interval for the trend estimate

  x.len <- length(trend.est)

  qval <- stats::qnorm(1 - sig.lvl / 2)

  J <- log2(x.len)

  cov.mat <- create.covmat(lacf.est, x.len)

  # calculate the wavelet transform matrix:

  W <- t(wavethresh::GenW(n = x.len, filter.number = filter.number, family = family))

  boundary_test <- c(rep(0, x.len - 1), 1)

  y_wd <- wavethresh::wd(boundary_test, family = family, filter.number = filter.number)

  boundary_coefs <- list()

  for (j in (J - 1):(J - max.scale)) {
    temp <- wavethresh::accessD(y_wd, level = j)

    boundary_coefs[[j]] <- (which(temp != 0))
  }

  boundary_vec <- rep(1, x.len)

  for (j in (J - 1):(J - max.scale)) {
    boundary_vec[(x.len - 2^(j + 1) + 2):(x.len - 2^(j + 1) + 1 + 2^j)][-boundary_coefs[[j]]] <- 0
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

  l <- list(trend.est = trend.est, trend.var = trend.sd.final^2, lower.CI = trend.est - trend.sd.final * qval, upper.CI = trend.est + trend.sd.final * qval)

  return(l)
}


#' @title Calculate Bootstrapped Confidence Intervals Of Wavelet-Based Trend Estimate
#' @description Internal function to calculate appropriate confidence intervals
#' for the nonlinear wavelet thresholding trend estimate, based on a given
#' trend estimate and spectral estimate of the time series.
#' @keywords internal
#' @noRd
trend.est.CI.bootstrap <- function(x, trend.est, spec.est, filter.number = 4, thresh.type = "soft",
                                   normal = TRUE, transform.type = c("dec", "nondec")[2],
                                   family = "DaubLeAsymm", max.scale = floor(log2(length(x)) * 0.7),
                                   boundary.handle = TRUE, reps = 199, sig.lvl = 0.05,
                                   confint.type = c("percentile", "normal")[1],
                                   diff = TRUE, ...) {
  trend.mat <- matrix(0, nrow = reps, ncol = length(x))

  spec <- spec.est$S

  spec$D[spec$D < 0] <- 0

  A <- wavethresh::ipndacw(
    J = -spec.est$max.scale, filter.number = spec$filter$filter.number,
    family = spec$filter$family
  )

  if (diff == TRUE) {
    if (is.null(spec.est$lag)) {
      spec.est$lag <- 1
    }
    if (is.null(spec.est$diff.number)) {
      spec.est$diff.number <- 1
    }
    A1 <- Atau.mat.calc(
      J = spec.est$max.scale, filter.number = spec$filter$filter.number,
      family = spec$filter$family, lag = spec.est$lag
    )
    inv.mat <- solve(2 * A - 2 * A1)
  }


  for (i in 1:reps) {
    rep.x <- trend.est + wavethresh::LSWsim(spec)[1:length(x)]

    if (diff == TRUE) {
      rep.spec <- suppressWarnings(ewspec.diff(rep.x,
        lag = spec.est$lag,
        filter.number = spec$filter$filter.number,
        family = spec$filter$family, binwidth = spec.est$binwidth,
        max.scale = spec.est$max.scale, boundary.handle = spec.est$boundary.handle,
        supply.inv.mat = TRUE, inv.mat = inv.mat,
        diff.number = spec.est$diff.number, ...
      ))

      rep.trend <- suppressWarnings(wav.diff.trend.est(
        x = rep.x, spec.est = rep.spec, filter.number = filter.number,
        family = family, max.scale = max.scale, transform.type = transform.type,
        boundary.handle = boundary.handle,
        thresh.type = thresh.type, normal = normal, T.CI = FALSE
      ))
    } else {
      rep.trend <- suppressWarnings(wav.trend.est(
        x = rep.x, filter.number = filter.number,
        family = family, max.scale = max.scale, transform.type = transform.type,
        boundary.handle = boundary.handle, T.CI = FALSE
      ))
    }

    trend.mat[i, ] <- rep.trend$T
  }

  if (confint.type == "percentile") {
    conf.int <- apply(trend.mat, 2, FUN = stats::quantile, probs = c(sig.lvl / 2, (1 - sig.lvl / 2)))
  } else {
    sd.est <- stats::qnorm(1 - sig.lvl / 2) * apply(trend.mat, 2, FUN = stats::sd)
    conf.int <- rbind(trend.est - sd.est, trend.est + sd.est)
  }
  return(conf.int)
}



#' @title Spectral Estimation Error Checks
#' @description Internal function for error checking for spectral estimation
#' @keywords internal
#' @noRd
ewspec.checks <- function(x, max.scale, binwidth, lag, boundary.handle, S.smooth, smooth.type) {
  if (any(is.na(x))) {
    stop("Data contains mising values.")
  }
  if (!is.numeric(x)) {
    stop("Data is not numeric")
  }
  stopifnot("Parameter S.boundary.handle must be logical variable" = is.logical(boundary.handle))
  stopifnot("Parameter S.smooth must be logical variable" = is.logical(S.smooth))
  stopifnot("Smoothing type must be one of 'mean', 'median', or 'epan'." = smooth.type == "mean" ||
    smooth.type == "median" || smooth.type == "epan")

  x.len <- length(x)

  if (max.scale %% 1 != 0) {
    stop("max.scale parameter must be an integer.")
  }
  if (max.scale < 1 || max.scale > floor(log2(x.len))) {
    warning("max.scale parameter is outside valid range. Resetting to default value.")
    max.scale <- floor(log2(x.len) * 0.7)
  }
  if (binwidth %% 1 != 0) {
    stop("binwidth parameter must be an integer.")
  }

  J <- wavethresh::IsPowerOfTwo(x.len)

  if (is.na(J) == TRUE) {
    warning("Data length is not power of two. Boundary correction has been applied for spectral estimation.")
    boundary.handle <- TRUE
    dyadic <- FALSE
    J <- floor(log2(x.len)) + 1
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
    x.len = x.len, max.scale = max.scale, boundary.handle = boundary.handle,
    J = J, dyadic = dyadic
  ))
}




#' @title Calculate Boundary Handled Spectrum
#' @description Internal function for spectral estimation used when there is
#' boundary handling
#' @keywords internal
#' @noRd

calc.final.spec <- function(spec, dyadic, x.len) {
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

    final.spec.J <- floor(log2(x.len)) + 1

    final_spec <- wavethresh::cns(2^final.spec.J,
      filter.number = spec$filter$filter.number,
      family = spec$filter$family
    )

    lower <- floor((2^est.spec.J - x.len) / 2)
    upper <- lower + x.len - 1


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

S.calc <- function(x.wd, max.scale, J, inv.mat, filter.number, family) {
  # access smoothed,uncorrected wavelet periodogram:

  uncor.spec <- x.wd$SmoothWavPer

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

  l <- list(S = S, WavPer = x.wd$WavPer, SmoothWavPer = uncor.spec)

  return(l)
}

#' @title Calculate oundary Handled Smoothed Periodogram Estimate
#' @description Internal function for calculating smoothed spectral estimate when
#' boundary handling is used
#' @keywords internal
#' @noRd
smooth.wav.per.calc <- function(x.wd, J, x.len, filter.number, family,
                                dyadic, max.scale) {
  temp <- locits::ewspec3(rep(0, 2^J), filter.number = filter.number, family = family)

  temp$SmoothWavPer <- calc.final.spec(x.wd$SmoothWavPer, dyadic = dyadic, x.len = x.len)
  temp$WavPer <- calc.final.spec(x.wd$WavPer, dyadic = dyadic, x.len = x.len)

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


#' @title Trend Estimation Error Checks
#' @description Internal function for checking wav.trend.est and wav.diff.trend.est
#' @keywords internal
#' @noRd
trend.est.checks <- function(x, max.scale, boundary.handle, transform.type,
                             T.CI, reps, sig.lvl, est.type) {
  stopifnot(
    "Parameter T.transform must be either 'dec' or 'nondec'" =
      transform.type == "dec" || transform.type == "nondec"
  )
  stopifnot("Parameter T.CI must be logical variable" = is.logical(T.CI))
  stopifnot("Parameter T.boundary.handle must be logical variable" = is.logical(boundary.handle))


  if (any(is.na(x))) {
    stop("Data contains mising values.")
  }
  if (!is.numeric(x)) {
    stop("Data is not numeric")
  }

  if (!is.numeric(reps)) {
    stop("Number of bootstrap replications should be a single positive integer.")
  }
  if ((length(reps) != 1) || (reps %% 1 != 0) || (reps < 1)) {
    stop("Number of bootstrap replications should be a single positive integer.")
  }
  stopifnot("Error: sig.lvl must be a number between 0 and 1." = sig.lvl >= 0 && sig.lvl <= 1)


  x.len <- length(x)

  if (max.scale %% 1 != 0) {
    stop("max.scale parameter must be an integer.")
  }
  if (max.scale < 1 || max.scale > floor(log2(x.len))) {
    warning("max.scale parameter is outside valid range. Resetting to default value.")
    max.scale <- floor(log2(x.len) * 0.7)
  }

  J <- wavethresh::IsPowerOfTwo(x.len)

  if (is.na(J) == TRUE) {
    warning("Data length is not power of two. Boundary correction has been applied for trend estimation.")
    boundary.handle <- TRUE
    dyadic <- FALSE
    J <- floor(log2(x.len)) + 1
  } else {
    dyadic <- TRUE
  }


  return(list(
    x.len = x.len, max.scale = max.scale, boundary.handle = boundary.handle,
    J = J, dyadic = dyadic
  ))
}

#' @title Epanechnikov Kernel
#' @description Internal function for Epanechnikov smoothing
#' @keywords internal
#' @noRd
epan.kern.f <- function(tt) {
  sqrt(2) * (1 - (tt^2) / 5)
}

#' @title Epanechnikov Kernel Calculation
#' @description Internal function for calculating Epanechnikov kernel filter
#' @keywords internal
#' @noRd
epan <- function(epan.len) {
  mytt <- seq(from = -sqrt(5), to = sqrt(5), length = epan.len)

  sf <- sum(epan.kern.f(mytt))

  return(epan.kern.f(mytt) / sf)
}


#' @title Plot a spectral estimate
#' @description Internal function for plotting spectral estimate
#' @keywords internal
#' @noRd
spec.plot <- function(x, xlabvals, xlabchars, ylabchars, first.level = 0, n,
                      main = "Spectral Estimate", scaling = "global",
                      rhlab = FALSE, sub = "", NotPlotVal = 0.005, xlab = "Time",
                      ylab = "Scale", aspect = "Identity", ...) {
  ctmp <- class(x)
  if (is.null(ctmp)) {
    stop("wd has no class")
  } else if (ctmp != "wd") {
    stop("wd is not of class wd")
  }

  levels <- wavethresh::nlevelsWT(x)


  nlevels <- levels - first.level
  type <- x$type

  if (aspect != "Identity") {
    sub <- paste(sub, "(", aspect, ")")
  }

  plot(c(0, 0, n, n), c(0, nlevels + 1, nlevels + 1, 0),
    type = "n",
    xlab = xlab, ylab = ylab, main = main, yaxt = "n", xaxt = "n",
    sub = sub, ...
  )

  if (missing(ylabchars)) {
    ylabchars <- 1:nlevels
  }

  axis(2, at = 1:(nlevels), labels = ylabchars)
  if (missing(xlabchars)) {
    if (missing(xlabvals)) {
      axx <- c(
        0, 2^(levels - 2), 2^(levels - 1),
        2^(levels - 1) + 2^(levels - 2), 2^levels
      )
      axis(1, at = axx)
    } else {
      lx <- pretty(xlabvals, n = 4)
      if (lx[1] < min(xlabvals)) {
        lx[1] <- min(xlabvals)
      }
      if (lx[length(lx)] > max(xlabvals)) {
        lx[length(lx)] <- max(xlabvals)
      }
      xix <- NULL
      for (i in 1:length(lx)) {
        u <- (xlabvals - lx[i])^2
        xix <- c(xix, (1:length(u))[u == min(u)])
      }
      axx <- xix
      axl <- signif(lx, digits = 2)
      axis(1, at = axx, labels = axl)
    }
  } else {
    axis(1, at = xlabvals, labels = xlabchars)
  }
  myxx <- 1:n
  height <- 1
  first.last.d <- x$fl.dbase$first.last.d
  axr <- NULL
  if (scaling == "global") {
    my <- 0
    for (i in ((levels - 1):first.level)) {
      y <- wavethresh::accessD(x, i, aspect = aspect)
      my <- max(c(my, abs(y)))
    }
  }
  if (scaling == "compensated") {
    my <- 0
    for (i in ((levels - 1):first.level)) {
      y <- wavethresh::accessD(x, i, aspect = aspect) * 2^(i / 2)
      my <- max(c(my, abs(y)))
    }
  }
  if (scaling == "super") {
    my <- 0
    for (i in ((levels - 1):first.level)) {
      y <- wavethresh::accessD(x, i, aspect = aspect) * 2^i
      my <- max(c(my, abs(y)))
    }
  }
  shift <- 1
  for (i in ((levels - 1):first.level)) {
    y <- wavethresh::accessD(x, i, aspect = aspect)
    if (type == "wavelet") {
      n <- 2^i
    } else {
      y <- y[c((n - shift + 1):n, 1:(n - shift))]
      shift <- shift * 2
    }
    xplot <- myxx
    ly <- length(y)
    if (scaling == "by.level") {
      my <- max(abs(y))
    }
    if (scaling == "compensated") {
      y <- y * 2^(i / 2)
    }
    if (scaling == "super") {
      y <- y * 2^i
    }
    if (my == 0) {
      y <- rep(0, length(y))
    } else {
      y <- (0.5 * y) / my
    }
    axr <- c(axr, my)
    if (max(abs(y)) > NotPlotVal) {
      segments(xplot, height, xplot, height + y)
    }
    if (i != first.level) {
      if (type == "wavelet") {
        x1 <- myxx[seq(1, n - 1, 2)]
        x2 <- myxx[seq(2, n, 2)]
        myxx <- (x1 + x2) / 2
      }
      height <- height + 1
    }
  }
  if (rhlab == TRUE) {
    axis(4, at = 1:length(axr), labels = signif(axr, digits = 3))
  }
}

#' @title Convert matrix to wd object
#' @description Internal function for plotting spectral estimate
#' @keywords internal
#' @noRd

mat.to.spec <- function(s.mat, filter.number = 1, family = "DaubExPhase") {
  J <- nrow(s.mat)

  spec <- wavethresh::cns(2^J, filter.number = filter.number, family = family)

  for (j in 1:J) {
    spec <- wavethresh::putD(spec, level = J - j, s.mat[j, ])
  }

  spec
}


#' @title LACF calculation
#' @description Internal function for estimating lacf, used inside TLSW function
#' for calculating confidence intervals if required.
#' @keywords internal
#' @noRd

TLSW.TLSWlacf <- function(x, filter.number = 4, family = "DaubExPhase",
                          spec.est = NULL, lag.max = NULL, ...) {
  stopifnot("Paramter lag.max should be a nonegative integer." = lag.max >= 0)

  if (is.null(spec.est)) {
    spec.est <- ewspec.trend(
      x = x, an.filter.number = filter.number, an.family = family,
      gen.filter.number = filter.number, gen.family = family, ...
    )
  }

  dsname <- deparse(substitute(x))

  S <- spec.est$S
  SmoothWP <- spec.est$SmoothWavPer

  J <- S$nlevels
  Smat <- matrix(S$D, nrow = 2^J, ncol = J)[1:length(x), ]
  Psi <- wavethresh::PsiJmat(-J, filter.number = filter.number, family = family)
  nc <- ncol(Psi)
  L <- (nc - 1) / 2
  dimnames(Psi) <- list(NULL, c(-L:0, 1:L))
  if (is.null(lag.max)) {
    lag.max <- floor(10 * (log10(length(x))))
  }
  if (L + 1 + lag.max > ncol(Psi)) {
    warning(paste(
      "lag.max too high. Have reset it to ",
      ncol(Psi) - L - 1, ". Higher lags are zero"
    ))
    lag.max <- ncol(Psi) - L - 1
  }
  the.lacf <- Smat %*% Psi[, (L + 1):(L + 1 + lag.max)]
  the.lacor <- sweep(the.lacf, 1, the.lacf[, 1], FUN = "/")
  l <- list(
    lacf = the.lacf, lacr = the.lacor, name = dsname,
    date = date(), SmoothWP = SmoothWP, S = S, J = J
  )
  class(l) <- "lacf"
  return(l)
}


#' @title Wavelet Periodogram smoothing
#' @description Internal function for smoothing the wavelet periodogram via
#' Epanechnikov or median smoothing.
#' @keywords internal
#' @noRd

WP.manual.smooth <- function(x.wd, smooth.type, max.scale, binwidth) {
  J2 <- wavethresh::nlevelsWT(x.wd$SmoothWavPer)

  if (smooth.type == "median") {
    for (j in 1:max.scale) {
      x.wd$SmoothWavPer <- suppressWarnings(wavethresh::putD(x.wd$SmoothWavPer,
        level = J2 - j,
        2.125 * stats::runmed(wavethresh::accessD(x.wd$WavPer, level = J2 - j), k = binwidth)
      ))
    }
  }
  if (smooth.type == "epan") {
    epan.filter <- epan(binwidth)
    for (j in 1:max.scale) {
      temp.dj <- wavethresh::accessD(x.wd$WavPer, level = J2 - j)
      temp.dj <- c(rev(temp.dj[1:(floor((binwidth - 1) / 2))]), temp.dj, rev(temp.dj[(length(temp.dj) - floor(binwidth / 2) + 1):length(temp.dj)]))

      temp <- stats::filter(temp.dj, epan.filter)
      temp <- temp[!is.na(temp)]
      x.wd$SmoothWavPer <- wavethresh::putD(x.wd$SmoothWavPer, level = J2 - j, temp)
    }
  }

  return(x.wd)
}
