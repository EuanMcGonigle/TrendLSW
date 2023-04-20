ewspec.diff <- function(data, lag = 1, filter.number = 1, family = "DaubExPhase",
                        binwidth = floor(2 * sqrt(length(data))), diff.number = 1,
                        max.scale = floor(log2(length(data)) * 0.7), WP.smooth = TRUE,
                        boundary.handle = FALSE, AutoReflect = FALSE,
                        supply.inv.mat = FALSE, inv = NULL) {

  # function that computes the spectral estimate of a time series that has a trend.

  # user chooses a maximum scale of the wavelet transform to analyse,
  # binwidth of the running mean smoother, and the method of differencing.


  if (sum(is.na(data)) != 0) {
    stop("Data contains mising values.")
  }
  if (!is.atomic(data)) {
    stop("Data is not atomic")
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

  if (boundary.handle == TRUE) {
    data <- get.boundary.timeseries(data, type = "LSW.diff")
  }

  # difference data to remove trend/seasonality and calculate the appropriate correction matrix
  # and its inverse:

  if (diff.number == 1) {
    diff.data <- c(diff(data, lag))
    diff.data <- c(diff.data, rep(0, lag))
  } else if (diff.number == 2) {
    diff.data <- c(diff(diff(data)), 0, 0)
  } else {
    diff.number <- 1
    diff.data <- c(diff(data, lag))
    diff.data <- c(diff.data, rep(0, lag))
    warning("Function only implements 1st or 2nd differences. Using 1st difference.")
  }

  if (supply.inv.mat == FALSE) {
    A <- wavethresh::ipndacw(J = -max.scale, filter.number = filter.number, family = family)
    A1 <- Atau.mat.calc(J = max.scale, filter.number = filter.number, family = family, lag = lag)

    if (diff.number == 1) {
      inv <- solve(2 * A - 2 * A1)
    } else if (diff.number == 2) {
      A2 <- Atau.mat.calc(J = max.scale, filter.number = filter.number, family = family, lag = 2)
      inv <- solve(6 * A - 8 * A1 + 2 * A2)
    }
  }

  calc.final.spec <- function(spec, dyadic, data.len) {
    if (dyadic == TRUE) {
      final_spec <- wavethresh::cns(2^(spec$nlevels - 2),
        filter.number <- spec$filter$filter.number,
        family <- spec$filter$family
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
        filter.number <- spec$filter$filter.number,
        family <- spec$filter$family
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



  # calculate raw wavelet periodogram which we need to correct:

  data.wd <- locits::ewspec3(diff.data,
    filter.number = filter.number, family = family,
    binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = WP.smooth)

  if (boundary.handle == TRUE) {
    temp <- locits::ewspec3(rep(0, 2^J), filter.number = filter.number, family = family)

    temp$SmoothWavPer <- calc.final.spec(data.wd$SmoothWavPer, dyadic = dyadic, data.len = data.len)
    temp$WavPer <- calc.final.spec(data.wd$WavPer, dyadic = dyadic, data.len = data.len)

    if (max.scale < J) {
      for (j in 0:(J - 1 - max.scale)) {
        temp$WavPer <- wavethresh::putD(temp$WavPer, level = j, rep(0, 2^J))
        temp$SmoothWavPer <- wavethresh::putD(temp$SmoothWavPer, level = j, rep(0, 2^J))
      }
    }


    data.wd <- temp
  }

  # access smoothed,uncorrected wavelet periodogram:

  uncor.spec <- data.wd$SmoothWavPer

  # perform correction: mutiply by inverse matrix, non-estimated scales are set to zero.

  uncor.spec.mat <- matrix(0, nrow = max.scale, ncol = 2^J)

  for (j in 1:max.scale) {
    uncor.spec.mat[j, ] <- wavethresh::accessD(uncor.spec, level = J - j)
  }

  # perform correction step:

  cor.spec.mat <- inv %*% uncor.spec.mat

  # now fill in wd object with final spectrum estimate.

  S <- wavethresh::cns(2^J, filter.number = filter.number, family = family)

  for (j in 1:max.scale) {
    S <- wavethresh::putD(S, level = J - j, cor.spec.mat[j, ])
  }

  # return final EWS estimate, along with smoothed and unsmoothed periodogram:

  l <- list(S = S, WavPer = data.wd$WavPer, SmoothWavPer = uncor.spec)

  return(l)
}
