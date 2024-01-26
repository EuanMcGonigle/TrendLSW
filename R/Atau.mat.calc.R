#' @title Lagged Autocorrelation Wavelet Inner Product Calculation
#' @description Internal function for computing the matrix of lagged autocorrelation wavelet inner
#' products. This is not intended for general use by regular users of the package.
#' @details Computes the lagged inner product matrix of the discrete
#' non-decimated autocorrelation wavelets. This matrix is used in the
#' calculation to correct the wavelet periodogram of the differenced time
#' series. With \code{lag} \eqn{= \tau}, the matrix returned is the matrix \eqn{A^\tau} in McGonigle et al. (2022).
#' @param J The dimension of the matrix required. Should be a positive integer.
#' @param filter.number The index of the wavelet used to compute the inner
#' product matrix.
#' @param family The family of wavelet used to compute the inner product
#' matrix.
#' @param lag The lag of matrix to calculate. A lag of 0 corresponds to the
#' matrix \eqn{A} defined in Nason et al. (2000).
#' @return A J-dimensional square matrix giving the lagged inner product
#' autocorrelation wavelet matrix.
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#'
#' Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and
#' adaptive estimation of the evolutionary wavelet spectrum. \emph{Journal of
#' the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{62(2)}, 271--292.
#' @seealso \link{TLSW}
#' @keywords internal
Atau.mat.calc <- function(J, filter.number = 1, family = "DaubExPhase", lag = 1) {
  if (!is.numeric(lag)) {
    stop("The lag parameter should be a positive integer.")
  }
  if ((length(lag) != 1) || (lag %% 1 != 0) || (lag < 1)) {
    stop("The lag parameter should be a positive integer.")
  }
  if (!is.numeric(J)) {
    stop("The parameter J should be a positive integer.")
  }
  if ((length(J) != 1) || (J %% 1 != 0) || (J < 1)) {
    stop("The parameter J should be a positive integer.")
  }

  lagged_A_mat <- matrix(0, nrow = J, ncol = J)

  Psi_mat <- wavethresh::PsiJmat(J = -J, filter.number = filter.number, family = family)

  Psi_mat <- cbind(Psi_mat, matrix(0, nrow = J, ncol = lag))

  P.ncol <- ncol(Psi_mat)

  for (row in 1:J) {
    for (column in row:J) {
      lagged_A_mat[row, column] <- sum(Psi_mat[row, 1:(P.ncol - lag)] * Psi_mat[column, (lag + 1):P.ncol])
      lagged_A_mat[column, row] <- lagged_A_mat[row, column]
    }
  }

  return(lagged_A_mat)
}
