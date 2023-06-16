#' @title Lagged Autocorrelation Wavelet Inner Product Calculation
#' @description Computes the matrix of lagged autocorrelation wavelet inner
#' products.
#' @details Computes the lagged inner product matrix of the discrete
#' non-decimated autocorrelation wavelets. This matrix is used in the
#' calculation to correct the wavelet periodogram of the differenced time
#' series. With \code{lag} \eqn{= \tau}, the matrix returned is the one called \eqn{A^\tau} in McGonigle et al. (2022).
#' @param J The dimension of the matrix required. Should be a positive integer.
#' @param filter.number The index of the wavelet used to compute the inner
#' product matrix.
#' @param family The family of wavelet used to compute the inner product
#' matrix.
#' @param lag The lag of matrix to calculate. A lag of 0 corresponds to the
#' standard A matrix.
#' @return A J-dimensional square matrix giving the lagged inner product
#' autocorrelation wavelet matrix.
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
#' time-varying first and second-order structure of time series via wavelets
#' and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
#' @examples
#' A1 <- Atau.mat.calc(J = 5, filter.number = 1, family = "DaubExPhase", lag = 1)
#' @export
#' @seealso \link{ewspec.diff}
Atau.mat.calc <- function(J, filter.number = 1, family = "DaubExPhase", lag = 1) {

  lagged_A_mat <- matrix(0, nrow = J, ncol = J)

  new_Psi_mat <- wavethresh::PsiJmat(J = -J - 1, filter.number = filter.number, family = family)[1:J, ]

  for (row in 1:J) {
    for (column in row:J) {
      lagged_A_mat[row, column] <- sum(new_Psi_mat[row, ] * do.shift(new_Psi_mat[column, ], lag))
      lagged_A_mat[column, row] <- lagged_A_mat[row, column]
    }
  }

  return(lagged_A_mat)
}
