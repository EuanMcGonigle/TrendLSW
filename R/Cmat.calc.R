#' @title Cross Autocorrelation Wavelet Inner Product Matrix Calculation
#' @description Internal function to compute the cross autocorrelation matrix of inner products.
#' This is not intended for general use by regular users of the package.
#' @details Computes the cross inner product matrix of the discrete
#' non-decimated autocorrelation wavelets. This matrix is used to correct the
#' wavelet periodogram analysed using a different wavelet to the wavelet that
#' is assumed to generate the time series. The matrix returned is the one
#' denoted \eqn{C^{(0,1)}} in McGonigle et al. (2022).
#' @param J The dimension of the matrix required. Should be a positive integer.
#' @param gen.filter.number The index of the generating wavelet used to compute
#' the inner product matrix.
#' @param an.filter.number The index of the analysing wavelet used to compute
#' the inner product matrix.
#' @param gen.family The family of generating wavelet used to compute the inner
#' product matrix.
#' @param an.family The family of analysing wavelet used to compute the inner
#' product matrix.
#' @return A J-dimensional square matrix giving the cross inner product
#' autocorrelation wavelet matrix.
#' @seealso \link{TLSW}
#' @references McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
#' locally stationary wavelet processes. \emph{Journal of Time Series
#' Analysis}, 43(6), 895-917.
#' @keywords internal
Cmat.calc <- function(J, gen.filter.number = 1, an.filter.number = 1,
                      gen.family = "DaubExPhase", an.family = "DaubExPhase") {
  if (!is.numeric(J)) {
    stop("The parameter J should be a positive integer.")
  }
  if ((length(J) != 1) || (J %% 1 != 0) || (J < 1)) {
    stop("The parameter J should be a positive integer.")
  }


  C_mat <- matrix(0, nrow = J, ncol = J)

  h.gen <- wavethresh::filter.select(filter.number = gen.filter.number, family = gen.family)$H # access the filter for chosen wavelet

  h.gen.len <- as.numeric(length(h.gen))

  h.an <- wavethresh::filter.select(filter.number = an.filter.number, family = an.family)$H # access the filter for chosen wavelet

  h.an.len <- as.numeric(length(h.an))

  swap <- FALSE

  # use the wavethresh function to calculate the autocorrelation wave matrix,
  # subsetted to get the same form as above calculated

  new_Psi_mat1 <- wavethresh::PsiJmat(J = -J - 1, filter.number = gen.filter.number, family = gen.family)[1:J, ]

  new_Psi_mat2 <- wavethresh::PsiJmat(J = -J - 1, filter.number = an.filter.number, family = an.family)[1:J, ]

  if (h.gen.len >= h.an.len) {
    mod_psi1 <- new_Psi_mat1
    mod_psi2 <- matrix(0, nrow = nrow(mod_psi1), ncol = ncol(mod_psi1))
    psi1_len <- length(new_Psi_mat1[J, ])
    psi2_len <- length(new_Psi_mat2[J, ])

    midpoint <- (psi1_len + 1) / 2

    mod_psi2[1:J, (midpoint - (psi2_len - 1) / 2):(midpoint + (psi2_len - 1) / 2)] <- new_Psi_mat2
  } else {
    mod_psi1 <- new_Psi_mat2
    swap <- TRUE
    mod_psi2 <- matrix(0, nrow = nrow(mod_psi1), ncol = ncol(mod_psi1))
    psi1_len <- length(new_Psi_mat2[J, ])
    psi2_len <- length(new_Psi_mat1[J, ])

    midpoint <- (psi1_len + 1) / 2

    mod_psi2[1:J, (midpoint - (psi2_len - 1) / 2):(midpoint + (psi2_len - 1) / 2)] <- new_Psi_mat1
  }

  # all the code is doing is making sure the inner product is calculated right
  # we sit the small wave matrix inside the big one so that the multiplications
  # will line up

  for (row in 1:J) {
    for (column in 1:J) {
      C_mat[row, column] <- sum(mod_psi1[row, ] * mod_psi2[column, ])
    }
  }

  if (swap == TRUE) {
    C_mat <- t(C_mat)
  }

  return(C_mat)
}
