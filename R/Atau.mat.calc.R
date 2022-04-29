Atau.mat.calc = function(J, filter.number = 1, family = "DaubExPhase", lag=1){

  # function that calculates the lagged A matrix, A^tau, the innter product matrix at lag tau

  #initialise matrix:

  lagged_A_mat = matrix(0, nrow = J, ncol = J)

  #use the wavethresh function to calculate the autocorrelation wave matrix,
  #subsetted to get the same form as above calculated, which is required
  #to do the final calculation
  do.shift = function (v, places, dir = "right")
  {
    vnew <- NULL
    d <- substring(dir, 1, 1)
    if (d == "r" & places < 0) {
      d <- "l"
    }
    else {
      if (d == "l" & places < 0) {
        d <- "r"
      }
    }
    n <- length(v)
    if (n == 1) {
      places <- 0
    }
    p <- abs(places)
    if (p == 0) {
      vnew <- v
    }
    else {
      if (d == "r") {
        vnew <- c(v[(n - p + 1):n], v[1:(n - p)])
      }
      else {
        vnew <- c(v[(p + 1):n], v[1:p])
      }
    }
    vnew
  }

  new_Psi_mat = wavethresh::PsiJmat(J = -J-1, filter.number = filter.number, family = family)[1:J,]

  for (row in 1:J){
    for (column in row:J){
      lagged_A_mat[row,column] = sum(new_Psi_mat[row,]*do.shift(new_Psi_mat[column,],lag))
      lagged_A_mat[column,row] = lagged_A_mat[row,column]
    }
  }

  return(lagged_A_mat)

}
