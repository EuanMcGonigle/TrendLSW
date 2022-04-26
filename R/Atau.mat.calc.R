Atau.mat.calc = function(J, filter.number = 1, family = "DaubExPhase", lag=1){

  # function that calculates the lagged A matrix, A^tau, the innter product matrix at lag tau

  #require(binhf)
  #require(wavethresh)

  #initialise matrix:

  lagged_A_mat = matrix(0, nrow = J, ncol = J)


  #use the wavethresh function to calculate the autocorrelation wave matrix,
  #subsetted to get the same form as above calculated, which is required
  #to do the final calculation

  new_Psi_mat = PsiJmat(J = -J-1, filter.number = filter.number, family = family)[1:J,]

  for (row in 1:J){
    for (column in row:J){
      lagged_A_mat[row,column] = sum(new_Psi_mat[row,]*shift(new_Psi_mat[column,],lag))
      lagged_A_mat[column,row] = lagged_A_mat[row,column]
    }
  }

  return(lagged_A_mat)

}
