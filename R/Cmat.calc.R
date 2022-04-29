Cmat.calc = function(J, gen.filter.number = 1, an.filter.number = 1,
                      gen.family = "DaubExPhase", an.family = "DaubExPhase"){

  # function that calculates the cross autocorrelation matrix C^(1,0)

  C_mat = matrix(0, nrow = J, ncol = J)

  h.gen = filter.select(filter.number = gen.filter.number, family = gen.family)$H #access the filter for chosen wavelet

  h.gen.len = as.numeric(length(h.gen))

  h.an = filter.select(filter.number = an.filter.number, family = an.family)$H #access the filter for chosen wavelet

  h.an.len = as.numeric(length(h.an))

  swap = FALSE

  #use the wavethresh function to calculate the autocorrelation wave matrix,
  #subsetted to get the same form as above calculated

  new_Psi_mat1 = wavethresh::PsiJmat(J = -J-1, filter.number = gen.filter.number, family = gen.family)[1:J,]

  new_Psi_mat2 = wavethresh::PsiJmat(J = -J-1, filter.number = an.filter.number, family = an.family)[1:J,]

  if (h.gen.len>=h.an.len){
    mod_psi1 = new_Psi_mat1
    mod_psi2 = matrix(0,nrow = nrow(mod_psi1), ncol = ncol(mod_psi1))
    psi1_len = length(new_Psi_mat1[J,])
    psi2_len = length(new_Psi_mat2[J,])

    midpoint = (psi1_len+1)/2

    mod_psi2[1:J,(midpoint-(psi2_len-1)/2):(midpoint+(psi2_len-1)/2)] = new_Psi_mat2
  }
  else{
    mod_psi1 = new_Psi_mat2
    swap = TRUE
    mod_psi2 = matrix(0,nrow = nrow(mod_psi1), ncol = ncol(mod_psi1))
    psi1_len = length(new_Psi_mat2[J,])
    psi2_len = length(new_Psi_mat1[J,])

    midpoint = (psi1_len+1)/2

    mod_psi2[1:J,(midpoint-(psi2_len-1)/2):(midpoint+(psi2_len-1)/2)] = new_Psi_mat1

  }

  #all the code is doing is making sure the inner product is calculated right
  #we sit the small wave matrix inside the big one so that the multiplications
  #will line up

  for (row in 1:J){
    for (column in 1:J){
      C_mat[row,column] = sum(mod_psi1[row,]*mod_psi2[column,])
    }
  }

  if (swap==TRUE){
    C_mat = t(C_mat)
  }

  return(C_mat)

}
