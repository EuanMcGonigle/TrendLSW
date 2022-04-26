create.covmat = function(lacf){

  lacf = lacf$lacf

  data_len = length(lacf[,1])

  max_lag = length(lacf[1,])

  cov_mat = matrix(0,nrow = data_len,ncol = data_len)

  for (row in 1:data_len){
    for (column in row:(min((max_lag-1+row),data_len))){
      cov_mat[row,column] = lacf[row,(abs(column-row)+1)]
      cov_mat[column,row] = cov_mat[row,column]
    }
  }

  cov_mat

}
