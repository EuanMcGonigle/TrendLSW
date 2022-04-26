trend.estCI = function(trend_est, cov_mat, filter.number = 10, family = "DaubLeAsymm", alpha = 0.95){

  #function to create confidence interval for the trend estimate

  data_len = length(trend_est)

  siz <- 1 - alpha
  qval <- qnorm(1 - siz/2)

  #calculate the wavelet transform matrix:

  W = t(GenW(n = data_len,filter.number = filter.number,family = family))

  scale = log2(data_len)

  boundary_test = c(rep(0,data_len-1),1)

  y_wd = wd(boundary_test,family = family, filter.number = filter.number)

  boundary_coefs = list(1:scale)

  for (i in (scale-1):1){

    temp = accessD(y_wd,level = i)

    boundary_coefs[[i]] = (which(temp!=0))

  }

  boundary_vec = c(1,rep(0,data_len-2),1)

  for(i in (scale-1):3){

    boundary_vec[(data_len-2^(i+1)+2):(data_len-2^(i+1)+1+2^i)][boundary_coefs[[i]]] = 1

  }


  #calculate the linear thresholding operation matrix:

  A = diag(boundary_vec)

  #calcualte the covariance matrix of the trend estimate:

  R = t(W)%*%A%*%W

  trend_cov = (R)%*%cov_mat%*%t(R)

  trend_var = diag(trend_cov)

  l = list(trend.est = trend_est, trend.var = trend_var, upper.conf = trend_est+qval*trend_var,
           lower.conf = trend_est-qval*trend_var)

  return(l)

}
