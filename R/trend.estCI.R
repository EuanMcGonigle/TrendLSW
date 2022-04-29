trend.estCI = function(trend.est, lacf.est, filter.number = 4, family = "DaubLeAsymm", alpha = 0.95){

  #function to create confidence interval for the trend estimate

  data.len = length(trend.est)

  size <- 1 - alpha
  qval <- qnorm(1 - size/2)

  create.covmat = function(lacf, data.len){

    lacf = lacf$lacf

    max.lag = length(lacf[1,])

    cov.mat = matrix(0,nrow = data.len,ncol = data.len)

    for (row in 1:data.len){
      for (column in row:(min((max.lag-1+row),data.len))){
        cov.mat[row,column] = lacf[row,(abs(column-row)+1)]
        cov.mat[column,row] = cov.mat[row,column]
      }
    }

    cov.mat

  }

  cov.mat = create.covmat(lacf.est, data.len)

  #calculate the wavelet transform matrix:

  W = t(wavethresh::GenW(n = data.len,filter.number = filter.number,family = family))

  scale = log2(data.len)

  boundary_test = c(rep(0,data.len-1),1)

  y_wd = wd(boundary_test,family = family, filter.number = filter.number)

  boundary_coefs = list(1:scale)

  for (i in (scale-1):1){

    temp = accessD(y_wd,level = i)

    boundary_coefs[[i]] = (which(temp!=0))

  }

  boundary_vec = c(1,rep(0,data.len-2),1)

  for(i in (scale-1):3){

    boundary_vec[(data.len-2^(i+1)+2):(data.len-2^(i+1)+1+2^i)][boundary_coefs[[i]]] = 1

  }


  #calculate the linear thresholding operation matrix:

  A = diag(boundary_vec)

  #calcualte the covariance matrix of the trend estimate:

  R = t(W)%*%A%*%W

  trend.cov = (R)%*%cov.mat%*%t(R)

  trend.sd = suppressWarnings(sqrt(diag(trend.cov)))

  v.nas = which(is.na(trend.sd))

  trend.sd.final = trend.sd

  if(length(v.nas>0)){
    for(i in 1:length(v.nas)){
      trend.sd.final[v.nas[i]] = trend.sd[which(!is.nan(trend.sd))][which.min(abs(v.nas[i]-which(!is.nan(trend.sd))))]
    }
  }

  l = list(trend.est = trend.est, trend.var = trend.sd.final^2, lower.conf = trend.est-trend.sd.final*qval, upper.conf = trend.est+trend.sd.final*qval)

  return(l)

}
