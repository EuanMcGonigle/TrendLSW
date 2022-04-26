bootstrap.trendCI = function(trend.est, spec, var.mat, nsims=100){

  trend.mat = matrix(0, nrow= nsims, ncol = length(trend.est))

  spec$D[spec$D<0] = 0

  for (i in 1:nsims){

    x = trend.est+LSWsim(spec)

    x.trend = wav.diff.trend.est(x,spec, policy = "manual", Tmat = var.mat)$trend.est

    trend.mat[i,] = x.trend

    cat(".")
  }

  trend.sd = rep(0, 1024)

  for (i in 1:1024){

    trend.sd[i] = sd(trend.mat[,i])

  }

  return(list(trend.sd = trend.sd, trend.mat = trend.mat))

}






