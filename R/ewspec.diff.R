ewspec.diff = function(data, lag=1, filter.number  = 1, family = "DaubExPhase",
                       binwidth = floor(2*sqrt(length(data))), diff.number = 1,
                       max.scale = floor(log2(length(data))*0.7), WP.smooth = TRUE,
                       AutoReflect = FALSE, supply.inv.mat = FALSE, inv = NULL){

  #require(wavethresh)
  #require(locits)

  #function that computes the spectral estimate of a time series that has a trend.

  #user chooses a maximum scale of the wavelet transform to analyse,
  #binwidth of the running mean smoother, and the method of differencing.


  if (sum(is.na(data))!=0){
    stop("Data contains mising values.")
  }
  if (!is.atomic(data)){
    stop("Data is not atomic")
  }

  data.len = length(data)

  if(max.scale%%1!=0){
    stop("max.scale parameter must be an integer.")
  }
  if(max.scale < 1 || max.scale > floor(log2(data.len))){
    warning("max.scale parameter is outside valid range. Resetting to default value.")
    max.scale = floor(log2(data.len)*0.7)
  }
  if(binwidth%%1!=0){
    stop("binwidth parameter must be an integer.")
  }

  J = wavethresh::IsPowerOfTwo(data.len)

  if(is.na(J)==TRUE){
    warning("Data length is not power of two. Boundary correction has been applied.")
    boundary.handle=TRUE
    dyadic=FALSE
    J = floor(log2(data.len))+1
  } else{
    dyadic=TRUE
  }

  #difference data to remove trend/seasonality and calculate the appropriate correction matrix
  #and its inverse:

  if (diff.number==1){
    diff.data = c(diff(data,lag))
    diff.data = c(diff.data, rep(0,lag))
  } else if(diff.number==2){
    diff.data = c(diff(diff(data)),0,0)
  } else{
    diff.number=1
    diff.data = c(diff(data,lag))
    diff.data = c(diff.data, rep(0,lag))
    warning("Function only implements 1st or 2nd differences. Using 1st difference.")
  }

  if(supply.inv.mat ==FALSE){

    A = wavethresh::ipndacw(J = -max.scale, filter.number = filter.number, family = family)
    A1 = Atau.mat.calc(J=max.scale, filter.number = filter.number, family = family,lag = lag)

    if (diff.number==1){
      inv = solve(2*A-2*A1)
    }
    else if (diff.number==2){
      A2 = Atau.mat.calc(J=max.scale, filter.number = filter.number, family = family,lag = 2)
      inv = solve(6*A-8*A1+2*A2)
    }
  }


  #calculate raw wavelet periodogram which we need to correct:

  data.wd = locits::ewspec3(diff.data, filter.number = filter.number, family = family,
                    binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = WP.smooth)

  #access smoothed,uncorrected wavelet periodogram:

  uncor.spec = data.wd$SmoothWavPer

  #perform correction: mutiply by inverse matrix, non-estimated scales are set to zero.

  uncor.spec.mat = matrix(0,nrow = max.scale, ncol = data.len)

  for (j in 1:max.scale){
    uncor.spec.mat[j,] = accessD(uncor.spec,level = J-j)
  }

  #perform correction step:

  cor.spec.mat = inv%*%uncor.spec.mat

  #now fill in wd object with final spectrum estimate.

  S = cns(data.len, filter.number = filter.number,family = family)

  for (j in 1:max.scale){
    S = wavethresh::putD(S, level = J-j, cor.spec.mat[j,])
  }

  #return final EWS estimate, along with smoothed and unsmoothed periodogram:

  l = list(S = S, WavPer = data.wd$WavPer, SmoothWavPer = uncor.spec)

  return(l)

}
