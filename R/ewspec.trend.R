ewspec.trend = function(data, an.filter.number  = 10, an.family = "DaubLeAsymm",
                        gen.filter.number = an.filter.number, gen.family = an.family,
                        binwidth = floor(2*sqrt(length(data))),
                        max.scale = floor(log2(length(data))*0.7), WP.smooth = TRUE,
                        AutoReflect = TRUE, supply.mat = FALSE, mat = NULL,
                        boundary.handle = TRUE){

  #function that computes the spectral estimate of a time series that has a smooth trend
  #that can be zeroed out by the wavelet coefficients.

  #user chooses a maximum scale of the wavelet transform to analyse, and
  #binwidth of the running mean smoother.

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


  #calculate the appropriate correction matrix and its inverse:

  if (supply.mat==FALSE){
    if (an.filter.number == gen.filter.number && an.family==gen.family){
      A = wavethresh::ipndacw(J = -max.scale, filter.number = an.filter.number, family = an.family)
      inv = solve(A)
    }
    else{
      C = Cmat.calc(J=max.scale, gen.filter.number = gen.filter.number,
                     an.filter.number = an.filter.number, gen.family = gen.family,
                     an.family = an.family)
      inv = solve(C)
    }
  }
  else{
    inv = mat
  }

  if (boundary.handle ==TRUE){
    data = get.boundary.timeseries(data)
  }
  #calculate raw wavelet periodogram which we need to correct:

  data.wd = locits::ewspec3(data, filter.number = an.filter.number, family = an.family,
                    binwidth = binwidth, AutoReflect = AutoReflect, WPsmooth = WP.smooth)

  calc.final.spec = function(spec, dyadic,data.len){

    if(dyadic==TRUE){
      final_spec = cns(2^(spec$nlevels-2), filter.number = spec$filter$filter.number,
                       family = spec$filter$family)

      lower = 2^(spec$nlevels-2)+2^(spec$nlevels-3)+1
      upper = 2^(spec$nlevels-1)+2^(spec$nlevels-3)


      for (j in 0:(spec$nlevels-3)){

        bh_d = accessD(spec,level = j+2)[lower:upper]

        final_spec = putD(final_spec,level = j, bh_d)


      }

      return(final_spec)
    } else{

      est.spec.J = (spec$nlevels)

      final.spec.J = (floor(log2(data.len))+1)

      final_spec = cns(2^final.spec.J, filter.number = spec$filter$filter.number,
                       family = spec$filter$family)

      lower = floor((2^est.spec.J-data.len)/2)
      upper = lower+data.len-1


      for (j in 0:(final.spec.J-1)){

        bh_d = c(accessD(spec,level = j+(est.spec.J-final.spec.J))[lower:upper], rep(0,2^final.spec.J+lower-upper-1))

        final_spec = putD(final_spec,level = j, bh_d)


      }

      return(final_spec)

    }
  }


  #access smoothed, uncorrected wavelet periodogram:

  if (boundary.handle==TRUE){

   temp = locits::ewspec3(rep(0,2^J),filter.number = an.filter.number, family = an.family)

   temp$SmoothWavPer = calc.final.spec(data.wd$SmoothWavPer,dyadic = dyadic,data.len = data.len)
   temp$WavPer = calc.final.spec(data.wd$WavPer,dyadic = dyadic,data.len = data.len)

    if(max.scale<J){
      for (j in 0:(J-1-max.scale)){
          temp$WavPer = wavethresh::putD(temp$WavPer,level = j, rep(0,2^J))
          temp$SmoothWavPer = wavethresh::putD(temp$SmoothWavPer,level = j, rep(0,2^J))
      }
    }


    data.wd = temp

  }

  uncor.spec = data.wd$SmoothWavPer

  #perform correction: mutiply by inverse matrix, non-estimated scales are set to zero.

  uncor.spec.mat = matrix(0,nrow = max.scale, ncol = 2^J)

  for (j in 1:max.scale){
    uncor.spec.mat[j,] = accessD(uncor.spec,level = J-j)
  }

  #perform correction step:

  cor.spec.mat = inv%*%uncor.spec.mat

  #now fill in wd object with final spectrum estimate.

  S = wavethresh::cns(2^J, filter.number = gen.filter.number,family = gen.family)

  for (j in 1:max.scale){
    S = wavethresh::putD(S, level = J-j, cor.spec.mat[j,])
  }

  #return final estimate, along with smoothed and unsmoothed periodogram:

  l = list(S = S, WavPer = data.wd$WavPer, SmoothWavPer = uncor.spec)

  return(l)

}

