wav.diff.trend.est = function(data, spec, filter.number = 4, thresh.type = "soft",
                              normal = TRUE,family = "DaubLeAsymm",
                              max.scale = floor(0.7*log2(length(data))),
                               boundary.handle = FALSE){


  data.len = length(data)

  if(max.scale%%1!=0){
    stop("max.scale parameter must be an integer.")
  }
  if(max.scale < 1 || max.scale > floor(log2(data.len))){
    warning("max.scale parameter is outside valid range. Resetting to default value.")
    max.scale = floor(log2(data.len)*0.7)
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

  #by default, do T.I. denoising:

  if (boundary.handle == TRUE){
    orig.data = data
    data = get.boundary.timeseries(data)
  }

  data.len = length(data)
  J = wavethresh::IsPowerOfTwo(data.len)


  data.wd = wavethresh::wd(data, filter.number = filter.number, family = family,type = "station")

  #calculate C to use in variance estimate of wavelet coefficients

  C = Cmat.calc(J=max.scale, an.filter.number = filter.number,
                gen.filter.number = spec$filter$filter.number,
                an.family = family, gen.family = spec$filter$family)


  replace.neg.values = function(var.mat,max.scale){

    for (j in 1:max.scale){

      var.row = var.mat[j,]

      if(sum(var.row<=0)>0){

        var.row[var.row<0] = 0
        var.row0 = which(var.row==0)
        var.row.non0 = var.row[which(var.row!=0)]

        for(i in 1:length(var.row0)){
          var.mat[j,(var.row0[i])] = var.row.non0[which.min(abs(var.row0[i]-which(var.row!=0)))]
        }
      }
    }

    var.mat

  }
  #find negative values and replace them with nearby positive values:


  spec.mat = matrix(0,nrow = max.scale,ncol = length(accessD(spec,level = 0)))

  for (j in 1:max.scale){
    spec.mat[j,] = accessD(spec,level = nlevelsWT(spec)-j)
  }

  #calculate variance estimate matrix

  var.mat = C%*%spec.mat

  var.mat = replace.neg.values(var.mat,max.scale)

  data.wd.thresh = data.wd


  if (boundary.handle==TRUE){

    #below code determines the boundary coefficients for a given wavelet

    boundary.test = c(rep(0,data.len-1),1)

    y.wd = wavethresh::wd(boundary.test,family = family, filter.number = filter.number,type = "station")

    #create EWS estimate matrix

    bc.var.mat = matrix(0,nrow = max.scale, ncol = data.len)

    lower = floor((data.len-length(orig.data))/2)
    upper = lower+length(orig.data)-1

    bc.var.mat[,lower:upper] = var.mat[,1:length(orig.data)]

    #boundary.coefs = list()

    #for (i in 1:max.scale){

      #temp = accessD(y.wd,level = nlevelsWT(y.wd)-i)

      #boundary.coefs[[i]] = (which(temp!=0))

      #bc.var.mat[i,boundary.coefs[[i]]] = 0
    #}



    for (j in 1:max.scale){

      dj = accessD(data.wd,level = J-j)

      if (normal==TRUE){
        thresh = sqrt(2*bc.var.mat[j,]*log(data.len))
      }
      else{
        thresh = sqrt(bc.var.mat[j,])*log(data.len)
      }

      temp1 = dj

      temp1[abs(dj)<thresh] = 0

      if (thresh.type =="soft"){
        temp2 = dj[abs(dj)>=thresh]
        temp1[abs(dj)>=thresh] = sign(temp2)*(abs(temp2)-thresh[abs(dj)>=thresh])
      }

      data.wd.thresh = putD(data.wd.thresh,level = J-j, temp1)
    }

    #invert the thresholded object to get trend estimate:

    trend.est = AvBasis(convert(data.wd.thresh))

    if (boundary.handle ==TRUE){
      if(dyadic==TRUE){
        lower = 2^(J-2)+2^(J-3)+1
        upper = 2^(J-1)+2^(J-3)
      } else{
        lower = floor((data.len-length(orig.data))/2)
        upper = lower+length(orig.data)-1
      }
      trend.est = trend.est[lower:upper]
    }

    return(trend.est)

  } else{

    #threshold the wavelet coefficients using the variance estimate and user
    #inputted rules

    for (j in 1:max.scale){

      dj = accessD(data.wd,level = J-j)

      if (normal==TRUE){
        thresh = sqrt(2*var.mat[j,]*log(data.len))
      }
      else{
        thresh = sqrt(var.mat[j,])*log(data.len)
      }

      temp1 = dj

      temp1[abs(dj)<thresh] = 0

      if (thresh.type =="soft"){
        temp2 = dj[abs(dj)>=thresh]
        temp1[abs(dj)>=thresh] = sign(temp2)*(abs(temp2)-thresh[abs(dj)>=thresh])
      }

      data.wd.thresh = putD(data.wd.thresh,level = J-j, temp1)
    }

    #invert the thresholded object to get trend estimate:

    trend.est = AvBasis(convert(data.wd.thresh))

    return(trend.est)

  }

}

