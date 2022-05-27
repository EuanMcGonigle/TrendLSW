wav.trend.est = function(data, filter.number = 4, family = "DaubLeAsymm",
                         max.scale = floor(log2(length(data))*0.7), type = "dec", boundary.handle =FALSE){

  #this function carries out wavelet thresholding of a time series to obtain a
  #trend estimate. All non-boundary wavelet coefficients up to a specified scale
  #are set to zero.

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
  if (boundary.handle == TRUE){
    orig.data = data
    data = get.boundary.timeseries(data)
  }
  data.len = length(data)
  J = wavethresh::IsPowerOfTwo(data.len)


  #below code determines the boundary coefficients for a given wavelet

  boundary.test = c(rep(0,data.len-1),1)

  if (type == "dec"){
    y.wd = wavethresh::wd(boundary.test,family = family, filter.number = filter.number)
  }
  else if (type == "nondec"){
    y.wd = wavethresh::wd(boundary.test,family = family, filter.number = filter.number,type = "station")
  }

  boundary.coefs = list()

  for (j in (J-1):(J-max.scale)){

    temp = wavethresh::accessD(y.wd,level = j)

    boundary.coefs[[j]] = (which(temp!=0))

  }

  #perform DWT on series
  if (type =="dec"){
    data.wd = wavethresh::wd(data,filter.number= filter.number, family = family)
  }
  else if (type =="nondec"){
    data.wd = wavethresh::wd(data,filter.number= filter.number, family = family, type = "station")
  }

  data.thresh = data.wd

  #set to zero the non-boundary coefficients

  for (j in (J-1):(J-max.scale)){

    temp = wavethresh::accessD(data.wd,level = j)

    temp[-boundary.coefs[[j]]] = 0

    data.thresh = wavethresh::putD(data.thresh, temp,level = j)

  }

  #perform inverse transform on thresholded coefficients
  if (type == "dec"){
    data_wr = wavethresh::wr(data.thresh)
  }
  else if (type == "nondec"){
    data_wr = wavethresh::AvBasis(convert(data.thresh))
  }

  #subset the longer estimate to get the ture estimate
  if (boundary.handle ==TRUE){
    if(dyadic==TRUE){
      lower = 2^(J-2)+2^(J-3)+1
      upper = 2^(J-1)+2^(J-3)
    } else{
      lower = floor((data.len-length(orig.data))/2)
      upper = lower+length(orig.data)-1
    }
    data_wr = data_wr[lower:upper]
  }

  data_wr

}
