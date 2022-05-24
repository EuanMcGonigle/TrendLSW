get.boundary.timeseries = function(data){

  #this function takes a time series and produces a 4*data.len length series that adds
  #boundary handling to both sides.

  data.len = length(data)
  J = wavethresh::IsPowerOfTwo(data.len)

  s = seq(from = 0 , to = (data.len-1)/data.len,length = data.len)

  l = lm(data~poly(s,3,raw = TRUE))

  bh_right = 2*predict(l, newdata = data.frame(s=1))

  bh_left = 2*predict(l, newdata = data.frame(s=-1/data.len))

  bh_series1 = c(-rev(data)+bh_left, data, -rev(data)+bh_right)
  if(is.na(J)==TRUE){
    l = 2^floor(log2(length(bh_series1)))
    k = floor((3*data.len-l)/2)
    if(data.len%%2==0){
      bh_series2 = bh_series1[(k+1):(3*data.len-k)]
    } else{
      bh_series2 = bh_series1[(k+1):(3*data.len-k-1)]
    }

    return(bh_series2)
  } else{
    l = length(bh_series1)

    bh_series2 = c((data)[(data.len/2+1):data.len]-abs(bh_right-bh_left),bh_series1, (data)[1:(data.len/2)]+abs(bh_right-bh_left))

    return(bh_series2)
  }



}
