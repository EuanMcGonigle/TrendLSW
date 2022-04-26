macg.vs.wav.est = function(data, filter.number = 4, family = "DaubLeAsymm", max.scale,
                           boundary.handle = FALSE){

  data.len = length(data)
  J = log2(data.len)

  data.wd = wd(data,filter.number = filter.number, family = family)
  data.t = data.wd

  for (j in 1:max.scale){

    m = floor((J-j-1)/2)
    l = 2^(m+1)
    #print(l)

    thresh = rep(0,2^(J-j))
    #print(length(thresh))

    temp = accessD(data.wd,level = J-j)
    #print(length(temp))

    len = length(temp)/(2^(m+1))
    #print(len)

    for(i in 1:len){
      thresh[(l*(i-1)+1):(l*i)] = mad(temp[(l*(i-1)+1):(l*i)])*sqrt(log(2^(J-j)))
    }
    #plot.ts(thresh)

    temp[abs(temp)<thresh] = 0

    data.t = putD(data.t,level = J-j,temp)

  }


  if (boundary.handle==TRUE){

    data2.wd = wd(data2, filter.number = filter.number,family = family)
    data2 = get.boundary.timeseries(data)

    for (j in 1:max.scale){

      n1 =data.len*2^(-j)

      n2 = (data.len*4*2^(-j)-n1)/2

      print(n1)
      print(n2)

      temp = accessD(data2.wd,level = J+2-j)

      temp2 = accessD(data.t,level = J-j)

      print(length(temp2))

      temp3 = c(temp[1:n2],temp2,temp[(4*length(temp2)-n2+1):(4*length(temp2))])

      print(length(temp3))

      plot.ts(temp3)

      data2.wd = putD(data2.wd,level = J+2-j, temp3)

    }

    plot(data2.wd)

    data.est = wr(data2.wd)

    lower = 2^(J)+2^(J-1)+1
    upper = 2^(J+1)+2^(J-1)
    data.est = data.est[lower:upper]
  }
  else{
    data.est = wr(data.t)
  }


  data.est

}


