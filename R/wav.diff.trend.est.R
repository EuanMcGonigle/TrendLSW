wav.diff.trend.est = function(data, spec, filter.number = 4, thresh.type = "soft",
                              normal = TRUE,family = "DaubLeAsymm",
                              max.scale = floor(0.7*log2(length(data))),
                              trans.type = "nondec", boundary.handle = FALSE){

  data.len = length(data)
  J = log2(data.len)

  #by default, do T.I. denoising:

  if (boundary.handle == TRUE){
    orig.data = data
    data = get.boundary.timeseries(data)
    max.scale = max.scale+2
  }

  if (max.scale <1){
    max.scale = 1
  }
  if (max.scale >J){
    max.scale = J-1
  }

  #below code determines the boundary coefficients for a given wavelet

  boundary.test = c(rep(0,data.len-1),1)

  if (trans.type == "dec"){
    y.wd = wavethresh::wd(boundary.test,family = family, filter.number = filter.number)
  }
  else if (trans.type == "nondec"){
    y.wd = wavethresh::wd(boundary.test,family = family, filter.number = filter.number,type = "station")
  }

  boundary.coefs = list()

  for (i in (J-1):(J-max.scale)){

    temp = accessD(y.wd,level = i)

    boundary.coefs[[i]] = (which(temp!=0))

  }

  if (trans.type=="nondec"){
    data.wd = wavethresh::wd(data, filter.number = filter.number, family = family,
                             type = "station")
  }
  else{
    data.wd = wavethresh::wd(data, filter.number = filter.number, family = family)
  }

  #calculate C to use in variance estimate of wavelet coefficients

  C = Cmat.calc(J=max.scale, an.filter.number = filter.number,
                gen.filter.number = spec$filter$filter.number,
                an.family = family, gen.family = spec$filter$family)


  #create EWS estimate matrix

  spec.mat = matrix(0,nrow = max.scale,ncol = data.len)

  for (j in 1:max.scale){
    spec.mat[j,] = accessD(spec,level = J-j)
  }

  #calulate variance estimate matrix

  var.mat = C%*%spec.mat

  #find negative values and replace them with nearby positive values:

  for (j in 1:max.scale){

      var.row = var.mat[j,]

      if(sum(var.row<0)>0){

        var.row[var.row<0] = 0
        var.row0 = which(var.row==0)
        var.row.non0 = var.row[which(var.row!=0)]

        for(i in 1:length(var.row0)){
          var.mat[j,(var.row0[i])]= var.row.non0[which.min(abs(var.row0[i]-which(var.row!=0)))]
          }
        }
  }

  if (boundary.handle==TRUE){

    rev.var.mat = var.mat[,ncol(var.mat):1]

    spec.mat = cbind(rev.var.mat[,1:(data.len/2)],rev.var.mat,var.mat,rev.var.mat,rev.var.mat[,1:(data.len/2)])
    J = J+2

  }

  data.wd.thresh = data.wd

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

