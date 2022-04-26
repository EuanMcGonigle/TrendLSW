per.var.est = function(data,spec, filter.number = 1, family = "DaubExPhase",
         max.scale = floor(0.7*log2(length(data)))){
  
  #function computes the variance of the raw wavelet periodogram, given an estimate of 
  #wavelet spectrum. 
  
  data.len = length(data)
  J = log2(data.len)
  
  
  if (max.scale <1){
    max.scale = 1
  }
  if (max.scale >J){
    max.scale = J-1
  }
  
  #bcalculate wavelet transform:
  
  data.wd = wavethresh::wd(data, filter.number = filter.number, family = family,
                             type = "station")

  
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
  
  
  return(var.mat)
  
  
}
