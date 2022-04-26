easy.spec.plot = function(spectrum, bylev = FALSE){

  scales = nlevelsWT(spectrum)
  if (bylev ==FALSE){
    plot(spectrum, ylabchars = -(1:scales), main = "", sub = "", ylab = "Scale")
  }
  else{
    plot(spectrum, ylabchars = -(1:scales), main = "", sub = "", ylab = "Scale", scaling = "by.level")
  }
}
