calc.final.spec = function(spec, dyadic=TRUE,data.len){

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
