LSWsim.anydist = function (spec, distribution = "Normal")
{
  #function to simulate LSW processes using innovations from
  #other distribution families

  if (any(spec$D < 0))
    stop("All spectral elements must be non-negative")
  nlev <- nlevelsWT(spec)
  len <- 2^nlev
  for (i in (nlev - 1):0) {
    v <- accessD(spec, level = i)
    if(distribution=="Poisson"){
      v <- sqrt(v) * 2^(nlev - i) * (rpois(len,1)-1)
    }
    else if(distribution=="Exponential"){
      v <- sqrt(v) * 2^(nlev - i) * (rexp(len,rate=1)-1)
    }
    else if(distribution=="Normal"){
      v <- sqrt(v) * 2^(nlev - i) * rnorm(len, mean = 0, sd = 1)
    }
    else if(distribution=="Chisquare"){
      v <- sqrt(v) * 2^(nlev - i) * (rchisq(n=len, df = 1/2)-1/2)
    }
    spec <- putD(spec, level = i, v = v)
  }
  AvBasis(convert(spec))
}
