spec.to.mat = function(spec){

  J = nlevelsWT(spec)

  spec.mat = matrix(0,ncol = 2^J, nrow = J)

  for (j in 1:J){

    spec.mat[j,] = accessD(spec,level = J-j)

  }

  spec.mat


}

mat.to.spec = function(s.mat){

  J = nrow(s.mat)

  spec = cns(2^J)

  for (j in 1:J){

    spec = putD(spec, level = J-j, s.mat[j,])

  }

  spec


}
