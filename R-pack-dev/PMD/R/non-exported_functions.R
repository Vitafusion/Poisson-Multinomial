#' @noRd
l.vec.compute=function(k, cn.vec, m)
{
  k=k-1
  l.vec=rep(0, m-1)
  for(i in 1:(m-1))
  {
    aa=k%%cn.vec[i]
    bb=(k-aa)/cn.vec[i]
    l.vec[i]=bb
    k=aa
  }
  l.vec=l.vec+1
  return(l.vec)
}



#' @noRd
pmd.by.demands = function(x_vec,pp,B=1000){
    x_vec = as.vector(x_vec)
    nn = nrow(pp)
    mm = ncol(pp)
    # if(sum(x_vec)!=nn)   stop("invalid x_vec.")
    res0=0
    #input simulation method here
    temp=pmd_simulation_singlepoint(pp, x_vec, B)
    res0=round(temp[[1]],10)
    return(res0)
}
