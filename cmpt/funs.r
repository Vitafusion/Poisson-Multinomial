#library(Rcpp)
#library(RcppArmadillo)
#options(stringsAsFactors = FALSE)
#sourceCpp("pm-fft.cpp")

################################################################################
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
################################################################################
dpmn.arma=function(kk=NULL, pp)
{
  mm=ncol(pp)
  nn=nrow(pp)
    
  nn.vec=rep(nn+1, mm-1)
  l.vec=rep(0, mm-1)
  cn.vec=cumprod(nn.vec)
  cn.vec=c(1, cn.vec[-(mm-1)])
  cn.vec=cn.vec[length(cn.vec):1]
  cn.vec=as.integer(cn.vec)
  
  nnt=prod(nn.vec)
    
  #browser()
  
  res0=pmn_mdfft_arma(nnt, pp, nn.vec, l.vec, cn.vec)
  
  #example an_array[k + 27 * (j + 12 * i)]
  #print(round(res0, 9))

  res=array(0, nn.vec)
  
  res.expr="res[idx[1]"
  if(mm>=3)
  {
    for(i in 2:(mm-1))
    {
       res.expr=paste0(res.expr, ", idx[", i, "]")
    }
  }
  res.expr=paste0(res.expr, "]=res0[i]")
  
  #browser()
  
  #print(nnt)
  
  for(i in 1:nnt)
  {
    idx=l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
    #print(idx)
    eval(parse(text=res.expr))
  }

  res=round(res, 10)
  return(res)
}
################################################################################

