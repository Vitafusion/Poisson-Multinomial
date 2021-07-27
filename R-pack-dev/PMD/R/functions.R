################################################################################
#'mu.cmpt
#'
#' @param p The p is a probability vmatrix which will be input by user
#' @return The probability mass of mean mu
#' @export
#'

mu.cmpt = function(p){
  n = nrow(p)
  m = ncol(p)-1
  if(n==1)  pp = t(as.matrix(p[,1:m])) else pp = as.matrix(p[,1:m])
  mu = matrix(0, nrow = 1, ncol = ncol(pp))
  for (i in 1:n) {
    mu = mu + pp[i,]
  }
  return(as.vector(mu))
}
################################################################################
#'sig.cmpt
#'
#' @param p The p is a probability vmatrix which will be input by user
#' @return The probability mass of variance sigma
#' @export
#'
sig.cmpt = function(p){
  n = nrow(p)
  m = ncol(p)-1
  if(n==1) pp = t(as.matrix(p[,1:m])) else pp = as.matrix(p[,1:m])
  sig = matrix(0,m,m)
  for (i in 1:n) {
    sig = sig + diag(pp[i,],nrow = m) - pp[i,]%*%t(pp[i,])
  }
  return(sig)
}
################################################################################
#' Title
#'
#' @param k input by given functions
#' @param cn.vec calculated by given fucntions
#' @param m dimensions
#'
#' @return A vector, which will be used in other functions
#' @export

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
#'pmd
#'
#' @param pp The pp is a probability matrix which will be input by user
#' @param vec result vec input by user
#' @param method method selected by user to compute the probability mass
#' @param t simulation repeat time
#' @return The probability mass of PMD
#' @export
#'
pmd <-function(pp,method="DFT-CF",vec=c(0,0,0,0,0),t=100)
{
  
  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }
  
  
  switch(method,
         "DFT-CF"={
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
           
         },
        # "simulation"=    {
           #simulation method
        # },
         "NA"=   {
           mm=ncol(pp) # m categories
           nn=nrow(pp) # n people
           if(sum(vec)>nn|any(vec<0)|length(vec)!=mm)
           {
             stop("invalid result vector")
           }
           mm = mm - 1
           x_vec = vec[1:(length(vec)-1)]
           lb = x_vec - 0.5
           ub = x_vec+0.5
           res = 0
           sig = sigma.calcu(pp)
           mu = mu.calcu(pp)
           res = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = sig)
         }
         
  )
  
  return(res)
}

################################################################################
#' Title
#'
#' @param n column dimension
#' @param m row dimension
#'
#' @return a randomly generated Poisson multinomial distribution probability matrix
#' @export
#'
#' @examples
#' pmatrix(2,2)
#' 
pmatrix <- function(n,m){
  p <- matrix(0,nrow = n,ncol = m,byrow = T)
  for (i in 1:n) {
    r <- runif(m)
    r <- r/sum(r) #generate row
    p[i,] <- r
  }
  return(p)
}