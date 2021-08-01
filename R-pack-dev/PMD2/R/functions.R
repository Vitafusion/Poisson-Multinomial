
################################################################################
#' l.vec
#'
#' @param k input by given functions
#' @param cn.vec calculated by given fucntions
#' @param m dimensions
#'
#' @return A vector, which will be used in other functions


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
#' probability density of Poisson multinomial distribution
#'
#' @param pp The pp is a probability matrix which will be input by user
#' @param vec result vec input by user
#' @param method method selected by user to compute the probability mass
#' @param t simulation repeat time
#' @return The probability mass of PMD
#' @examples
#' aa=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2, .5, .1, .1, .3), nrow=4, byrow=TRUE)
#' pp=aa[1:3,]
#' dpmn(pp)
#' @export
#'
dpmn <-function(pp,method="DFT-CF",vec=c(0,0,0,0,0),t=100)
{
  
  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }
  
  
  switch(method,
         "DFT-CF"={
           mm=ncol(pp) # ncol of pp
           nn=nrow(pp) # nrow of pp
           
           nn.vec=rep(nn+1, mm-1)
           l.vec=rep(0, mm-1)
           cn.vec=cumprod(nn.vec)
           cn.vec=c(1, cn.vec[-(mm-1)])
           cn.vec=cn.vec[length(cn.vec):1]
           cn.vec=as.integer(cn.vec) #((n+1)^(m-2),...,(n+1)^2,(n+1),1)
           
           nnt=prod(nn.vec) # (n+1)^(m-1) density points
           
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
           
           # asymptotic sigma
           n = nn
           m = mm-1
           if(n==1) pp = t(as.matrix(pp[,1:m])) else pp = as.matrix(pp[,1:m])
           sig = matrix(0,m,m)
           for (i in 1:n) {
             sig = sig + diag(pp[i,],nrow = m) - pp[i,]%*%t(pp[i,])
           }
           
           # asymptotic mu
           mu = matrix(0, nrow = 1, ncol = ncol(pp))
           for (i in 1:n) {
             mu = mu + pp[i,]
           }
           mu = as.vector(mu)
           res = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = sig)
         }
         
  )
  
  return(res)
}

################################################################################
#' pmatrix
#'
#' @param n column dimension
#' @param m row dimension
#'
#' @return a randomly generated Poisson multinomial distribution probability matrix
#' @export
#'
#' @examples
#' pp = pmatrix(2,2)
#' pp
#' @export
#'
pmatrix = function(n,m){
  pp = matrix(0,nrow = n,ncol = m,byrow = T)
  for (i in 1:n) {
    r = runif(m)
    r = r/sum(r) #generate row
    pp[i,] = r
  }
  return(pp)
}


########################################################################################
#' cumulative distribution function of PMN
#'
#' @param pp input matrix of probabilities
#' @param x input vector
#' @param method method
#' @param t repeating time
#' @return prob
#' @examples
#' aa=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2, .5, .1, .1, .3), nrow=4, byrow=TRUE)
#' pp=aa[1:3,]
#' ppmn(pp,c(3,2,1,3))
#' @export
#'
ppmn = function(pp,x,method="DFT-CF",t=1000){
  if(any(pp<0)|any(pp>1)){
    stop("invalid values in pp.")
  }
  nn = nrow(pp)
  mm = ncol(pp)
  if(any(x<0)){
    stop("invalid values in x.")
  }
  if(length(x)!=mm){
    stop("invalid format of x.")
  }
  #idx formed
  nn.vec=rep(nn+1, mm-1)
  l.vec=rep(0, mm-1)
  cn.vec=cumprod(nn.vec)
  cn.vec=c(1, cn.vec[-(mm-1)])
  cn.vec=cn.vec[length(cn.vec):1]
  cn.vec=as.integer(cn.vec)
  
  nnt=prod(nn.vec)
  idx = as.data.frame(matrix(0,nrow = nnt,ncol=mm))
  idx0 = as.data.frame(matrix(0,nrow = nnt,ncol=(mm-1)))
  #transfer l.vec to result vector(idx)
  for(i in 1:nnt)
  {
    idx[i,1:(mm-1)]=l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
    idx0[i,] = l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
    idx[i,1:(mm-1)] = idx[i,1:(mm-1)]-1
    idx[i,mm] = nn - sum(idx[i,1:(mm-1)])
    #print(idx)
  }
  
  # filter density points
  conditions = c()
  expr0 = 'which('
  expr1 = ')'
  for(i in 1:mm){
    conditions[i] = paste0('idx$V',i,'<=',x[i])
  }
  cond = conditions[1]
  for(i in 1:(mm-1)){
  cond = paste0(cond,'&',conditions[i+1])
  }
  expr = paste0(expr0,cond,expr1)
  index = eval(parse(text=expr))
  points = idx[index,]
  switch(method,
         "DFT-CF" = {
           res = dpmn(pp)
           temp.index = idx0[index,]
           prob  = 0
           res.expr="prob = prob + res[temp[1]"
           if(mm>=3)
           {
             for(i in 2:(mm-1))
             {
               res.expr=paste0(res.expr, ", temp[", i, "]")
             }
           }
           res.expr=paste0(res.expr, "]")
           
           for(i in 1:nrow(temp.index)){
            temp =  as.numeric(temp.index[i,])
            eval(parse(text=res.expr))
           }
           
         },
         "simulation" = {
             T=t
             points.pos = points[which(points[,mm]>=0),]
             prob = 0
             for(i in 1:nrow(points.pos)){
               prob = prob + pmd.by.demands(as.numeric(points.pos[i,]),pp,T)
           }
         },
         "NA" = {
           prob = 0
           for(i in 1:nrow(points.pos)){
             prob = prob + dpmn(pp,method="NA",vec = points.pos[i,])
           }
         })
  return(prob)
}
########################################################################################
#' generate random number from PMD
#'
#' @param pp input matrix of probabilities
#' @return the random number vector generated from PMD
#' @examples
#' aa=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2, .5, .1, .1, .3), nrow=4, byrow=TRUE)
#' rpmd(pp)
#' @export
#'
rpmd = function(pp){
  if(any(pp<0)|any(pp>1)){
    stop("invalid values in pp.")
  }
  
  rnd = rpmd_arma(pp)
  
  return(rnd)
}

###############################################################################
#' PMN density calculated by input vector
#'
#' @param pp The pp is a probability matrix which will be input by user
#' @param x_vec result vec input by user
#' @param t simulation repeat time
#' @return The probability mass of PMD
#' @examples
#' aa=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2, .5, .1, .1, .3), nrow=4, byrow=TRUE)
#' pp=aa[1:3,]
#' pmd.by.demands(c(1,2,0,0),pp,t=10^5)
#' @export
#'
pmd.by.demands = function(x_vec,pp,t=1000){
    x_vec = as.vector(x_vec)
    nn = nrow(pp)
    mm = ncol(pp)
    # if(sum(x_vec)!=nn)   stop("invalid x_vec.")
    res0=0
    #input simulation method here
    temp=pm_simulation_arma(pp, x_vec, t)
    #temp = .C("pmd_simulation_vec",as.double(res0), as.integer(nn), as.integer(mm), as.double(pp), as.integer(x_vec), as.double(t), PACKAGE = "poissonmulti")
    res0=round(temp[[1]],10)
    return(res0)
}

