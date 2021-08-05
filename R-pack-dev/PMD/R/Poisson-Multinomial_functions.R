
################################################################################
#' @title Probability Mass of Poisson-Multinomial Distributions
#'
#' @description
#' Computation of probability mass for Poisson-Multinomial Distributions using exact, simulation, approximation methods. Users are allowed to specified a method and can choose to compute single mass point or all mass points. For simulation method, users can also choose the repeating time to enhance the accuracy of outcomes.
#' 
#' @param pp         A matrix of probabilities.
#'                   Each row of pp should add up to 1.
#' @param vec        Result vector(probability mass point) specified by user.
#'                   eg. pp is 4 by 3 matrix then a user might be interested in the probability of getting result: vec=c(0,0,1,2).
#' @param method     Method selected by user to compute the probability mass. There are totally 4 methods.
#'                  \code{"DFT-CF"}: An exact method to calculate all probability mass points of Poisson-Multinomial Distributions via FFT algorithm.
#'                   \code{"simulation"}: A simulation method calculating all probability mass points.
#'                   \code{"NA by demands"}: An approximation method using Normal approximation to compute the probability for the 'vec' vector input by user.
#'                   \code{"simulation by demands"}: The same simulation method as above just to compute single probability mass point as input by user.
#' @param B Simulation repeat time.
#' 
#' @return For a single mass point, \code{dpmd} returns a probability.
#'         For all probability mass points of a given \code{pp}, it returns a multi-dimensional array. For instance, for the \code{pp} matrix in the following example, the value of the array element \eqn{a_{1,2,1}} = 0.90 means the value of probability mass point (0,1,0,2) is 0.90.
#'    
#' @examples
#' 
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow=3, byrow=TRUE)
#' dpmd(pp)
#' dpmd(pp,"simulation",B=10^3)
#' dpmd(pp,"NA by demands", vec = c(0,0,1,2))
#' dpmd(pp,"simulation by demands", vec = c(0,0,1,2), B=10^3)
#' @export
#'
dpmd <-function(pp,method="DFT-CF",vec=c(0,0,0,0,0),B=100)
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
         "simulation"={
             mm=ncol(pp) # ncol of pp
             nn=nrow(pp) # nrow of pp
             nn.vec=rep(nn+1, mm-1)
             l.vec=rep(0, mm-1)
             cn.vec=cumprod(nn.vec)
             cn.vec=c(1, cn.vec[-(mm-1)])
             cn.vec=cn.vec[length(cn.vec):1]
             cn.vec=as.integer(cn.vec) #((n+1)^(m-2),...,(n+1)^2,(n+1),1)
             nnt=prod(nn.vec) # (n+1)^(m-1) density points
             
             res0 = pmd_simulation_allpoints(pp, nnt, l.vec, cn.vec, B)
             
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
         "NA by demands"=   {
           mm=ncol(pp) # m categories
           nn=nrow(pp) # n people
           if(sum(vec)>nn|any(vec<0)|length(vec)!=mm)
           {
             stop("invalid result vector")
           }
           mm = mm - 1
           x_vec = vec[1:(length(vec)-1)]
           lb = as.numeric(x_vec - 0.5)
           ub = as.numeric(x_vec+0.5)
           res = 0
           
           # asymptotic sigma
           n = nn
           m = mm
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
           res = res[[1]]
         },
         "simulation by demands" = {
             res = pmd.by.demands(vec,pp,B)
         }
         
  )
  
  return(res)
}



########################################################################################
#' @title cumulative mass function of PMN
#'
#' @description  By an input vector x = \eqn{(x_{1},x_{2},...)}, this function compute \eqn{P(X_{1} \leq x_{1}, X_{2} \leq x_{2}, ...)}. 
#' @param pp input matrix of probabilities
#' @param x input result vector
#' @param method method selected by users to compute the cumulative mass probabilities.
#' @param B repeating time
#' @return prob
#' @examples
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow=3, byrow=TRUE)
#' ppmd(pp,c(3,2,1,3))
#' ppmd(pp,c(3,2,1,3),"simulation",B=10^3)
#' ppmd(pp,c(3,2,1,3),"NA")
#' @export
#'
ppmd = function(pp,x,method="DFT-CF",B=1000){
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
           res = dpmd(pp)
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
             T=B
             points.pos = points[which(points[,mm]>=0),]
             prob = 0
             for(i in 1:nrow(points.pos)){
               prob = prob + pmd.by.demands(as.numeric(points.pos[i,]),pp,T)
           }
         },
         "NA" = {
           prob = 0
           points.pos = points[which(points[,mm]>=0),]
           for(i in 1:nrow(points.pos)){
             prob = prob + dpmd(pp,method="NA by demands",vec = points.pos[i,])
           }
         })
  return(prob)
}
########################################################################################
#' generate random number from PMD
#'
#' @param pp input matrix of probabilities
#' @return the random number vector generated from PMD.
#' @examples
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow=3, byrow=TRUE)
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



