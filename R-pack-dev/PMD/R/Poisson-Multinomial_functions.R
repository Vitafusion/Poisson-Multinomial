#' @title Probability Mass Function of Poisson-Multinomial Distributions
#' 
#' @description Probability mass function of Poisson-Multinomial distributions 
#' specified by input matrix and computed through selected method. This function 
#' is capable for computation of the whole probability mass function as well as
#' of one single probability mass point. 
#' 
#' @param pmat       A matrix of probabilities. Each row of pmat should add up 
#'                   to 1.
#' @param method     Character string stands for the method selected by user to 
#'                   compute the probability mass. The method can only be one of 
#'                   the following four: 
#'                   \code{"DFT-CF"},
#'                   \code{"NA"},
#'                   \code{"SIM"},
#'                   \code{"SIM-ALL"}.
#' @param x          Result vector(probability mass point) specified by user when 
#'                   the selected method is "SIM" or "NA". The vector 
#'                   \eqn{x = (x_{1}, x_{2}, \ldots)} is used for computing 
#'                   \eqn{P(X_{1}=x_{1}, X_{2}=x_{2}, \ldots)}.
#' @param B          Simulation repeating time. Will be ignored if users do not
#'                   choose \code{"SIM-ALL"} or \code{"SIM"}
#'                   as method.
#'                   
#' @details 
#' For the methods we applied in \code{dpmd}, \code{"DFT-CF"} is an exact method 
#' to calculate all probability mass points of Poisson-Multinomial Distributions
#' via FFT algorithm. When users select \code{"DFT-CF"}, \code{dpmd} will ignore
#' \code{vec} and output the whole probability mass function.
#' 
#' \code{"SIM-ALL"} is a simulation method using naive simulation scheme to 
#' calculate the whole probability mass function, under this selection the input 
#' of \code{vec} will be ignore. Notice the accuracy and running time will be
#' effected by user choice of \code{B}. Usually \code{B}=10^5 or 10^6 will be 
#' accurate enough. Increasing \code{B} to larger than 10^8 will heavily aggravate 
#' computation burden of a CPU or GPU. 
#' 
#' Given \code{pmat} with dimension \eqn{n \times m}, the number of total probability 
#' mass points is \eqn{(n+1)^{m-1}}. Thus when the dimension of \code{pmat} 
#' increases and the users selected method is one of \code{"DFT-CF"} and 
#' \code{"SIM-ALL"}, the computation burden of \code{dpmd} might challenge the 
#' capability of a computer because both of the methods calculate all probability 
#' mass points of Poisson-Multinomial distributions.
#' 
#'  
#' \code{"SIM"} is as same as \code{"SIM-ALL"} except that it only computes a 
#' single probability mass point specified by \code{vec}.
#' 
#' \code{"NA"} specifies an approximation method using Normal approximation to 
#' compute the probability mass point of the \code{vec} vector
#' input by user.
#' 
#' @return           
#' For a single probability mass point, \code{dpmd} returns a probability value. 
#' 
#' For all probability mass points of a given \code{pmat}, it returns a 
#' multi-dimensional array. For instance, for the \code{pmat} matrix in the 
#' following example, the value of the array element \eqn{a_{1,2,1}} = 0.90 means 
#' the value of probability mass point (0,1,0,2) is 0.90. 
#'                    
#' @examples
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' 
#' 
#' dpmd(pmat = pp)
#' dpmd(pmat = pp, method = "SIM-ALL", B = 10^3)
#' dpmd(pmat = pp, x = c(0,0,1,2), method = "NA" )
#' dpmd(pmat = pp, x = c(0,0,1,2), method = "SIM", B = 10^3)
#' 
#' @export
dpmd <-function(pmat, x = c(0,0,0,0), method="DFT-CF", B=10^3)
{
  if(is.matrix(pmat)==F){
    stop("pmat is not a matrix.")
  }
  if(any(pmat<0)|any(pmat>1))
  {
    stop("Invalid values in pmat.")
  }
  for(i in 1:nrow(pmat)){
    if(abs(sum(pmat[i,])-1)>10^(-16))
      stop("Existing a row that doesn't sum up to 1.")
  }
  
  
  switch(method,
         "DFT-CF"={
           mm=ncol(pmat) # ncol of pmat
           nn=nrow(pmat) # nrow of pmat
           
           nn.vec=rep(nn+1, mm-1)
           l.vec=rep(0, mm-1)
           cn.vec=cumprod(nn.vec)
           cn.vec=c(1, cn.vec[-(mm-1)])
           cn.vec=cn.vec[length(cn.vec):1]
           cn.vec=as.integer(cn.vec) #((n+1)^(m-2),...,(n+1)^2,(n+1),1)
           
           nnt=prod(nn.vec) # (n+1)^(m-1) probability mass points
           
           #browser()
           
           res0=pmn_mdfft_arma(nnt, pmat, nn.vec, l.vec, cn.vec)
           
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
         "SIM-ALL"={
             mm=ncol(pmat) # ncol of pmat
             nn=nrow(pmat) # nrow of pmat
             nn.vec=rep(nn+1, mm-1)
             l.vec=rep(0, mm-1)
             cn.vec=cumprod(nn.vec)
             cn.vec=c(1, cn.vec[-(mm-1)])
             cn.vec=cn.vec[length(cn.vec):1]
             cn.vec=as.integer(cn.vec) #((n+1)^(m-2),...,(n+1)^2,(n+1),1)
             nnt=prod(nn.vec) # (n+1)^(m-1) probability mass points
             
             res0 = pmd_simulation_allpoints(pmat, nnt, l.vec, cn.vec, B)
             
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
         "NA"=   {
           mm=ncol(pmat) # m categories
           nn=nrow(pmat) # n people
           if(sum(x)>nn|any(x<0)|length(x)!=mm)
           {
             stop("Invalid value or length of x.")
           }
           mm = mm - 1
           x_vec = x[1:(length(x)-1)]
           lb = as.numeric(x_vec - 0.5)
           ub = as.numeric(x_vec+0.5)
           res = 0
           
           # asymptotic sigma
           n = nn
           m = mm
           if(n==1) pmat = t(as.matrix(pmat[,1:m])) else pmat = as.matrix(pmat[,1:m])
           sig = matrix(0,m,m)
           for (i in 1:n) {
             sig = sig + diag(pmat[i,],nrow = m) - pmat[i,]%*%t(pmat[i,])
           }
           
           # asymptotic mu
           mu = matrix(0, nrow = 1, ncol = ncol(pmat))
           for (i in 1:n) {
             mu = mu + pmat[i,]
           }
           mu = as.vector(mu)
           res = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = n*sig)
           res = res[[1]]
         },
         "SIM" = {
           mm=ncol(pmat) # m categories
           nn=nrow(pmat) # n people
           if(sum(x)>nn|any(x<0)|length(x)!=mm)
           {
             stop("Invalid value or length of x.")
           }
             res = pmd.by.demands(x,pmat,B)
         }
         
  )
  
  return(res)
}



########################################################################################
#' @title Cumulative Distribution Function of Poisson-Multinomial Distribution
#' 
#' @description This function computes cumulative distribution function of 
#' Poisson-Multinomial distributions that specified by input probability matrix 
#' via given method.
#'  
#' @param pmat       A matrix of probabilities. Each row of pmat should add up 
#'                   to 1.
#' @param method     Character string stands for the method selected by user to 
#'                   compute the probability mass. The method can only be one of 
#'                   the following three: 
#'                   \code{"DFT-CF"},
#'                   \code{"NA"},
#'                   \code{"SIM-ALL"}. 
#' @param B          Simulation repeating time. Will be ignored if users do not
#'                   choose \code{"SIM-ALL"} as method.
#' @param x          Vector \eqn{x = (x_{1},x_{2},\ldots)} for computing 
#'                   \eqn{P(X_{1} \leq x_{1},X_{2} \leq x_{2},\ldots)}.
#' 
#' @details 
#' Three methods are same as listed in the details of \code{dpmd} but \code{"NA"}
#' stands for normal approximation. 
#' 
#' \code{ppmd} computes the cumulative distribution function by adding all probability 
#' mass points within hyper-dimensional space limited by \code{x}. 
#' 
#' @return 
#' The value of \eqn{P(X_{1} \leq x_{1},X_{2} \leq x_{2},\ldots)} of given 
#' \eqn{x = (x_{1},x_{2},\ldots)}.
#' 
#' @examples
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' 
#' ppmd(pmat = pp, x = c(3,2,1,3))
#' ppmd(pmat = pp, x = c(3,2,1,3), method = "NA")
#' ppmd(pmat = pp, x = c(3,2,1,3), method = "SIM-ALL", B = 10^3)
#' @export
ppmd = function(pmat,x,method="DFT-CF",B=10^3){
  if(is.matrix(pmat)==F){
    stop("pmat is not a matrix.")
  }
  if(any(pmat<0)|any(pmat>1)){
    stop("Invalid values in pmat.")
  }
  nn = nrow(pmat)
  mm = ncol(pmat)
  if(any(x<0)|length(x)!=mm)
  {
    stop("Invalid value or length of x.")
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
  
  # filter probability mass points
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
           res = dpmd(pmat)
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
         "SIM-ALL" = {
             T=B
             points.pos = points[which(points[,mm]>=0),]
             prob = 0
             for(i in 1:nrow(points.pos)){
               prob = prob + pmd.by.demands(as.numeric(points.pos[i,]),pmat,T)
           }
         },
         "NA" = {
           prob = 0
           points.pos = points[which(points[,mm]>=0),]
           for(i in 1:nrow(points.pos)){
             prob = prob + dpmd(pmat, x = points.pos[i,], method="NA")
           }
         })
  return(prob)
}
########################################################################################
#' @title Poisson-Multinomial Distribution Random Number Generator
#' @description Generating random samples of given a Poisson-Multinomial distribution.
#'  
#' @param pmat       A matrix of probabilities. Each row of pmat should add up 
#'                   to 1.
#' @param n          Number of samples to be generated.
#' 
#' @return 
#' A matrix of samples, each row stands for one sample.
#' 
#' @examples 
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' 
#' 
#' rpmd(pmat = pp, n = 5)
#' 
#' @export
rpmd = function(pmat, n=1){
  if(is.matrix(pmat)==F){
    stop("pmat is not a matrix.")
  }
  if(any(pmat<0)|any(pmat>1)){
    stop("Invalid values in pmat.")
  }
  mm = ncol(pmat)
  rnd = matrix(NA,nrow = n,ncol = mm)
  for(i in 1:n){
    rnd[i,] = t(rpmd_arma(pmat))
  }
  return(rnd)
}



