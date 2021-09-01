#' @title Probability Mass Function of Poisson-Multinomial Distributions
#' 
#' @description Probability mass function of Poisson-Multinomial distributions 
#' specified by input matrix and computed through selected method. This function 
#' is capable for computation of the whole probability mass function as well as
#' of one single probability mass point. 
#' 
#' @param pmat      An \eqn{n \times m} matrix of probabilities. \eqn{n} is the number of independent trials.
#'                  \eqn{m} is the number of categories. Also called success probability matrix.
#'                  Each row of pmat describes the success probability for the corresponding
#'                  trial and it should add up to 1.
#' @param method     Character string stands for the method selected by user to 
#'                   compute the probability mass. The method can only be one of 
#'                   the following four: 
#'                   \code{"DFT-CF"},
#'                   \code{"NA"},
#'                   \code{"SIM"},
#'                   \code{"SIM-ALL"}.
#' @param xmat       Result matrix of length \eqn{m} (probability mass point) specified by user. Each row of the matrix should has the form
#'                   \eqn{x = (x_{1}, \ldots, x_{m})} which is used for computing 
#'                   \eqn{P(X_{1}=x_{1}, \ldots, X_{m} = x_{m})}.
#' @param B          Number of repetitions in the simulation method. Will be ignored if users do not
#'                   choose \code{"SIM"} method.
#'                   
#' @details
#' Consider \eqn{n} independent trials and each trial leads to a success for exactly one of \eqn{m} categories. Each category has varying success probabilities from different trials. The Poisson multinomial distribution (PMD) gives the probability of any particular combination of numbers of successes for the \eqn{m} categories. The success probabilities form an \eqn{n \times m} matrix, which is called the success probability matrix and denoted by \code{pmat}. 
#' The total number of outcomes is \eqn{(n+1)^{m-1}}. 
#' For the methods we applied in \code{dpmd}, \code{"DFT-CF"} is an exact method 
#' to calculate all mass points of Poisson-Multinomial Distributions
#' via FFT algorithm. When users select \code{"DFT-CF"}, \code{dpmd} will ignore
#' \code{vec} and return the probability mass function for all outcomes.
#' 
#' \code{"SIM"} is a simulation method using a naive simulation scheme to 
#' calculate the whole probability mass function if the input of \eqn{xmat} is not specified. Notice that the accuracy and running time will be
#' affected by user choice of \code{B}. Usually \code{B}=1e5 or 1e6 will be 
#' accurate enough. Increasing \code{B} to larger than 1e8 will heavily aggravate 
#' computational burden of a CPU or GPU. 
#' 
#' When the dimension of \code{pmat} increases, the computation burden of \code{"DFT-CF"} as well as 
#' \code{"SIM"}  method might challenge the 
#' capability of a computer because both of the methods calculate all
#' mass points of Poisson-Multinomial distributions.
#' 
#' \code{"SIM"} is as same as \code{"SIM-ALL"} except that it only computes the
#'  probability mass function at a single outcome specified by \code{vec}.
#' 
#' \code{"NA"} is an approximation method using Normal approximation to 
#' compute the probability mass function of \code{vec} vector
#' specified by user.
#' 
#' @return           
#' For a single mass point, \code{dpmd} returns the probability mass function at that point. 
#' 
#' For all mass points of a given \code{pmat}, it returns a 
#' multi-dimensional array. For instance, for the \code{pmat} matrix in the 
#' following example, the value of the array element \eqn{a_{1,2,1}} = 0.90 means 
#' the value of probability mass point (0,1,0,2) is 0.90. 
#'                    
#' @examples
#' pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' x <- matrix(c(0,0,1,2), nrow=1) 
#' x1 <- matrix(c(0,0,1,2,2,1,0,0),nrow=2,byrow=TRUE)
#'
#' dpmd(pmat = pp)
#' dpmd(pmat = pp, x = x1)
#' dpmd(pmat = pp, x = x)
#'
#' dpmd(pmat = pp, x = x, method = "NA" )
#' dpmd(pmat = pp, x = x1, method = "NA" )
#'
#' dpmd(pmat = pp, x = x, method = "SIM", B = 1e3)
#' dpmd(pmat = pp, x = x1, method = "SIM", B = 1e3)
#' 
#' @export
dpmd <-function(pmat, x = NULL, method="DFT-CF", B=1e3)
{
  chck = pmat.check(pmat,x)
  if(chck!=1){stop(chck)}
  # x should be matrix
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

           if(!is.null(x))
           {
            nrow.x=nrow(x)
            res.x=rep(NA,nrow.x)
            if(nrow.x!=1)
            {
              for(j in 1:nrow.x)
              {
                idx.x=x[j,1:(mm-1)]+1
                res.x.expr="res.x[j]=res[idx.x[1]"
                if(mm>=3)
                {
                  for(i in 2:(mm-1))
                  {
                    res.x.expr=paste0(res.x.expr, ", idx.x[", i, "]")
                  }
                }
                res.x.expr=paste0(res.x.expr, "]")
                eval(parse(text=res.x.expr))
              }
              res=matrix(res.x,ncol=1)
            }
            else
            {
              res.x.expr="res.x=res[idx.x[1]"
              idx.x=x[1:(mm-1)]+1
                if(mm>=3)
                {
                  for(i in 2:(mm-1))
                  {
                    res.x.expr=paste0(res.x.expr, ", idx.x[", i, "]")
                  }
                }
                res.x.expr=paste0(res.x.expr, "]")
                eval(parse(text=res.x.expr))
                res=res.x
            }
           }
           
           
         },
         "SIM"={
            if(is.null(x))
            {
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
            }
            else
            {
              nrow.x = nrow(x)
              res = rep(NA,nrow.x)
              if(nrow.x!=1)
              {
                for (i in 1:nrow.x) {
                  res[i] = pmd.by.demands(x[i,],pmat,B)
                }
                res=matrix(res,ncol=1)
              }
              else
              {
                res = pmd.by.demands(x,pmat,B)
              }
            }
         },
         "NA"=   {
           mm=ncol(pmat) # m categories
           nn=nrow(pmat) # n people
           
           if(is.null(x))
           {
            stop("Value of x is not assigned.")
           }
           
           nrow.x = nrow(x)
           
           for (i in 1:nrow.x) 
           {
            if(sum(x[i,])>nn)
            {
              stop("Sum of a row of x greater than n.")
            }
           }
           
           mm = mm - 1
           
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


           res = rep(NA,nrow.x)


           if(nrow.x!=1)
           {
            for (i in 1:nrow.x) 
            {
              x_vec = x[i,1:mm]
              lb = as.numeric(x_vec - 0.5)
              ub = as.numeric(x_vec+0.5)
              res0 = 0
           
              res0 = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = sig)
              res[i] = res0[[1]]
            }
              res=matrix(res,ncol=1)
            }
            else
            {
              x_vec = x[1:mm]
              lb = as.numeric(x_vec - 0.5)
              ub = as.numeric(x_vec+0.5)
              res0 = 0
              res0 = mvtnorm::pmvnorm(lower=lb,upper = ub, mean = mu, sigma = sig)
              res = res0[[1]]
            }
         },
         
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
#' @param pmat      An \eqn{n \times m} matrix of probabilities. \eqn{n} is the number of independent trials.
#'                  \eqn{m} is the number of categories.
#'                  Each row of pmat describes the success probability for the corresponding
#'                  trial and it should add up to 1.
#' @param method     Character string stands for the method selected by user to 
#'                   compute the probability mass. The method can only be one of 
#'                   the following three: 
#'                   \code{"DFT-CF"},
#'                   \code{"NA"},
#'                   \code{"SIM-ALL"}. 
#' @param B          Number of repetitions in the simulation method. Will be ignored if users do not
#'                   choose \code{"SIM-ALL"} method.
#' @param x          A length \eqn{m} vector \eqn{x = (x_{1},\ldots,x_{m})} for computing 
#'                   \eqn{P(X_{1} \leq x_{1},\ldots, X_{m} \leq x_{m})}.
#' 
#' @details 
#' See Details in \code{dpmd} for the definition of the PMD and the introduction of notations.
#' \code{ppmd} computes the cumulative distribution function by adding all probability 
#' mass points within hyper-dimensional space limited by \code{x}. 
#' 
#' \code{"DFT-CF"} is an exact method 
#' to calculate all mass points of Poisson-Multinomial Distributions
#' via FFT algorithm.
#' \code{"SIM-ALL"} is a simulation method using a naive simulation scheme to 
#' calculate the whole probability mass function.
#' \code{"NA"} is an approximation method using Normal approximation method.
#' 
#' @return 
#' The value of \eqn{P(X_{1} \leq x_{1},\ldots, X_{m} \leq x_{m})} of given 
#' \eqn{x = (x_{1},\ldots, x_{m})}.
#' 
#' @examples
#' pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#' x <- matrix(c(3,2,1,3),nrow=1)
#'
#' ppmd(pmat = pp, x = x)
#' ppmd(pmat = pp, x = x, method = "NA")
#' ppmd(pmat = pp, x = x, method = "SIM-ALL", B = 1e3)
#' @export
ppmd = function(pmat,x,method="DFT-CF",B=1e3){
  chck = pmat.check(pmat,x)
  if(chck!=1){stop(chck)}
  
  nn = nrow(pmat)
  mm = ncol(pmat)
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
           if(nrow(points.pos)!=0){
            for(i in 1:nrow(points.pos)){
              prob = prob + dpmd(pmat, x = as.matrix(points.pos[i,]), method="NA")
            }
           }
         })
  return(prob)
}
########################################################################################
#' @title Poisson-Multinomial Distribution Random Number Generator
#' @description Generating random samples from Poisson-Multinomial distribution based on a given success probability matrix.
#'  
#' @param pmat      The \eqn{n \times m} success probability matrix, where \eqn{n} is the number of independent trials and \eqn{m} is the number of categories.
#'                  Each row of pmat describes the success probability for the corresponding
#'                  trial, which adds up to 1.
#' @param s         The number of samples to be generated.
#' 
#' @return 
#' An \eqn{s \times m} matrix of samples, each row stands for one sample from the PMD with success probability matrix \code{pmat}.
#' 
#' @examples 
#' pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
#'  
#' rpmd(pmat = pp, s = 5)
#' 
#' @export
rpmd = function(pmat, s=1){
  chck = pmat.check(pmat)
  if(chck!=1){stop(chck)}
  mm = ncol(pmat)
  rnd = matrix(NA,nrow = s,ncol = mm)
  for(i in 1:s){
    rnd[i,] = t(rpmd_arma(pmat))
  }
  return(rnd)
}



