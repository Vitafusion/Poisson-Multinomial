#' @useDynLib pmd
load.bino = function(){
  bino = read.table('bino.txt')
  bino = bino[,-1]
  colnames(bino) = c('n','MAE','TAE')
  return(bino)
}


bino.plot = function(bino){
  p = list()
  scaleFUN <- function(x) sprintf("%.6f", x)
  p[[1]] = ggplot(data = bino,aes(x=n,y=MAE)) +
    geom_path()+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN)
  p[[2]] = ggplot(data = bino,aes(x=n,y=TAE)) +
    geom_path()+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN)
  return(p)
}


data.generator = function(cov_num,N,m,beta){
  parm = matrix(beta,nrow=m-1,ncol=cov_num)
  x1 = rnorm(N)
  x2 = rnorm(N)
  x = as.matrix(cbind(x1,x2))
  P = exp(x%*%parm)/(rowSums(exp(x%*%parm))+1)
  P = cbind(P,1/(rowSums(exp(x%*%parm))+1))
  choice = apply(P,1,rmultinom,n=1,size=1)
  choice = t(choice)
  y = apply(choice,1,function(x) which(x==1))
  y = as.data.frame(y)
  return(list(x,y,P))

}

sigma.calcu = function(p){
  n = nrow(p)
  m = ncol(p)-1
  if(n==1) pp = t(as.matrix(p[,1:m])) else pp = as.matrix(p[,1:m])
  sig = matrix(0,m,m)
  for (i in 1:n) {
    sig = sig + diag(pp[i,],nrow = m) - pp[i,]%*%t(pp[i,])
  }
  return(sig)
}

mu.calcu = function(p){
  n = nrow(p)
  m = ncol(p)-1
  if(n==1)  pp = t(as.matrix(p[,1:m])) else pp = as.matrix(p[,1:m])
  mu = matrix(0, nrow = 1, ncol = ncol(pp))
  for (i in 1:n) {
    mu = mu + pp[i,]
  }
  return(as.vector(mu))
}


original.data.trans = function(){
  ori_m3 = as.data.frame(ori_m3[,2:3])
  colnames(ori_m3) = c("n","MAE")
  ori_m3$type = as.factor(3)
  ori_m4 = as.data.frame(ori_m4[,2:3])
  colnames(ori_m4) = c("n","MAE")
  ori_m4$type = as.factor(4)
  ori_m5 = as.data.frame(ori_m5[,2:3])
  colnames(ori_m5) = c("n","MAE")
  ori_m5$type = as.factor(5)
  ori_m6 = as.data.frame(ori_m6[,2:3])
  colnames(ori_m6) = c("n","MAE")
  ori_m6$type = as.factor(6)
  ori_m7 = as.data.frame(ori_m7[,2:3])
  colnames(ori_m7) = c("n","MAE")
  ori_m7$type = as.factor(7)
  ori_m8 = as.data.frame(ori_m8[,2:3])
  colnames(ori_m8) = c("n","MAE")
  ori_m8$type = as.factor(8)
  dat = rbind(ori_m3,ori_m4,ori_m5,ori_m6,ori_m7,ori_m8)
  return(dat)
}

test.accuracy.dat.trans = function(dat,m){
  dat = as.data.frame(dat[,-1])
  colnames(dat) = c('n','(p_min,p_max)','MAE','TAE')
  dat$m = rep(m,nrow(dat))
  dat$n = as.numeric(dat$n)
  dat$MAE = as.numeric(dat$MAE)
  dat$TAE = as.numeric(dat$TAE)
  dat$m = as.factor(dat$m)
  dat$type = as.factor(dat$m)
  dat$MAE_trans = dat$MAE*(sqrt(2*pi)^m)
  return(dat)
}


test.accuracy.dat.trans.r = function(dat,m){
  dat = as.data.frame(dat[,-1])
  colnames(dat) = c('n','MAE','TAE')
  dat$m = rep(m,nrow(dat))
  dat$n = as.numeric(dat$n)
  dat$MAE = as.numeric(dat$MAE)
  dat$TAE = as.numeric(dat$TAE)
  dat$m = as.factor(dat$m)
  dat$type = as.factor(dat$m)
  dat$MAE_trans = dat$MAE*((sqrt(2*pi))^m)
  return(dat)
}

compar.all.ori.mae.plot = function(m){
  expr=paste0('m','=',m)
  scaleFUN <- function(x) sprintf("%.4f", x)
  all_mae_dat %>% filter(type==m) %>%
    ggplot(aes(x=n,y=MAE,colour = method,group = method)) +
    geom_path(size=1,aes(linetype=method))+ ggtitle(expr) +
    scale_y_continuous(trans = log_trans(),labels = scaleFUN)
}


compar.ori.mae.plot = function(m,meth = 'norm'){
  expr=paste0('m','=',m,'(method=',meth,')')
  scaleFUN <- function(x) sprintf("%.4f", x)
  all_mae_dat %>% filter(type==m,method==meth|method=='origi') %>%
  ggplot(aes(x=n,y=MAE,colour = method,group = method)) + ggtitle(expr) +
    geom_path(size=1,aes(linetype=method))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN)
}


simu.acc.mae.plot = function(){

  scaleFUN <- function(x) sprintf("%.4f", x)
  ggplot(data = simu_acc_all_dat,aes(x=n,y=MAE,colour = m,group = type)) +
    geom_path(size=1,aes(linetype=m))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = "Simulation Method Accuracy Compared With DF-FFT",
         subtitle = "Simulation Repeat = 10^4, use MAE as criterion",
         caption = "Different colors and types of lines representing different m. The MAE is under log transformation")

}


simu.acc.tae.plot = function(){

  scaleFUN <- function(x) sprintf("%.4f", x)
  ggplot(data = simu_acc_all_dat,aes(x=n,y=TAE,colour = m,group = type)) +
    geom_path(size=1,aes(linetype=m))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = "Simulation Method Accuracy Compared With DF-FFT",
         subtitle = "Simulation Repeat = 10^4, use TAE as criterion",
         caption = "Different colors and types of lines representing different m. The TAE is under log transformation")

}


norm.acc.mae.plot = function(){

  scaleFUN <- function(x) sprintf("%.4f", x)
  ggplot(data = norm_acc_all_dat,aes(x=n,y=MAE,colour = m,group = type)) +
    geom_path(size=1,aes(linetype=m))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = "Normal Method Accuracy Compared With DF-FFT",
         subtitle = "Use MAE as criterion",
         caption = "Different colors and types of lines representing different m. The MAE is under log transformation")

}


norm.acc.tae.plot = function(){

  scaleFUN <- function(x) sprintf("%.4f", x)
  ggplot(data = norm_acc_all_dat,aes(x=n,y=TAE,colour = m,group = type)) +
    geom_path(size=1,aes(linetype=m))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = "Normal Method Accuracy Compared With DF-FFT",
         subtitle = "Use TAE as criterion",
         caption = "Different colors and types of lines representing different m. The TAE is under log transformation")

}

norm.acc.mae.trans.plot = function(){

  scaleFUN <- function(x) sprintf("%.4f", x)
  ggplot(data = norm_acc_all_dat,aes(x=n,y=MAE_trans,colour = m,group = type)) +
    geom_path(size=1,aes(linetype=m))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = "Normal Method Accuracy Compared With DF-FFT",
         subtitle = "Use transfered MAE as criterion",
         caption = "Different colors and types of lines representing different m. The MAE_transfered is under log transformation")

}

simu.acc.mae.trans.plot = function(){

  scaleFUN <- function(x) sprintf("%.4f", x)
  ggplot(data = simu_acc_all_dat,aes(x=n,y=MAE_trans,colour = m,group = type)) +
    geom_path(size=1,aes(linetype=m))+
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = "Simulation Method Accuracy Compared With DF-FFT",
         subtitle = "Use transfered MAE as criterion",
         caption = "Different colors and types of lines representing different m. The MAE_transfered is under log transformation")

}

comparison.plot = function(cat){
  str = paste0("Simulation Method vs Normal approximation(m=",cat,")")
  scaleFUN <- function(x) sprintf("%.4f", x)
  all_dat %>% filter(m==cat) %>% ggplot() +
    geom_path(aes(x=n,y=MAE,colour=method,group=method,linetype=method)) +
    scale_y_continuous(trans = log_trans(),labels = scaleFUN) +
    labs(title = str,
         subtitle = "Simulation Repeat = 10^4, use MAE as criterion",
         caption = "Different colors and types of lines representing different methods. The MAE is under log transformation")

}


make.data.equal.groups = function(dat, groups_number,category_name=c("setosa","versicolor","virginica"),category_column = 'Species'){
  n = nrow(dat)
  group = list()
  if(n %% groups_number != 0)
    stop("invalid number for grouping.")
  out = matrix(NA,nrow = groups_number,ncol = length(category_name))
  ss <- sample(1:groups_number,size=n,replace=TRUE,prob=rep(1/groups_number,groups_number))
  for (i in 1:groups_number) {
    group[[i]] = dat[ss==i,]
    expr = paste(category_column)
    for(j in 1:length(category_name))
      out[i,j] = nrow(group[[i]][which(group[[i]][,expr]==category_name[j]),])
  }
  return(list(group,out))
}



make.data.groups = function(dat, groups_number,category_name=c("setosa","versicolor","virginica"),category_column = 'Species'){
  n = nrow(dat)
  group = list()
  k=1
  if(n %% groups_number != 0)
    stop("invalid number for grouping.")
  out = matrix(NA,nrow = groups_number,ncol = length(category_name))
  ss <- sample(1:groups_number,size=n,replace=TRUE,prob=rep(1/groups_number,groups_number))
  for (i in 1:groups_number) {
    if(nrow(dat[ss==i,])!=0){
      group[[k]] = dat[ss==i,]
      expr = paste(category_column)
      for(j in 1:length(category_name)){
        out[k,j] = nrow(group[[k]][which(group[[k]][,expr]==category_name[j]),])
      }
      k = k + 1
    }
  }
  out = na.omit(out)
  return(list(group,out))
}

cal_pmatrix = function(parm, x_mat, category_number){
  x_mat = as.matrix(x_mat)
  n = ncol(x_mat)
  m = category_number
  parm = matrix(parm,n,m-1)
  P = matrix(0, nrow = nrow(x_mat), ncol = category_number)
  P = exp(x_mat%*%parm)/(rowSums(exp(x_mat%*%parm))+1)
  P = cbind(P,1/(rowSums(exp(x_mat%*%parm))+1))
  return(P)
}

####################################### Using groups to analysis likelihoods
point.loglik.calcu = function(pp,count_result){
  #browser()
  res = pmd(pp)
  l.vec = compute_l.vec(count_result)
  l.vec = paste(l.vec,collapse = ",")
  expr0 = "res["
  expr = paste0(expr0,l.vec,"]")
  likel = eval(parse(text=expr))
  log_lik = log(likel)
  #browser()
  return(log_lik)
}

toltal.loglik.calcu = function(parm,result_mat = count_result,group = groups,cat_number = category_number,covariate_num) {
  minus_log_lik = 0
  for (i in 1:nrow(result_mat)) {
    x = group[[i]]
    x_mat = x[,1:covariate_num]
    n = nrow(x)
    result = result_mat[i,]
    P = cal_pmatrix(parm, x_mat, cat_number)
    minus_log_lik = minus_log_lik - point.loglik.calcu(P,result)
    #browser()
  }
  #browser()
  return(minus_log_lik)
}

toltal.loglik.calcu.ai4i = function(parm,result_mat = count_result,group = groups,cat_number = category_number,
                                    covariate_name = c("Air temperature [K]","Process temperature [K]","Rotational speed [rpm]","Torque [Nm]")) {
  minus_log_lik = 0
  for (i in 1:nrow(result_mat)) {
    x = group[[i]]
    x_mat = x[,covariate_name]
    n = nrow(x)
    result = result_mat[i,]
    P = cal_pmatrix(parm, x_mat, cat_number)
    minus_log_lik = minus_log_lik - point.loglik.calcu(P,result)
    #browser()
  }
  #browser()
  return(minus_log_lik)
}


toltal.loglik.calcu.driver = function(parm,result_mat = count_result,group = groups,cat_number = category_number,
                                    covariate_name = c("event.start","reaction.start","impact.time","event.end")) {
  minus_log_lik = 0
  for (i in 1:nrow(result_mat)) {
    x = group[[i]]
    x_mat = x[,covariate_name]
    n = nrow(x)
    result = result_mat[i,]
    P = cal_pmatrix(parm, x_mat, cat_number)
    minus_log_lik = minus_log_lik - point.loglik.calcu(P,result)
    #browser()
  }
  #browser()
  return(minus_log_lik)
}


load.driver = function(){
  driver = read.csv(file = 'driver.csv')
  driver = driver[,-1]
  return(driver)
}

driver.prep = function(){
  driver = read.csv(file = 'shrp.csv')
  driver = driver[,c(4,6,7,8,9)]
  colnames(driver) = c('event severity','event start','reaction start','impact time','event end')
  driver$`event severity`=as.factor(driver$`event severity`)
  for(i in 2:5){
    driver[,i] = as.character(driver[,i])
    driver[,i] = as.numeric(driver[,i])
    driver[,i] = scale(driver[,i])
  }
  driver = na.omit(driver)
  driver = driver[1:8000,]
  write.csv(driver,file = 'driver.csv')
}

ai4i.prep = function(){
  ai4i = read.table(file = 'dat_ai4i.csv',sep = ',')
  colnames(ai4i) = c('NA','Product ID','Type','Air temperature [K]','Process temperature [K]','Rotational speed [rpm]','Torque [Nm]')
  #colnames(ai4i) = ai4i[1,]
  ai4i = ai4i[-1,2:7]
  ai4i$Type=as.factor(ai4i$Type)
  for(i in 3:6){
    ai4i[,i] = as.character(ai4i[,i])
    ai4i[,i] = as.numeric(ai4i[,i])
    ai4i[,i] = scale(ai4i[,i])
  }
  return(ai4i)
}

compute_l.vec = function(result){
  m = length(result)
  res = result[1:(m-1)]
  l.vec = rep(0,m-1)
  for (i in 1:(m-1)) {
    l.vec[i] = as.numeric(res[i]) + 1
  }
  return(l.vec)
}

#' @title a function to compute l.vec
#' @param m,cn.vec,k categories and cn.vec
#' @return l.vec a vector
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

#' @title a function to create pp matrix
#' @param n,m dimension
#' @return pp a result array
#' @export
pmatrix <- function(n,m){
  p <- matrix(0,nrow = n,ncol = m,byrow = T)
  for (i in 1:n) {
    r <- runif(m)
    r <- r/sum(r) #generate row
    p[i,] <- r
  }
  return(p)
}
#' @title main function
#' @param pp a matrix of input
#' @return res a result array
#' @export
pmd <-function(pp,method="DFT-CF",vec=c(0,0,0,0,0),t=100)
{

  if(any(pp<0)|any(pp>1))
  {
    stop("invalid values in pp.")
  }


  switch(method,
         "DFT-CF"={
           mm=ncol(pp) # m categories
           nn=nrow(pp) # n people

           nn.vec=rep(nn+1, mm-1)
           l.vec=rep(0, mm-1)
           cn.vec=cumprod(nn.vec)
           cn.vec=c(1, cn.vec[-(mm-1)])
           cn.vec=cn.vec[length(cn.vec):1]
           cn.vec=as.integer(cn.vec)

           nnt=prod(nn.vec)
           res0=double(nnt)

           #browser()

           tmp=.C("pmn_mdfft", as.double(res0), as.integer(nnt), as.integer(nn), as.integer(mm), as.double(pp), as.integer(nn.vec), as.integer(l.vec), as.integer(cn.vec))
           res0=tmp[[1]]
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

         },
         "simulation"=    {
           mm=ncol(pp) # m categories
           nn=nrow(pp) # n people
           nn.vec=rep(nn+1, mm-1)
           l.vec=rep(0, mm-1)
           cn.vec=cumprod(nn.vec)
           cn.vec=c(1, cn.vec[-(mm-1)])
           cn.vec=cn.vec[length(cn.vec):1]
           cn.vec=as.integer(cn.vec)
           nnt=prod(nn.vec)
           res0=double(nnt)
           temp=.C("pmd_simulation", as.double(res0),
                   as.integer(nnt) , as.integer(nn), as.integer(mm), as.double(pp),
                   as.integer(nn.vec), as.integer(l.vec), as.integer(cn.vec), as.integer(t))
           res=round(temp[[1]],10)
         },
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
           #temp=.C("pmd_normal",as.double(res), as.integer(nn), as.integer(mm), as.double(pp),as.integer(x_vec))
           #res = temp[[1]]

         }

  )

  return(res)
}
