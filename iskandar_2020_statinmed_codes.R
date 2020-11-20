##############################################################################
#This code is to conduct p-box analyses (1D, parallel version)
#author: Rowan Iskandar
#email: rowan.iskandar@gmail.com
#last edited: 20112020
#The code for the numerical exercise is available under a GNU GPL license and can be found at https://github.com/rowaniskandar/CM_ME_plosone.
#loading relevant libraries
#library(arules)
library(nloptr) #optimization
library(doParallel)  
library(future)
library(future.apply)
library(purrr)
library(truncnorm)
library(gmailr)
library(zip)
library(support)
#library(plotly)
library(ggplot2)
currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  

##############################################################################
#apply nlopt for a given set of params (find pbox for each sampled param pdf)
fn_bbmod1_simp_par_inner <- function(params_pdf_in, params_c_in, index_slice,pbox.int,index.array){
  dim_slice <- length(index_slice)
  result_list <- as.list(1:(dim_slice*3))
  dim(result_list) <- c(dim_slice,3)
  for (i in 1:dim_slice){
    k <- index_slice[i]
    lbound <- c(pbox.int[3,index.array[k,1],1],pbox.int[5,index.array[k,2],1])
    ubound <- c(pbox.int[3,index.array[k,1],2],pbox.int[5,index.array[k,2],2])
    #max
    optim.max <- directL(fn_bbmod1_simp_max,lbound,ubound ,
                         nl.info = TRUE, control=list(xtol_rel=1e-6, maxeval=1000), 
                         params_pdf = params_pdf_in,
                         params_c = params_c_in,
                         randomized=TRUE)
    #optim.results[k,1] <- optim.max$value
    #min
    optim.min <- directL(fn_bbmod1_simp_min,lbound,ubound ,
                         nl.info = TRUE, control=list(xtol_rel=1e-6, maxeval=1000),
                         params_pdf = params_pdf_in,
                         params_c = params_c_in,
                         randomized=TRUE)
    #optim.results[k,2] <- optim.min$value
    #weight
    #weight
    w1 <- pbox.int[3,index.array[k,1],3]
    w2 <- pbox.int[5,index.array[k,2],3]
    w <- w1*w2
    #optim.results[k,3] <- w1*w2
    result_list[i,] <- c(optim.min$value,optim.max$value, w)
  }
  return(result_list)
}
#main cohort model (take two different sets of parameters:  pdf, constant) or no pbox version
fn_bbmod1_simp <-function(params_pdf, params_c){
  #parameter values for transition rates in the 4-state model (note: these are rates, not probabilities)
  c12 <- params_c[1]
  c13 <- params_pdf[[1]]
  c14 <- params_pdf[[2]]
  c23 <- params_pdf[[3]]
  c24 <- params_pdf[[4]]
  c34 <- params_pdf[[5]]
  
  popsize <- 10000
  ns <- 4 #number of states
  tpoints <- seq(0,50,1)
  
  params <- c(c12 ,
              c13 ,
              c14 ,
              c23 ,
              c24 ,
              c34 )
  #params <- as.vector(params)
  init <-c(popsize,0,0,0)
  p11 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    p<-exp(-s1*t)
    return(p)}
  p12 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    s2 <-params[4]+params[5]
    p<-params[1]*(exp(-s2*t)-exp(-s1*t))/(s1-s2)
    return(p)}
  p13 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    s2 <-params[4]+params[5]
    s3 <-params[6]
    p1<-params[2]*(exp(-s3*t)-exp(-s1*t))/(s1-s3)
    p2<-(params[1]*params[4]*((s1-s2)*exp(-s3*t)-(s1-s3)*exp(-s2*t)+(s2-s3)*exp(-s1*t))) /((s1-s2)*(s1-s3)*(s2-s3))
    p<-p1+p2
    return(p)}
  p14 <- function(params,t) {
    p1 <- p11(params,t)
    p2 <- p12(params,t)
    p3 <- p13(params,t)
    p <- 1-p1-p2-p3
    return(p)}
  
  p21 <- 0
  p22 <-function(params,t) {
    s2 <-params[4]+params[5]
    p<-exp(-s2*t)
    return(p)}
  
  p23 <-function(params,t) {
    s2 <-params[4]+params[5]
    s3 <- params[6]
    p<-params[4]*(exp(-s3*t)-exp(-s2*t))/(s2-s3)
    return(p)}
  p24 <- function(params,t) {
    p1 <- 0
    p2 <- p22(params,t)
    p3 <- p23(params,t)
    p <- 1-p1-p2-p3
    return(p)}
  
  p31 <- 0
  p32 <- 0
  
  p33 <-function(params,t) {
    s3 <-params[6]
    p<-exp(-s3*t)
    return(p)}
  
  p34 <- function(params,t) {
    p <- 1-p33(params,t)
    return(p)}
  
  p41 <- 0
  p42 <- 0
  p43 <- 0
  p44 <- 1
  
  prob <-function(params,t){
    P <- matrix(c(p11(params,t) ,p21,p31,p41,p12(params,t) ,p22(params,t) ,p32,p42,p13(params,t) ,p23(params,t) ,p33(params,t) ,p43,p14(params,t) ,p24(params,t) ,p34(params,t) ,p44 ),nrow=4,ncol=4)
    return(P)
  }
  CM.trace <- matrix(rep(0,length(tpoints)*4),nrow=length(tpoints),ncol=4)
  p.temp <- 0
  for(i in 1:length(tpoints)){
    # for(j in 1:i){
    #   if(j==1){
    #     p.temp <- prob(params,1)
    #   }
    #   else{
    #     p.temp <- p.temp%*%prob(params,j)
    #   }
    # }
    #CM.trace[i,] <- init%*%p.temp
    CM.trace[i,] <- init%*%prob(params,i)
  }
  
  CM.trace_d <- CM.trace
  CM.trace_d[1,] <-c(popsize, 0 ,0 ,0)
  for (t in 2:length(tpoints)){
    CM.trace_d[t,] <- CM.trace[t-1,]
  }
  
  data.CM <- as.data.frame(CM.trace_d)
  
  LE <- sum(data.CM$V1+data.CM$V2+data.CM$V3)/popsize
  return(LE)
}

#main cohort model (take three different sets of parameters: pbox, pdf, constant)
fn_bbmod1_simp_min <-function(params_pbox, params_pdf, params_c){
  #parameter values for transition rates in the 4-state model (note: these are rates, not probabilities)
  c12 <- params_c[1]
  c13 <- params_pdf[[1]]
  c14 <- params_pbox[1]
  c23 <- params_pdf[[2]]
  c24 <- params_pbox[2]
  c34 <- params_pdf[[3]]
  
  popsize <- 10000
  ns <- 4 #number of states
  tpoints <- seq(0,50,1)
  
  params <- c(c12 ,
              c13 ,
              c14 ,
              c23 ,
              c24 ,
              c34 )
  #params <- as.vector(params)
  init <-c(popsize,0,0,0)
  p11 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    p<-exp(-s1*t)
    return(p)}
  p12 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    s2 <-params[4]+params[5]
    p<-params[1]*(exp(-s2*t)-exp(-s1*t))/(s1-s2)
    return(p)}
  p13 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    s2 <-params[4]+params[5]
    s3 <-params[6]
    p1<-params[2]*(exp(-s3*t)-exp(-s1*t))/(s1-s3)
    p2<-(params[1]*params[4]*((s1-s2)*exp(-s3*t)-(s1-s3)*exp(-s2*t)+(s2-s3)*exp(-s1*t))) /((s1-s2)*(s1-s3)*(s2-s3))
    p<-p1+p2
    return(p)}
  p14 <- function(params,t) {
    p1 <- p11(params,t)
    p2 <- p12(params,t)
    p3 <- p13(params,t)
    p <- 1-p1-p2-p3
    return(p)}
  
  p21 <- 0
  p22 <-function(params,t) {
    s2 <-params[4]+params[5]
    p<-exp(-s2*t)
    return(p)}
  
  p23 <-function(params,t) {
    s2 <-params[4]+params[5]
    s3 <- params[6]
    p<-params[4]*(exp(-s3*t)-exp(-s2*t))/(s2-s3)
    return(p)}
  p24 <- function(params,t) {
    p1 <- 0
    p2 <- p22(params,t)
    p3 <- p23(params,t)
    p <- 1-p1-p2-p3
    return(p)}
  
  p31 <- 0
  p32 <- 0
  
  p33 <-function(params,t) {
    s3 <-params[6]
    p<-exp(-s3*t)
    return(p)}
  
  p34 <- function(params,t) {
    p <- 1-p33(params,t)
    return(p)}
  
  p41 <- 0
  p42 <- 0
  p43 <- 0
  p44 <- 1
  
  prob <-function(params,t){
    P <- matrix(c(p11(params,t) ,p21,p31,p41,p12(params,t) ,p22(params,t) ,p32,p42,p13(params,t) ,p23(params,t) ,p33(params,t) ,p43,p14(params,t) ,p24(params,t) ,p34(params,t) ,p44 ),nrow=4,ncol=4)
    return(P)
  }
  CM.trace <- matrix(rep(0,length(tpoints)*4),nrow=length(tpoints),ncol=4)
  p.temp <- 0
  for(i in 1:length(tpoints)){
    CM.trace[i,] <- init%*%prob(params,i)
  }
  
  CM.trace_d <- CM.trace
  CM.trace_d[1,] <-c(popsize, 0 ,0 ,0)
  for (t in 2:length(tpoints)){
    CM.trace_d[t,] <- CM.trace[t-1,]
  }
  
  data.CM <- as.data.frame(CM.trace_d)
  
  LE <- sum(data.CM$V1+data.CM$V2+data.CM$V3)/popsize
  return(LE)
}

fn_bbmod1_simp_max <-function(params_pbox, params_pdf, params_c){
  #parameter values for transition rates in the 4-state model (note: these are rates, not probabilities)
  c12 <- params_c[1]
  c13 <- params_pdf[[1]]
  c14 <- params_pbox[1]
  c23 <- params_pdf[[2]]
  c24 <- params_pbox[2]
  c34 <- params_pdf[[3]]
  
  popsize <- 10000
  ns <- 4 #number of states
  tpoints <- seq(0,50,1)
  
  params <- c(c12 ,
              c13 ,
              c14 ,
              c23 ,
              c24 ,
              c34 )
  #params <- as.vector(params)
  init <-c(popsize,0,0,0)
  p11 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    p<-exp(-s1*t)
    return(p)}
  p12 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    s2 <-params[4]+params[5]
    p<-params[1]*(exp(-s2*t)-exp(-s1*t))/(s1-s2)
    return(p)}
  p13 <- function(params,t) {
    s1 <-params[1]+params[2]+params[3]
    s2 <-params[4]+params[5]
    s3 <-params[6]
    p1<-params[2]*(exp(-s3*t)-exp(-s1*t))/(s1-s3)
    p2<-(params[1]*params[4]*((s1-s2)*exp(-s3*t)-(s1-s3)*exp(-s2*t)+(s2-s3)*exp(-s1*t))) /((s1-s2)*(s1-s3)*(s2-s3))
    p<-p1+p2
    return(p)}
  p14 <- function(params,t) {
    p1 <- p11(params,t)
    p2 <- p12(params,t)
    p3 <- p13(params,t)
    p <- 1-p1-p2-p3
    return(p)}
  
  p21 <- 0
  p22 <-function(params,t) {
    s2 <-params[4]+params[5]
    p<-exp(-s2*t)
    return(p)}
  
  p23 <-function(params,t) {
    s2 <-params[4]+params[5]
    s3 <- params[6]
    p<-params[4]*(exp(-s3*t)-exp(-s2*t))/(s2-s3)
    return(p)}
  p24 <- function(params,t) {
    p1 <- 0
    p2 <- p22(params,t)
    p3 <- p23(params,t)
    p <- 1-p1-p2-p3
    return(p)}
  
  p31 <- 0
  p32 <- 0
  
  p33 <-function(params,t) {
    s3 <-params[6]
    p<-exp(-s3*t)
    return(p)}
  
  p34 <- function(params,t) {
    p <- 1-p33(params,t)
    return(p)}
  
  p41 <- 0
  p42 <- 0
  p43 <- 0
  p44 <- 1
  
  prob <-function(params,t){
    P <- matrix(c(p11(params,t) ,p21,p31,p41,p12(params,t) ,p22(params,t) ,p32,p42,p13(params,t) ,p23(params,t) ,p33(params,t) ,p43,p14(params,t) ,p24(params,t) ,p34(params,t) ,p44 ),nrow=4,ncol=4)
    return(P)
  }
  CM.trace <- matrix(rep(0,length(tpoints)*4),nrow=length(tpoints),ncol=4)
  p.temp <- 0
  for(i in 1:length(tpoints)){
    # for(j in 1:i){
    #   if(j==1){
    #     p.temp <- prob(params,1)
    #   }
    #   else{
    #     p.temp <- p.temp%*%prob(params,j)
    #   }
    # }
    #CM.trace[i,] <- init%*%p.temp
    CM.trace[i,] <- init%*%prob(params,i)
  }
  
  CM.trace_d <- CM.trace
  CM.trace_d[1,] <-c(popsize, 0 ,0 ,0)
  for (t in 2:length(tpoints)){
    CM.trace_d[t,] <- CM.trace[t-1,]
  }
  
  data.CM <- as.data.frame(CM.trace_d)
  
  LE <- sum(data.CM$V1+data.CM$V2+data.CM$V3)/popsize
  return(-LE)
}

#pbox function
pbox_minmax <-function(x,type,min,max){
  if (type=="lower"){
    if (x<max){
      p=0
    }
    else{
      p=1
    }
  }
  else{
    if (x<min){
      p=0
    }
    else{
      p=1
    }
  }
  return(p)
}
#given min, max, mean
pbox_mean <- function(x,type,min,max,mean){
  if (type=="lower"){
    if (x<mean){
      p=0
    }
    else if(mean <= x & x<max){
      p=(x-mean)/(x-min)
    }
    else {
      p=1
    }
  }
  else{
    if (x<mean){
      p=(max-mean)/(max-x)
    }
    else{
      p=1
    } 
  }
  return(p)
}
#given min, max, median
pbox_median <- function(x,type,min,max,median){
  if (type=="lower"){
    if (x<median){
      p=0
    }
    else if(max > x & x>=median){
      p=1/2
    }
    else{
      p=1
    }
  }
  else{
    if (x<min){
      p=0
    }
    else if(min<=x & x<median){
      p=1/2
    }
    else{
      p=1
    }
  }
  return(p)
}
#inverse p-box functions for simulating p-box
pbox_mean_inv <- function(u,type,min,max,mean){
  crit <- (max-mean)/(max-min)
  if (type=="upper"){
    #need to double check
    if(u==1){
      x=max
      #x=mean
    }
    else if(crit<u && u<1){
      x=(u*max-max+mean)/u
    }
    else if(u==crit){
      x=min
    }
    else if(u<crit){
      x=min
    }
  }
  else{
    if (u==1){
      x=max
    }
    else if (0<u && u< crit){
      x=(u*min-mean)/(u-1)
    }
    else if(u==0){
      x=min #not sure, min <= x < mean will yield u=0
    }
    else if(crit <=u && u<1){ #not sure about the condition
      x=max
    }
  }
  return(x)
}


pbox_std_inv <- function(u,type, min, max, mean, std){
  F1 <- (std^2)/((mean-min)^2+std^2)
  print(F1)
  F2 <- (max-mean)^2/((max-mean)^2 + std^2)
  print(F2)
  if (type=="lower"){
    if (u < F1 ){
      x=min
    }
    else if ( F1 <=u  && u < F2 ){
      x=mean-std*sqrt((1-u)/u)
    }
    else {
      num <- (max-mean)*(mean-min)-std^2
      denom <- u*(max-min)+mean-max
      x <- max - num/denom
    }
  }
  else {
    if (u <= F1 ){
      num <- (max-mean)*(mean-min)-std^2
      denom <- max-mean-u*(max-min)
      x<- min + num/denom
    }
    else if ( F1 < u  && u <= F2 ){
      x=mean+std*sqrt(u/(1-u))
    }
    else {
      x=max
    }
  }
  return(x)
}
##############################################################################
#PBA starts here
#4-states model
#6 parameters
c12 <- 0.05 #const x1
c13 <- 0.01 #pdf x2
c14 <- 0.001 #pbox x3
c23 <- 0.1 #pdf x4
c24 <- 0.05 #pbox x5
c34 <- 1 #pdf x6

n.MC <- 1 #number of MC sampling of x2, x4, x6
n.sp <- 200 #number of support point sampling of x2, x4, x6
num.pbox.int <- 5 #number of p-box discretization for UP
mean.std <- 1/3

sp.indicator <- 0 #use support point (1) or MC (0)
all.MC <- 1 #1: no pbox
###########################################################################
#pbox set up
#different types of data
#data for all params
#columns: min, max, median, mean, std
params.table <- matrix(rep(0,6*5),nrow=6,ncol=5)
# as of 2/10/2020, we compare two scenarios (x1: constant, x2-x6: known cdf)
# vs/ (x1: constant, x3,x5: pbox, others: cdf)
params.table[1,] = c(0,10,9999,c12,9999)
params.table[2,] = c(0,10,9999,c13,mean.std*c13)
params.table[3,]= c(0,10,9999,c14,mean.std*c14)
params.table[4,]= c(0,10,9999,c23,mean.std*c23)
params.table[5,]= c(0,10,9999,c24,mean.std*c24)
params.table[6,]= c(0,10,9999,c34,mean.std*c34)

#1. uncertainty quantification
#discretization
#over the cdf values (u)
#columns: (1) lower bound - interval (2) upper bound - interval (3) width of interval
inv.int <- matrix(rep(0,(num.pbox.int-1)*3),nrow=num.pbox.int-1,ncol=3) 
dummy <- seq(0,1,length.out=num.pbox.int)
for (i in 1:(num.pbox.int-1)){
  inv.int[i,1]=dummy[i]
  inv.int[i,2]=dummy[i+1]
  inv.int[i,3]=dummy[i+1]-dummy[i]
}
#over x values, using inverse pbox functions
#columns: (1) variable, (2) intervals, (3) lower or upper bounds or weight
pbox.int <-array(rep(0,6*(num.pbox.int-1)*3), dim=c(6, num.pbox.int-1, 3))
#pbox.int <- matrix(rep(0,(num.pbox.int-1)*2),nrow=num.pbox.int-1,ncol=3)
#x5
for (i in 1:(num.pbox.int-1)){
  pbox.int[5,i,1]=pbox_std_inv(inv.int[i,1],"lower",params.table[5,1],params.table[5,2],params.table[5,4],params.table[5,5])
  pbox.int[5,i,2]=pbox_std_inv(inv.int[i,2],"upper",params.table[5,1],params.table[5,2],params.table[5,4],params.table[5,5])
  pbox.int[5,i,3]=inv.int[i,3]
}
#x3
for (i in 1:(num.pbox.int-1)){
  pbox.int[3,i,1]=pbox_std_inv(inv.int[i,1],"lower",params.table[3,1],params.table[3,2],params.table[3,4],params.table[3,5])
  pbox.int[3,i,2]=pbox_std_inv(inv.int[i,2],"upper",params.table[3,1],params.table[3,2],params.table[3,4],params.table[3,5])
  pbox.int[3,i,3]=inv.int[i,3]
}

#2. uncertainty propagation
#create multi-indices (for each p-box parameter)
num.param.pba <- 2 #number of model parameters with p-box
n.indices <- (num.pbox.int-1)**num.param.pba #number of unique indices
input_data <- 1:n.indices
#array to store multi indices
index.array <- matrix(rep(0,n.indices),nrow=n.indices,ncol=num.param.pba+1) #last column to store total weight
k=1
for (i in 1:(num.pbox.int-1)){
  for (j in 1:(num.pbox.int-1)){
    index.array[k,1]=i
    index.array[k,2]=j
    k=k+1
  }
}
#pdf (truncated normal)
# x2 <- rtruncnorm(n.MC, a=params.table[2,1], b=params.table[2,2], mean = params.table[2,4], sd = params.table[2,5])
# x4 <- rtruncnorm(n.MC, a=params.table[4,1], b=params.table[4,2], mean = params.table[4,4], sd = params.table[4,5])
# x6 <- rtruncnorm(n.MC, a=params.table[6,1], b=params.table[6,2], mean = params.table[6,4], sd = params.table[6,5])

#pdf of x2, x4, x6 (all gamma)
#MC sampling method
#note that is x2-x6 for All MC scenario (03102020)
set.seed(710)
x2 <- rgamma(n.MC, shape=(params.table[2,4])^2/params.table[2,5]^2, scale=params.table[2,5]^2/(params.table[2,4]))
x3 <- rgamma(n.MC, shape=(params.table[3,4])^2/params.table[3,5]^2, scale=params.table[3,5]^2/(params.table[3,4]))
x4 <- rgamma(n.MC, shape=(params.table[4,4])^2/params.table[4,5]^2, scale=params.table[4,5]^2/(params.table[4,4]))
x5 <- rgamma(n.MC, shape=(params.table[5,4])^2/params.table[5,5]^2, scale=params.table[5,5]^2/(params.table[5,4]))
x6 <- rgamma(n.MC, shape=(params.table[6,4])^2/params.table[6,5]^2, scale=params.table[6,5]^2/(params.table[6,4]))

#support point method
#note that is x2-x6 for All MC scenario (03102020)
set.seed(710)
dist.param <- vector("list",5)
dist.param[[1]] <- c(params.table[2,4]^2/params.table[2,5]^2,params.table[2,5]^2/(params.table[2,4]))
dist.param[[2]] <- c(params.table[3,4]^2/params.table[3,5]^2,params.table[3,5]^2/(params.table[3,4]))
dist.param[[3]] <- c(params.table[4,4]^2/params.table[4,5]^2,params.table[4,5]^2/(params.table[4,4]))
dist.param[[4]] <- c(params.table[5,4]^2/params.table[5,5]^2,params.table[5,5]^2/(params.table[5,4]))
dist.param[[5]] <- c(params.table[6,4]^2/params.table[6,5]^2,params.table[6,5]^2/(params.table[6,4]))
D <- sp(n.sp,5,dist.str=rep("gamma",5),dist.param=dist.param)
x2.sp <- D$sp[,1]
x3.sp <- D$sp[,2]
x4.sp <- D$sp[,3]
x5.sp <- D$sp[,4]
x6.sp <- D$sp[,5]
list_test <- list()
list_test.sp <- list()
#list for no pbox
#MC list
for (i in 1:n.MC){
  list_test[[i]] <- list(x2[i],x3[i],x4[i],x5[i],x6[i])
}
#sp list
for (i in 1:n.sp){
  list_test.sp[[i]] <- list(x2.sp[i],x3.sp[i],x4.sp[i],x5.sp[i],x6.sp[i])
}
#Pbox scenarios
# dist.param.pbox <- vector("list",3)
# dist.param.pbox[[1]] <- c(params.table[2,4]^2/params.table[2,5]^2,params.table[2,5]^2/(params.table[2,4]))
# dist.param.pbox[[2]] <- c(params.table[4,4]^2/params.table[4,5]^2,params.table[4,5]^2/(params.table[4,4]))
# dist.param.pbox[[3]] <- c(params.table[6,4]^2/params.table[6,5]^2,params.table[6,5]^2/(params.table[6,4]))
# D.pbox <- sp(n.sp,3,dist.str=rep("gamma",3),dist.param=dist.param)
# x2.sp.pbox <- D.pbox$sp[,1]
# x4.sp.pbox <- D.pbox$sp[,2]
# x6.sp.pbox <- D.pbox$sp[,3]

#pbox scenario (only 3 params are varied MC-lly)
list_test.pbox <- list()
list_test.sp.pbox <- list()
#list for pbox
for (i in 1:n.MC){
  list_test.pbox[[i]] <- list(x2[i],x4[i],x6[i])
}
for (i in 1:n.sp){
  list_test.sp.pbox[[i]] <- list(x2.sp[i],x4.sp[i],x6.sp[i])
}

#choose: sp or MC 
if (sp.indicator==1){
  params_pdf <- list_test.sp
  params_pdf.pbox <- list_test.sp.pbox
  
} else{
  params_pdf <- list_test
  params_pdf.pbox <- list_test.pbox
}
params_c <- c12

#non-parallel run if needed
# ptm1 <- proc.time()
# result_list_outer <- lapply(params_pdf,fn_bbmod1_simp_par_inner, params_c, input_data,pbox.int,index.array)
# time1 <- (proc.time()-ptm1)
# ela1 <-time1[1]
##############################################################################
#Run with MC
# num.cores <- detectCores()
# availableCores()
# #access VM - rstudio docker through shell
# cl <- makeClusterPSOCK(
#   # Public IP number of EC2 instance
#   workers=availableCores()-1,
#   dryrun = FALSE,
#   connectTimeout = 120,
#   verbose=TRUE,
#   outfile=NULL,
#   rshlogfile=TRUE
# )
# 
# plan(list(tweak(cluster, workers = cl), multiprocess))
# set.seed(710)
# ptm2 <- proc.time()
# result_list_outer.pbox <- future.apply::future_lapply(params_pdf.pbox,fn_bbmod1_simp_par_inner, params_c, input_data,pbox.int,index.array)
# result_list_outer <- future.apply::future_lapply(params_pdf,fn_bbmod1_simp, params_c)
# time2 <- (proc.time()-ptm2)
# ela2 <-time2[3]
# 
# #saving workspace #always specify (1) computer (AWS, etc) (2) number of PBA intervals
# #(3) number of MCs (4) date
# file_name <- "c4.2xlarge_PBA4_sp25_18082020.rds"
# save.image(file=file_name)
# zip_name <- 'c4.2xlarge_PBA4_sp25_18082020.zip'
# zip(zipfile=zip_name,
#     files=file_name)
# #gmailr
# gm_auth_configure(path = "credentials.json")
# html_msg <-
#   gm_mime() %>%
#   gm_to("rowan.iskandar@gmail.com") %>%
#   gm_from("pbox.output@gmail.com") %>%
#   gm_subject(file_name) %>%
#   gm_text_body("done!")
# msg <- html_msg %>%
#   gm_subject(file_name) %>%
#   gm_attach_file(zip_name)
# gm_create_draft(msg)
# gm_send_message(msg)
# 
# #run this to trigger authentication for the first time using an instance
# gm_auth_configure(path = "credentials.json")
# html_msg <-
#   gm_mime() %>%
#   gm_to("rowan.iskandar@gmail.com") %>%
#   gm_from("pbox.output@gmail.com") %>%
#   gm_subject("test") %>%
#   gm_text_body("done!")
# msg <- html_msg %>%
#   gm_subject("test") 
# gm_create_draft(msg)
# gm_send_message(msg)
# 
# stopCluster(cl)
# #rearranging results
# #columns: (1) min of outcomes, (2) max of outcomes, (3) weights
# #row correspond to a multi-index
# #for comparison MC vs sp
# #first: MC results (all MC and pbox)
# n.data <- n.MC
# results.MC <- array(rep(0,n.data*n.indices*3), c(n.data,n.indices,3))
# results.MC.pbox <- array(rep(0,n.data*n.indices*3), c(n.data,n.indices,3))
# 
# for (i in 1:n.data){
#   z <- 0
#   for (j in 1:3){
#     for (k in 1:n.indices){
#       z <- z+1
#       results.MC[i,k,j] <- result_list_outer[[i]][[z]]
#       results.MC.pbox[i,k,j] <- result_list_outer.pbox[[i]][[z]]
#     }
#   }
# }
# results.MC <- abs(results.MC)
# results.MC.pbox <- abs(results.MC.pbox)
# 
# #derive empirical cdf of outcomes
# min.x <- 0
# max.x <- 40
# n.x <-100
# x.values <- seq(min.x,max.x,length.out=n.x)
# cdf.x.upper.MC <- array(seq(0,n.data*n.x), c(n.data,n.x))
# cdf.x.lower.MC <- array(seq(0,n.data*n.x), c(n.data,n.x))
# cdf.x.upper.MC.pbox <- array(seq(0,n.data*n.x), c(n.data,n.x))
# cdf.x.lower.MC.pbox <- array(seq(0,n.data*n.x), c(n.data,n.x))
# #ubound
# for (i in 1:n.data){
#   df.x <- as.data.frame(results.MC[i, , ])
#   df.x.pbox <- as.data.frame(results.MC.pbox[i, , ])
#   for (j in 1:n.x){
#     cdf.x.lower.MC[i,j]<- sum(df.x$V3[df.x$V2<=x.values[j]]) #note reverse, see Ferson 2015
#     cdf.x.upper.MC[i,j]<- sum(df.x$V3[df.x$V1<=x.values[j]])
#     cdf.x.lower.MC.pbox[i,j]<- sum(df.x.pbox$V3[df.x.pbox$V2<=x.values[j]]) #note reverse, see Ferson 2015
#     cdf.x.upper.MC.pbox[i,j]<- sum(df.x.pbox$V3[df.x.pbox$V1<=x.values[j]])
#   }
# }
# #summary of cdf
# x.lower.summary.MC <- array(rep(0,n.x*6), c(n.x,6))
# x.upper.summary.MC <- array(rep(0,n.x*6), c(n.x,6))
# 
# x.lower.summary.MC.pbox <- array(rep(0,n.x*6), c(n.x,6))
# x.upper.summary.MC.pbox <- array(rep(0,n.x*6), c(n.x,6))
# 
# for(i in 1:n.x){
#   x.lower.summary.MC[i,] <- summary(cdf.x.lower.MC[,i])
#   x.upper.summary.MC[i,] <- summary(cdf.x.upper.MC[,i])
#   
#   x.lower.summary.MC.pbox[i,] <- summary(cdf.x.lower.MC.pbox[,i])
#   x.upper.summary.MC.pbox[i,] <- summary(cdf.x.upper.MC.pbox[,i])
# }
# 
# x.min.lower.MC.df <- as.data.frame(cbind(x.values,x.lower.summary.MC[,1]))
# x.Q1.lower.MC.df <- as.data.frame(cbind(x.values,x.lower.summary.MC[,2]))
# x.median.lower.MC.df <- as.data.frame(cbind(x.values,x.lower.summary.MC[,3]))
# x.mean.lower.MC.df <- as.data.frame(cbind(x.values,x.lower.summary.MC[,4]))
# x.Q3.lower.MC.df <- as.data.frame(cbind(x.values,x.lower.summary.MC[,5]))
# x.max.lower.MC.df <- as.data.frame(cbind(x.values,x.lower.summary.MC[,6]))
# 
# x.min.lower.MC.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.MC.pbox[,1]))
# x.Q1.lower.MC.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.MC.pbox[,2]))
# x.median.lower.MC.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.MC.pbox[,3]))
# x.mean.lower.MC.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.MC.pbox[,4]))
# x.Q3.lower.MC.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.MC.pbox[,5]))
# x.max.lower.MC.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.MC.pbox[,6]))
# 
# colnames(x.min.lower.MC.df) <- c("x","y")
# colnames(x.Q1.lower.MC.df) <- c("x","y")
# colnames(x.median.lower.MC.df) <- c("x","y")
# colnames(x.mean.lower.MC.df) <- c("x","y")
# colnames(x.Q3.lower.MC.df) <- c("x","y")
# colnames(x.max.lower.MC.df) <- c("x","y")
# 
# colnames(x.min.lower.MC.df.pbox) <- c("x","y")
# colnames(x.Q1.lower.MC.df.pbox) <- c("x","y")
# colnames(x.median.lower.MC.df.pbox) <- c("x","y")
# colnames(x.mean.lower.MC.df.pbox) <- c("x","y")
# colnames(x.Q3.lower.MC.df.pbox) <- c("x","y")
# colnames(x.max.lower.MC.df.pbox) <- c("x","y")
# 
# x.min.upper.MC.df <- as.data.frame(cbind(x.values,x.upper.summary.MC[,1]))
# x.Q1.upper.MC.df <- as.data.frame(cbind(x.values,x.upper.summary.MC[,2]))
# x.median.upper.MC.df <- as.data.frame(cbind(x.values,x.upper.summary.MC[,3]))
# x.mean.upper.MC.df <- as.data.frame(cbind(x.values,x.upper.summary.MC[,4]))
# x.Q3.upper.MC.df <- as.data.frame(cbind(x.values,x.upper.summary.MC[,5]))
# x.max.upper.MC.df <- as.data.frame(cbind(x.values,x.upper.summary.MC[,6]))
# 
# x.min.upper.MC.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.MC.pbox[,1]))
# x.Q1.upper.MC.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.MC.pbox[,2]))
# x.median.upper.MC.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.MC.pbox[,3]))
# x.mean.upper.MC.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.MC.pbox[,4]))
# x.Q3.upper.MC.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.MC.pbox[,5]))
# x.max.upper.MC.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.MC.pbox[,6]))
# 
# colnames(x.min.upper.MC.df) <- c("x","y")
# colnames(x.Q1.upper.MC.df) <- c("x","y")
# colnames(x.median.upper.MC.df) <- c("x","y")
# colnames(x.mean.upper.MC.df) <- c("x","y")
# colnames(x.Q3.upper.MC.df) <- c("x","y")
# colnames(x.max.upper.MC.df) <- c("x","y")
# 
# colnames(x.min.upper.MC.df.pbox) <- c("x","y")
# colnames(x.Q1.upper.MC.df.pbox) <- c("x","y")
# colnames(x.median.upper.MC.df.pbox) <- c("x","y")
# colnames(x.mean.upper.MC.df.pbox) <- c("x","y")
# colnames(x.Q3.upper.MC.df.pbox) <- c("x","y")
# colnames(x.max.upper.MC.df.pbox) <- c("x","y")

#############################################################
#run with sp
sp.indicator <- 1 #use support point (1) or MC (0)

if (sp.indicator==1){
  params_pdf <- list_test.sp
  params_pdf.pbox <- list_test.sp.pbox
  
} else{
  params_pdf <- list_test
  params_pdf.pbox <- list_test.pbox
}
#params constant
params_c <- c12
#non-parallel run if needed
# ptm1 <- proc.time()
# result_list_outer <- lapply(params_pdf,fn_bbmod1_simp_par_inner, params_c, input_data,pbox.int,index.array)
# time1 <- (proc.time()-ptm1)
# ela1 <-time1[1]
##############################################################################
num.cores <- detectCores()
availableCores()
#access VM - rstudio docker through shell
cl <- makeClusterPSOCK(
  # Public IP number of EC2 instance
  workers=availableCores(),
  dryrun = FALSE,
  connectTimeout = 120,
  verbose=TRUE,
  outfile=NULL,
  rshlogfile=TRUE
)

plan(list(tweak(cluster, workers = cl), multiprocess))
set.seed(710)
ptm2 <- proc.time()
result_list_outer.sp.pbox <- future.apply::future_lapply(params_pdf.pbox,fn_bbmod1_simp_par_inner, params_c, input_data,pbox.int,index.array)
result_list_outer.sp <- future.apply::future_lapply(params_pdf,fn_bbmod1_simp, params_c)
time2 <- (proc.time()-ptm2)
ela2 <-time2[3]

#support point results
n.data <- n.sp
results.sp <- array(rep(0,n.data), c(n.data))
#results.sp <- array(rep(0,n.data*n.indices*3), c(n.data,n.indices,3))

results.sp.pbox <- array(rep(0,n.data*n.indices*3), c(n.data,n.indices,3))

for (i in 1:n.data){
  results.sp[i] <- result_list_outer.sp[[i]]
  z <- 0
  for (j in 1:3){
    for (k in 1:n.indices){
      z <- z+1
      #results.sp[i,k,j] <- result_list_outer.sp[[i]][[z]]
      results.sp.pbox[i,k,j] <- result_list_outer.sp.pbox[[i]][[z]]
    }
  }
}
results.sp <- abs(results.sp)
results.sp.pbox <- abs(results.sp.pbox)

#derive empirical cdf of outcomes
min.x <- 0
max.x <- 40
n.x <-1000
x.values <- seq(min.x,max.x,length.out=n.x)
cdf.x.upper.sp <- array(seq(0,n.data*n.x), c(n.data,n.x))
#cdf.x.lower.sp <- array(seq(0,n.data*n.x), c(n.data,n.x))
cdf.x.lower.sp <- array(seq(0,n.x), c(n.x))

cdf.x.upper.sp.pbox <- array(seq(0,n.data*n.x), c(n.data,n.x))
cdf.x.lower.sp.pbox <- array(seq(0,n.data*n.x), c(n.data,n.x))
#ubound
#this is for MC/sp run, no bounds
#second columns is for weight (1/n.data)
wt <- c(rep(1/n.data,n.data))
df.x <- as.data.frame(cbind(results.sp,wt)) 
colnames(df.x) <- c("V1","V2")

for (i in 1:n.data){
  #df.x <- as.data.frame(results.sp[i, , ])
  df.x.pbox <- as.data.frame(results.sp.pbox[i, , ])
  
  for (j in 1:n.x){
    #cdf.x.lower.sp[i,j]<- sum(df.x$V3[df.x$V2<=x.values[j]]) #note reverse, see Ferson 2015
    #cdf.x.upper.sp[i,j]<- sum(df.x$V3[df.x$V1<=x.values[j]])
    
    cdf.x.lower.sp.pbox[i,j]<- sum(df.x.pbox$V3[df.x.pbox$V2<=x.values[j]]) #note reverse, see Ferson 2015
    cdf.x.upper.sp.pbox[i,j]<- sum(df.x.pbox$V3[df.x.pbox$V1<=x.values[j]])
  }
}

for (j in 1:n.x){
  cdf.x.lower.sp[j]<- sum(df.x$V2[df.x$V1<=x.values[j]]) #note reverse, see Ferson 2015
}

#summary of cdf
#x.lower.summary.sp <- array(rep(0,n.x*6), c(n.x,6))
#x.upper.summary.sp <- array(rep(0,n.x*6), c(n.x,6))

x.lower.summary.sp.pbox <- array(rep(0,n.x*6), c(n.x,6))
x.upper.summary.sp.pbox <- array(rep(0,n.x*6), c(n.x,6))


for(i in 1:n.x){
  #x.lower.summary.sp[i,] <- summary(cdf.x.lower.sp[,i])
  #x.upper.summary.sp[i,] <- summary(cdf.x.upper.sp[,i])
  
  x.lower.summary.sp.pbox[i,] <- summary(cdf.x.lower.sp.pbox[,i])
  x.upper.summary.sp.pbox[i,] <- summary(cdf.x.upper.sp.pbox[,i])
}

#x.min.lower.sp.df <- as.data.frame(cbind(x.values,x.lower.summary.sp[,1]))
#x.Q1.lower.sp.df <- as.data.frame(cbind(x.values,x.lower.summary.sp[,2]))
#x.median.lower.sp.df <- as.data.frame(cbind(x.values,x.lower.summary.sp[,3]))
#x.mean.lower.sp.df <- as.data.frame(cbind(x.values,x.lower.summary.sp[,4]))
x.mean.lower.sp.df <- as.data.frame(cbind(x.values,cdf.x.lower.sp))
#x.Q3.lower.sp.df <- as.data.frame(cbind(x.values,x.lower.summary.sp[,5]))
#x.max.lower.sp.df <- as.data.frame(cbind(x.values,x.lower.summary.sp[,6]))

x.min.lower.sp.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.sp.pbox[,1]))
x.Q1.lower.sp.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.sp.pbox[,2]))
x.median.lower.sp.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.sp.pbox[,3]))
x.mean.lower.sp.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.sp.pbox[,4]))
x.Q3.lower.sp.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.sp.pbox[,5]))
x.max.lower.sp.df.pbox <- as.data.frame(cbind(x.values,x.lower.summary.sp.pbox[,6]))

#colnames(x.min.lower.sp.df) <- c("x","y")
#colnames(x.Q1.lower.sp.df) <- c("x","y")
#colnames(x.median.lower.sp.df) <- c("x","y")
colnames(x.mean.lower.sp.df) <- c("x","y")
#colnames(x.Q3.lower.sp.df) <- c("x","y")
#colnames(x.max.lower.sp.df) <- c("x","y")

colnames(x.min.lower.sp.df.pbox) <- c("x","y")
colnames(x.Q1.lower.sp.df.pbox) <- c("x","y")
colnames(x.median.lower.sp.df.pbox) <- c("x","y")
colnames(x.mean.lower.sp.df.pbox) <- c("x","y")
colnames(x.Q3.lower.sp.df.pbox) <- c("x","y")
colnames(x.max.lower.sp.df.pbox) <- c("x","y")

# x.min.upper.sp.df <- as.data.frame(cbind(x.values,x.upper.summary.sp[,1]))
# x.Q1.upper.sp.df <- as.data.frame(cbind(x.values,x.upper.summary.sp[,2]))
# x.median.upper.sp.df <- as.data.frame(cbind(x.values,x.upper.summary.sp[,3]))
# x.mean.upper.sp.df <- as.data.frame(cbind(x.values,x.upper.summary.sp[,4]))
# x.Q3.upper.sp.df <- as.data.frame(cbind(x.values,x.upper.summary.sp[,5]))
# x.max.upper.sp.df <- as.data.frame(cbind(x.values,x.upper.summary.sp[,6]))

x.min.upper.sp.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.sp.pbox[,1]))
x.Q1.upper.sp.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.sp.pbox[,2]))
x.median.upper.sp.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.sp.pbox[,3]))
x.mean.upper.sp.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.sp.pbox[,4]))
x.Q3.upper.sp.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.sp.pbox[,5]))
x.max.upper.sp.df.pbox <- as.data.frame(cbind(x.values,x.upper.summary.sp.pbox[,6]))

# colnames(x.min.upper.sp.df) <- c("x","y")
# colnames(x.Q1.upper.sp.df) <- c("x","y")
# colnames(x.median.upper.sp.df) <- c("x","y")
# colnames(x.mean.upper.sp.df) <- c("x","y")
# colnames(x.Q3.upper.sp.df) <- c("x","y")
# colnames(x.max.upper.sp.df) <- c("x","y")

colnames(x.min.upper.sp.df.pbox) <- c("x","y")
colnames(x.Q1.upper.sp.df.pbox) <- c("x","y")
colnames(x.median.upper.sp.df.pbox) <- c("x","y")
colnames(x.mean.upper.sp.df.pbox) <- c("x","y")
colnames(x.Q3.upper.sp.df.pbox) <- c("x","y")
colnames(x.max.upper.sp.df.pbox) <- c("x","y")

###########################################################
#plot of cdf
c.color <- c(LBF="red",UBF="red",Precise="blue")

currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  
plot. <- ggplot()+ 
  #geom_line(data = x.min.lower.df, aes(x=x, y=y, color = "min.lower"), alpha=0.5, size=1) +
  #geom_line(data = x.Q1.lower.df, aes(x=x, y=y, color = "Q1.lower"), alpha=0.5, size=1) +
  #geom_line(data = x.median.lower.df, aes(x=x, y=y, color = "median.lower"), alpha=0.5, size=1) +
  #geom_line(data = x.mean.lower.MC.df, aes(x=x, y=y, color = "mean.lower.MC"), alpha=0.5, size=1) +
  geom_line(data = x.mean.lower.sp.df, aes(x=x, y=y, color = "Precise"), alpha=0.5, size=1) +
  #geom_line(data = x.mean.lower.MC.df.pbox, aes(x=x, y=y, color = "mean.lower.MC.pbox"), alpha=0.5, size=1) +
  geom_line(data = x.mean.lower.sp.df.pbox, aes(x=x, y=y, color = "LBF"), alpha=0.5, size=1) +
  #geom_line(data = x.Q3.lower.df, aes(x=x, y=y, color = "Q3.lower"), alpha=0.5, size=1) +
  #geom_line(data = x.max.lower.df, aes(x=x, y=y, color = "max.lower"), alpha=0.5, size=1) +
  #geom_line(data = x.min.upper.df, aes(x=x, y=y, color = "min.upper"), alpha=0.5, size=1) +
  #geom_line(data = x.Q1.upper.df, aes(x=x, y=y, color = "Q1.upper"), alpha=0.5, size=1) +
  #geom_line(data = x.median.upper.df, aes(x=x, y=y, color = "median.upper"), alpha=0.5, size=1) +
  #geom_line(data = x.mean.upper.MC.df, aes(x=x, y=y, color = "mean.upper.MC"), alpha=0.5, size=1) +
  #geom_line(data = x.mean.upper.sp.df, aes(x=x, y=y, color = "mean.upper.sp"), alpha=0.5, size=1) +
  #geom_line(data = x.mean.upper.MC.df.pbox, aes(x=x, y=y, color = "mean.upper.MC.pbox"), alpha=0.5, size=1) +
  geom_line(data = x.mean.upper.sp.df.pbox, aes(x=x, y=y, color = "UBF"), alpha=0.5, size=1) +
  #geom_line(data = x.Q3.upper.df, aes(x=x, y=y, color = "Q3.upper"), alpha=0.5, size=1) +
  #geom_line(data = x.max.upper.df, aes(x=x, y=y, color = "max.upper"), alpha=0.5, size=1) +
  scale_x_continuous(name = "QoI",expand=c(0,0),limits=c(0,max.x)) +
  scale_y_continuous(name = "CDF",expand=c(0,0),limits=c(0, 1)) +
  scale_color_manual(values = c.color,
                     name = "CDF types")+
  #ggtitle("P-box of quantity of interest (QoI)" ) +
  theme(axis.line = element_line(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12,face="bold"),
    legend.title = element_text(color = "black", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(color = "black", size = 10)) +
  guides(color = guide_legend(title="CDF types"), linetype = guide_legend(title="CDF types"))
  ggsave(paste(currpath,"/result/pbox_comparison_std_npbox_5_sp_", n.sp , ".tiff",sep="") , units="in", width=6, height=4, dpi=300, compression = 'lzw')

  
##############################################################################
#end of code
