########################################################################################
########################################################################################
#p-box functions and their inverses
#started 31072020
#updated 29082020
#rowan iskandar
########################################################################################
########################################################################################
#part 1
#pbox functions (given different data)
#given min, max
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

#given min, max, mean, std
pbox_mean_std <- function(x,type, min, max, mean, std){
  xi1 <- mean - (std^2)/(max-mean)
  xi2 <- mean + std^2/(mean-min)
  #print(xi1)
  #print(xi2)
  if (type=="lower"){
    if (x <= xi1 ){
      p=0
    }
    else if (xi1 < x && x <= xi2 ){
      num <- std^2+(max-mean)*(x-mean)
      denom <- (max-min)*(x-min)
      p <- num/denom
    }
    else {
      num <- (x-mean)^2
      denom <- (x-mean)^2 + std^2
      p <- num/denom
    }
  }
  else {
    if (x <= xi1 ){
      num <- (std)^2
      denom <- (mean-x)^2 + std^2
      p <- num/denom
    }
    else if (xi1 < x && x <= xi2 ){
      num <- (max-mean)*(max-min+mean-x)-std^2
      denom <- (max-min)*(max-x)
      p <- num/denom
    }
    else {
      p=1
    }
  }
  return(p)
}

pbox_median_std <- function(x, type, min, max, median, mean, std){
  if (type=="upper"){ #
    min.median<-pbox_median(x,type,min,max,median)
    min.std<-pbox_mean_std(x,type,min,max,mean,std)
    return(min(min.median,min.std))
  }
  else{
    max.median<-pbox_median(x,type,min,max,median)
    max.std<-pbox_mean_std(x,type,min,max,mean,std)
    return(max(max.median,max.std))
  }
}
########################################################################################
#inverse p-box functions for simulating p-box
#inverse for mean
#new version as of 04112020 #incorrect!
# pbox_mean_inv <- function(u,type,min,max,mean){
#   #crit <- (max-mean)/(max-min)
#   crit <- mean
#   if (type=="lower"){
#     #need to double check
#     if(u<0){
#       #x=max
#       x=mean
#       #x=min
#     }
#     else if(0<=u && u<crit){ #changed from crit < u
#       x=(u*max-max+mean)/u
#     }
#     # else if(u==crit){
#     #   x=min
#     # }
#     else if(u>crit){
#       x=max
#     }
#   }
#   else{
#     if (u==1){
#       x=mean
#     }
#     else if (0<=u && u< crit){ #changed from 0<u
#       x=min
#     }
#     # else if(u==0){
#     #   #x=min #not sure, min <= x < mean will yield u=0
#     # }
#     else if(crit <=u && u<1){ #not sure about the condition
#       x=(u*min-mean)/(u-1)
#     }
#   }
#   return(x)
# }
#older (may 19 2021)
# pbox_mean_inv <- function(u,type,min,max,mean){
#   crit <- (max-mean)/(max-min)
#   print(crit)
#   if (type=="upper"){
#     #need to double check
#     if(u==1){
#       x=max
#       #x=mean
#     }
#     else if(crit<u && u<1){
#       x=(u*max-max+mean)/u
#     }
#     else if(u==crit){
#       x=min
#     }
#     else if(u<crit){
#       x=min
#     }
#   }
#   else{
#     if (u==1){
#       x=max
#     }
#     else if (0<u && u< crit){
#       x=(u*min-mean)/(u-1)
#     }
#     else if(u<=0){
#       x=min #not sure, min <= x < mean will yield u=0
#     }
#     else if(crit <=u && u<1){ #not sure about the condition
#       x=max
#     }
#   }
#   return(x)
# }

#newest version may 19 2021 correct version!!!!!
pbox_mean_inv <- function(u,type,min,max,mean){
  crit <- (max-mean)/(max-min)
  #print(crit)
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
    else if(u<=0){
      x=(u*min-mean)/(u-1) #not sure, min <= x < mean will yield u=0
    }
    else if(crit <=u && u<1){ #not sure about the condition
      x=max
    }
  }
  return(x)
}

#inverse for minmax
#newer version may 19 2021
pbox_minmax_inv <-function(u,type,min,max){
  if (type=="upper"){
    if(u<1){
      x=min
    }
    else{x=max}
  }
  else{
    if(u==0){
      x=min}
    else{
      x=max
    }
  }
  return (x)
}
#older version
# pbox_minmax_inv <-function(u,type,min,max){
#   if (type=="upper"){
#     x=min
#   }
#   else{
#     x=max
#   }
#   return (x)
# }

#inverse for median
pbox_median_inv <-function(u,type,min,max,median){
  if (type=="upper"){ 
    if (u<0.5){
      x=min
    }
    else if(u >= 0.5 && u<1){
      x=median
    }
    else{ #different from ver 1
      x=max
    }
  }
  else{
    if (u<0.5){
      x=median
    }
    else if(0.5 <= u){
      x=max
    }
    else{
      x=min
    }
  }
  return (x)
}

#newer version, May 19 2021
pbox_std_inv <- function(u,type, min, max, mean, std){
  F1 <- (std^2)/((mean-min)^2+std^2)
  #print(F1)
  F2 <- (max-mean)^2/((max-mean)^2 + std^2)
  #print(F2)
  if (type=="upper"){
    if (u < F1 ){
      x=min
    }
    else if ( F1 <=u  && u < F2 ){
      x=mean-std*sqrt((1-u)/u)
    }
    else if (u >= F2){ 
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
    else if (u> F2) {
      x=max
    }
  }
  return(x)
}

#older
# pbox_std_inv <- function(u,type, min, max, mean, std){
#   F1 <- (std^2)/((mean-min)^2+std^2)
#   print(F1)
#   F2 <- (max-mean)^2/((max-mean)^2 + std^2)
#   print(F2)
#   if (type=="upper"){
#     if (u < F1 ){
#       x=min
#     }
#     else if ( F1 <=u  && u < F2 ){
#       x=mean-std*sqrt((1-u)/u)
#     }
#     else {
#       num <- (max-mean)*(mean-min)-std^2
#       denom <- u*(max-min)+mean-max
#       x <- max - num/denom
#     }
#   }
#   else {
#     if (u <= F1 ){
#       num <- (max-mean)*(mean-min)-std^2
#       denom <- max-mean-u*(max-min)
#       x<- min + num/denom
#     }
#     else if ( F1 < u  && u <= F2 ){
#       x=mean+std*sqrt(u/(1-u))
#     }
#     else {
#       x=max
#     }
#   }
#   return(x)
# }

#newest 19 May 2021
pbox_median_std_inv <-function(u,type,min,max, median, mean, std){
  F1 <- (std^2)/((mean-min)^2+std^2)
  # print("---------------------------------")
  # print(F1)
  F2 <- (max-mean)^2/((max-mean)^2 + std^2)
  # print(F2)
  # print("---------------------------------")
  if (type=="upper"){ #
    if(u < F1){
      x=min
    }
    else if (F1 <=u && u< 0.5){
      x=mean-std*sqrt((1-u)/u)
    }
    else if (0.5 <= u && u< F2){
      x=median
    }
    else if ( F2 <= u && u <1){
      num <- (max-mean)*(mean-min)-std^2
      denom <- u*(max-min)+mean-max
      xtemp <- max - num/denom
      if (xtemp < median){
        x=median
      }
      else{x<-xtemp}
    }
    else if(u==1){
      x=max
    }
  }
  else{
    if(u<F1){
      num <- (max-mean)*(mean-min)-std^2
      denom <- max-mean-u*(max-min)
      xtemp<- min + num/denom
      if (xtemp >median){
        x<-median
      } else {
        x <-xtemp
      }
    }
    else if (F1 <= u && u <0.5){
      x<-median
    }
    else if(0.5 <= u && u <F2){
      x=mean+std*sqrt(u/(1-u))
    }
    else if (F2<=u){
      x=max
    }
  }
  return(x)
}
# 
# #older 19 May 2021
# pbox_median_std_inv_2 <-function(u,type,min,max, median, mean, std){
#   if (type=="upper"){ #
#     min.median<-pbox_median_inv(u,type,min,max,median)
#     min.std<-pbox_std_inv(u,type,min,max,mean,std)
#     return(min(min.median,min.std))
#   }
#   else{
#     max.median<-pbox_median_inv(u,type,min,max,median)
#     max.std<-pbox_std_inv(u,type,min,max,mean,std)
#     return(max(max.median,max.std))
#   }
# }



