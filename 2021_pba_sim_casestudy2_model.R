#R file for case study 2
#Probability bound analysis: A novel approach for quantifying parameter uncertainty in decision-analytic modeling and cost-effectiveness analysis
#author: Rowan Iskandar
#model file (used in 2021_pba_sim_casestudy2_analysis)
#based on
#Dong, Hengjin, and Martin Buxton. "Early assessment of the likely cost-effectiveness of a new technology: a Markov model with probabilistic sensitivity analysis of computer-assisted total knee replacement." International Journal of Technology Assessment in Health Care 22.2 (2006): 191-202.
#last updated: 13 August 2021

source("pbox.R") #make sure this file is at the same folder level

discount <- function(t,disc){
  weight <- 1/(1+disc/12)^t
  return(weight)
}

TKR_model_average <- function(
                      rr.CAS,
                      c.CAS,
                      p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                      p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                      p.tkr_minor.tkr_normal,p.tkr_minor.tkr_serious, p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                      p.tkr_serious.tkr_minor,p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                      p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                      p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,p.revision_simple.other_trt,
                      p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                      p.other_trt.tkr_minor,p.other_trt.tkr_serious,p.other_trt.tkr_normal_after_revision,p.other_trt.revision_simple,
                      p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                      c.tkr,
                      c.tkr_normal,
                      c.tkr_minor,
                      c.tkr_serious,
                      c.tkr_normal_after_revision,
                      c.revision_simple,
                      c.revision_complex,
                      c.other_trt,
                      c.death,
                      u.tkr,
                      u.tkr_normal,
                      u.tkr_minor,
                      u.tkr_serious,
                      u.tkr_normal_after_revision,
                      u.revision_simple,
                      u.revision_complex,
                      u.other_trt,
                      disc ,
                      num.states,
                      time.horizon,
                      init,
                      n,
                      WTP)
  {
  
  #1 tkr_init
  #2 tkr_normal
  #3 tkr_minor
  #4 tkr_serious
  #5 tkr_normal_after_revision
  #6 revision_simple
  #7 revision_complex
  #8 other_trt
  #9 death
  
  mat.p <- matrix(rep(0,num.states*num.states),nrow=num.states,ncol=num.states)
  
  mat.p.CAS <- matrix(rep(0,num.states*num.states),nrow=num.states,ncol=num.states)
  
  p.tkr_init.tkr_serious.CAS <- (rr.CAS)*p.tkr_init.tkr_serious
  p.tkr_normal.tkr_serious.CAS <- (rr.CAS)*p.tkr_normal.tkr_serious
  p.tkr_minor.tkr_serious.CAS <- (rr.CAS)*p.tkr_minor.tkr_serious
  p.revision_simple.tkr_serious.CAS <- (rr.CAS)*p.revision_simple.tkr_serious
  p.revision_complex.tkr_serious.CAS <- (rr.CAS)*p.revision_complex.tkr_serious
  p.other_trt.tkr_serious.CAS <- (rr.CAS)*p.other_trt.tkr_serious
  
  mat.p[1,2] <- p.tkr_init.tkr_normal
  mat.p[1,3] <- p.tkr_init.tkr_minor
  mat.p[1,4] <- p.tkr_init.tkr_serious

  mat.p[2,3] <- p.tkr_normal.tkr_minor
  mat.p[2,4] <- p.tkr_normal.tkr_serious
  mat.p[2,9] <- p.death + p.tkr_normal.death
  mat.p[2,2] <- 1-sum(mat.p[2, -c(2)]) #p.tkr_normal.tkr_normal
  
  
  mat.p[3,2] <- p.tkr_minor.tkr_normal
  mat.p[3,4] <- p.tkr_minor.tkr_serious
  mat.p[3,6] <- p.tkr_minor.revision_simple
  mat.p[3,8] <- p.tkr_minor.other_trt
  mat.p[3,9] <- p.death + p.tkr_normal.death
  mat.p[3,3] <- 1-sum(mat.p[3, -c(3)]) #p.tkr_minor.tkr_minor
  
  mat.p[4,3] <- p.tkr_serious.tkr_minor
  mat.p[4,4] <- p.tkr_serious.other_trt
  mat.p[4,6] <- p.tkr_serious.revision_simple
  mat.p[4,7] <- p.tkr_serious.revision_complex
  # mat.p[4,8] <- rowSums(mat.p[4, -4])
  mat.p[4,9] <- p.death + p.tkr_normal.death
  
  mat.p[5,6] <- p.tkr_normal_after_revision.revision_simple
  mat.p[5,7] <- p.tkr_normal_after_revision.revision_complex
  mat.p[5,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[5,5] <- 1-sum(mat.p[5, -c(5)])
  
  
  mat.p[6,3] <- p.revision_simple.tkr_minor
  mat.p[6,4] <- p.revision_simple.tkr_serious
  # mat.p[6,6] <- rowSums(mat.p[6, -6])
  mat.p[6,8] <- p.revision_simple.other_trt
  mat.p[6,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[6,5] <- 1-sum(mat.p[6, -c(5)])
  
  mat.p[7,4] <- p.revision_complex.tkr_serious
  #mat.p[7,7] <- rowSums(mat.p[7, -7])
  mat.p[7,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[7,5] <- 1-sum(mat.p[7, -c(5)])
  
  mat.p[8,3] <- p.other_trt.tkr_minor
  mat.p[8,4] <- p.other_trt.tkr_serious
  mat.p[8,7] <- p.other_trt.revision_simple
  #mat.p[8,8] <- p.other_trt.revision_simple
  mat.p[8,9] <- p.death + p.tkr_normal.death
  mat.p[8,5] <- 1-sum(mat.p[8, -c(5)])
  
  mat.p.CAS[1,4] <- p.tkr_init.tkr_serious.CAS
  p.residual <- 1-(p.tkr_init.tkr_serious.CAS+p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)
  mat.p.CAS[1,2] <- p.tkr_init.tkr_normal + p.tkr_init.tkr_normal/(p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)*p.residual
  mat.p.CAS[1,3] <- p.tkr_init.tkr_minor + p.tkr_init.tkr_minor/(p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)*p.residual

  mat.p.CAS[2,3] <- p.tkr_normal.tkr_minor
  mat.p.CAS[2,4] <- p.tkr_normal.tkr_serious.CAS
  mat.p.CAS[2,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[2,2] <- 1-sum(mat.p.CAS[2, -c(2)]) #p.tkr_normal.tkr_normal
  
  
  mat.p.CAS[3,2] <- p.tkr_minor.tkr_normal
  mat.p.CAS[3,4] <- p.tkr_minor.tkr_serious.CAS
  mat.p.CAS[3,6] <- p.tkr_minor.revision_simple
  mat.p.CAS[3,8] <- p.tkr_minor.other_trt
  mat.p.CAS[3,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[3,3] <- 1-sum(mat.p.CAS[3, -c(3)]) #p.tkr_minor.tkr_minor
  
  mat.p.CAS[4,3] <- p.tkr_serious.tkr_minor
  mat.p.CAS[4,4] <- p.tkr_serious.other_trt
  mat.p.CAS[4,6] <- p.tkr_serious.revision_simple
  mat.p.CAS[4,7] <- p.tkr_serious.revision_complex
  # mat.p.CAS[4,8] <- rowSums(mat.p.CAS[4, -4])
  mat.p.CAS[4,9] <- p.death + p.tkr_normal.death
  
  mat.p.CAS[5,6] <- p.tkr_normal_after_revision.revision_simple
  mat.p.CAS[5,7] <- p.tkr_normal_after_revision.revision_complex
  mat.p.CAS[5,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[5,5] <- 1-sum(mat.p.CAS[5, -c(5)])
  
  mat.p.CAS[6,3] <- p.revision_simple.tkr_minor
  mat.p.CAS[6,4] <- p.revision_simple.tkr_serious.CAS
  # mat.p.CAS[6,6] <- rowSums(mat.p.CAS[6, -6])
  mat.p.CAS[6,8] <- p.revision_simple.other_trt
  mat.p.CAS[6,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[6,5] <- 1-sum(mat.p.CAS[6, -c(5)])
  
  mat.p.CAS[7,4] <- p.revision_complex.tkr_serious.CAS
  #mat.p.CAS[7,7] <- rowSums(mat.p.CAS[7, -7])
  mat.p.CAS[7,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[7,5] <- 1-sum(mat.p.CAS[7,-c(5)])
  
  mat.p.CAS[8,3] <- p.other_trt.tkr_minor
  mat.p.CAS[8,4] <- p.other_trt.tkr_serious.CAS
  mat.p.CAS[8,7] <- p.other_trt.revision_simple
  #mat.p.CAS[8,8] <- p.other_trt.revision_simple
  mat.p.CAS[8,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[8,5] <- 1-sum(mat.p.CAS[8, -c(5)])
  
  #check matrix - rowsum == 1
  mat <- rowSums(mat.p.CAS)
  print(mat)

  
  #utilities
  utility <- c(u.tkr,
               u.tkr_normal,
               u.tkr_minor,
               u.tkr_serious,
               u.tkr_normal_after_revision,
               u.revision_simple,
               u.revision_complex,
               u.other_trt,
               0)
  utility <- as.matrix(utility)
  utility <-t(matrix(rep(utility,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  cost <- c(c.tkr,
          c.tkr_normal,
          c.tkr_minor,
          c.tkr_serious,
          c.tkr_normal_after_revision,
          c.revision_simple,
          c.revision_complex,
          c.other_trt,
          c.death)
  cost <- as.matrix(cost)
  cost <-t(matrix(rep(cost,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  
  cost.CAS <- c(c.tkr+c.CAS,
           c.tkr_normal,
           c.tkr_minor,
           c.tkr_serious,
           c.tkr_normal_after_revision,
           c.revision_simple+c.CAS,
           c.revision_complex+c.CAS,
           c.other_trt,
           c.death)
  cost.CAS <- as.matrix(cost.CAS)
  cost.CAS <-t(matrix(rep(cost.CAS,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  
  #discount weight
  disc.weight <- lapply(c(0:time.horizon),discount,disc=disc)
  disc.weight <- as.matrix(unlist(disc.weight))
  disc.weight <-matrix(rep(disc.weight,num.states),nrow=time.horizon+1,ncol=num.states)
  
  #running
  cohort <- matrix(rep(NA,num.states*(time.horizon+1)),nrow=time.horizon+1, ncol=num.states)
  cohort.CAS <- matrix(rep(NA,num.states*(time.horizon+1)),nrow=time.horizon+1, ncol=num.states)
  
  cohort[1,] <- init 
  cohort.CAS[1,] <- init 

  for(i in 2:(time.horizon+1)){
    cohort[i,] <-cohort[i-1,]%*%mat.p
    cohort.CAS[i,] <- cohort.CAS[i-1,]%*%mat.p.CAS
  }

  #calculate QoI
  LE <- sum(cohort[,1:8])/n
  mean.QALY <- sum((cohort*utility)*disc.weight)/n
  mean.cost <- sum((cohort*cost)*disc.weight)/n
  # QALY <- cohort%*%utility*disc.weight
  # cost <- cohort%*%cost*disc.weight
  QALY <- cohort*utility*disc.weight
  cost <- cohort*cost*disc.weight
  minor <- cohort[,3]
  serious <- cohort[,4]
  simple <- cohort[,6]
  complex <- cohort[,7]
  death <- cohort[,9]
  
  LE.CAS <- sum(cohort.CAS[,1:8])/n
  mean.QALY.CAS <- sum((cohort.CAS*utility)*disc.weight)/n
  QALY.CAS <- (cohort.CAS*utility)*disc.weight
  mean.cost.CAS <- sum((cohort.CAS*cost.CAS)*disc.weight)/n
  QALY.CAS <- cohort.CAS*utility*disc.weight
  cost.CAS <- cohort.CAS*cost.CAS*disc.weight
  
  minor.CAS <- cohort.CAS[,3]
  serious.CAS <- cohort.CAS[,4]
  simple.CAS <- cohort.CAS[,6]
  complex.CAS <- cohort.CAS[,7]
  death.CAS <- cohort.CAS[,9]
  
  NMB <- mean.QALY*WTP - mean.cost
  NMB.CAS <- mean.QALY.CAS*WTP - mean.cost.CAS
  INMB <- NMB.CAS-NMB

  # if(objective.index==1){
  #   return(INMB)
  # }
  # else{
  #   return(-INMB)
  # }
  
  output <- list(INMB, LE, cost, mean.cost, QALY, mean.QALY, LE.CAS, cost.CAS, mean.cost.CAS, QALY.CAS, mean.QALY.CAS,
                 minor, serious, simple, complex, minor.CAS, serious.CAS, simple.CAS, complex.CAS,death,death.CAS, QALY, QALY.CAS)
  return(output)
}

# p.tkr_minor.tkr_serious
# p.tkr_serious.tkr_minor
# p.revision_simple.other_trt
# p.other_trt.tkr_minor
# p.other_trt.tkr_serious
# p.other_trt.revision_simple


TKR_model <- function(list_params,
                      objective.index,
                      #rr.CAS,
                      c.CAS,
                      p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                      p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                      p.tkr_minor.tkr_normal,
                      #p.tkr_minor.tkr_serious, 
                      p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                      #p.tkr_serious.tkr_minor,
                      p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                      p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                      p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,
                      #p.revision_simple.other_trt,
                      p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                      #p.other_trt.tkr_minor,p.other_trt.tkr_serious,
                      p.other_trt.tkr_normal_after_revision,
                      #p.other_trt.revision_simple,
                      p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                      c.tkr,
                      c.tkr_normal,
                      c.tkr_minor,
                      c.tkr_serious,
                      c.tkr_normal_after_revision,
                      c.revision_simple,
                      c.revision_complex,
                      c.other_trt,
                      c.death,
                      u.tkr,
                      u.tkr_normal,
                      u.tkr_minor,
                      u.tkr_serious,
                      u.tkr_normal_after_revision,
                      u.revision_simple,
                      u.revision_complex,
                      u.other_trt,
                      disc ,
                      num.states,
                      time.horizon,
                      init,
                      n,
                      WTP)
{
  
  #1 tkr_init
  #2 tkr_normal
  #3 tkr_minor
  #4 tkr_serious
  #5 tkr_normal_after_revision
  #6 revision_simple
  #7 revision_complex
  #8 other_trt
  #9 death
  
  rr.CAS <- list_params[1]
  p.tkr_minor.tkr_serious <- list_params[2]
  p.tkr_serious.tkr_minor <- list_params[3]
  p.revision_simple.other_trt <- list_params[4]
  p.other_trt.tkr_minor <- list_params[5]
  p.other_trt.tkr_serious <- list_params[6]
  p.other_trt.revision_simple <- list_params[7]
  # print("---------")
  # print("---------")
  # print(rr.CAS)
  # print("---------")
  # print(p.tkr_init.tkr_serious)
  # print("---------")
  
  mat.p <- matrix(rep(0,num.states*num.states),nrow=num.states,ncol=num.states)
  
  mat.p.CAS <- matrix(rep(0,num.states*num.states),nrow=num.states,ncol=num.states)
  
  p.tkr_init.tkr_serious.CAS <- (1-rr.CAS)*p.tkr_init.tkr_serious
  p.tkr_normal.tkr_serious.CAS <- (1-rr.CAS)*p.tkr_normal.tkr_serious
  p.tkr_minor.tkr_serious.CAS <- (1-rr.CAS)*p.tkr_minor.tkr_serious
  p.revision_simple.tkr_serious.CAS <- (1-rr.CAS)*p.revision_simple.tkr_serious
  p.revision_complex.tkr_serious.CAS <- (1-rr.CAS)*p.revision_complex.tkr_serious
  p.other_trt.tkr_serious.CAS <- (1-rr.CAS)*p.other_trt.tkr_serious
  
  mat.p[1,2] <- p.tkr_init.tkr_normal
  mat.p[1,3] <- p.tkr_init.tkr_minor
  mat.p[1,4] <- p.tkr_init.tkr_serious
  
  mat.p[2,3] <- p.tkr_normal.tkr_minor
  mat.p[2,4] <- p.tkr_normal.tkr_serious
  mat.p[2,9] <- p.death + p.tkr_normal.death
  mat.p[2,2] <- 1-sum(mat.p[2, -c(2)]) #p.tkr_normal.tkr_normal
  
  
  mat.p[3,2] <- p.tkr_minor.tkr_normal
  mat.p[3,4] <- p.tkr_minor.tkr_serious
  mat.p[3,6] <- p.tkr_minor.revision_simple
  mat.p[3,8] <- p.tkr_minor.other_trt
  mat.p[3,9] <- p.death + p.tkr_normal.death
  mat.p[3,3] <- 1-sum(mat.p[3, -c(3)]) #p.tkr_minor.tkr_minor
  
  mat.p[4,3] <- p.tkr_serious.tkr_minor
  mat.p[4,4] <- p.tkr_serious.other_trt
  mat.p[4,6] <- p.tkr_serious.revision_simple
  mat.p[4,7] <- p.tkr_serious.revision_complex
  # mat.p[4,8] <- rowSums(mat.p[4, -4])
  mat.p[4,9] <- p.death + p.tkr_normal.death
  
  mat.p[5,6] <- p.tkr_normal_after_revision.revision_simple
  mat.p[5,7] <- p.tkr_normal_after_revision.revision_complex
  mat.p[5,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[5,5] <- 1-sum(mat.p[5, -c(5)])
  
  
  mat.p[6,3] <- p.revision_simple.tkr_minor
  mat.p[6,4] <- p.revision_simple.tkr_serious
  # mat.p[6,6] <- rowSums(mat.p[6, -6])
  mat.p[6,8] <- p.revision_simple.other_trt
  mat.p[6,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[6,5] <- 1-sum(mat.p[6, -c(5)])
  
  mat.p[7,4] <- p.revision_complex.tkr_serious
  #mat.p[7,7] <- rowSums(mat.p[7, -7])
  mat.p[7,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[7,5] <- 1-sum(mat.p[7, -c(5)])
  
  mat.p[8,3] <- p.other_trt.tkr_minor
  mat.p[8,4] <- p.other_trt.tkr_serious
  mat.p[8,7] <- p.other_trt.revision_simple
  #mat.p[8,8] <- p.other_trt.revision_simple
  mat.p[8,9] <- p.death + p.tkr_normal.death
  mat.p[8,5] <- 1-sum(mat.p[8, -c(5)])
  
  mat.p.CAS[1,4] <- p.tkr_init.tkr_serious.CAS
  p.residual <- 1-(p.tkr_init.tkr_serious.CAS+p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)
  mat.p.CAS[1,2] <- p.tkr_init.tkr_normal + p.tkr_init.tkr_normal/(p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)*p.residual
  mat.p.CAS[1,3] <- p.tkr_init.tkr_minor + p.tkr_init.tkr_minor/(p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)*p.residual
  
  mat.p.CAS[2,3] <- p.tkr_normal.tkr_minor
  mat.p.CAS[2,4] <- p.tkr_normal.tkr_serious.CAS
  mat.p.CAS[2,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[2,2] <- 1-sum(mat.p.CAS[2, -c(2)]) #p.tkr_normal.tkr_normal
  
  
  mat.p.CAS[3,2] <- p.tkr_minor.tkr_normal
  mat.p.CAS[3,4] <- p.tkr_minor.tkr_serious.CAS
  mat.p.CAS[3,6] <- p.tkr_minor.revision_simple
  mat.p.CAS[3,8] <- p.tkr_minor.other_trt
  mat.p.CAS[3,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[3,3] <- 1-sum(mat.p.CAS[3, -c(3)]) #p.tkr_minor.tkr_minor
  
  mat.p.CAS[4,3] <- p.tkr_serious.tkr_minor
  mat.p.CAS[4,4] <- p.tkr_serious.other_trt
  mat.p.CAS[4,6] <- p.tkr_serious.revision_simple
  mat.p.CAS[4,7] <- p.tkr_serious.revision_complex
  # mat.p.CAS[4,8] <- rowSums(mat.p.CAS[4, -4])
  mat.p.CAS[4,9] <- p.death + p.tkr_normal.death
  
  mat.p.CAS[5,6] <- p.tkr_normal_after_revision.revision_simple
  mat.p.CAS[5,7] <- p.tkr_normal_after_revision.revision_complex
  mat.p.CAS[5,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[5,5] <- 1-sum(mat.p.CAS[5, -c(5)])
  
  mat.p.CAS[6,3] <- p.revision_simple.tkr_minor
  mat.p.CAS[6,4] <- p.revision_simple.tkr_serious.CAS
  # mat.p.CAS[6,6] <- rowSums(mat.p.CAS[6, -6])
  mat.p.CAS[6,8] <- p.revision_simple.other_trt
  mat.p.CAS[6,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[6,5] <- 1-sum(mat.p.CAS[6, -c(5)])
  
  mat.p.CAS[7,4] <- p.revision_complex.tkr_serious.CAS
  #mat.p.CAS[7,7] <- rowSums(mat.p.CAS[7, -7])
  mat.p.CAS[7,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[7,5] <- 1-sum(mat.p.CAS[7,-c(5)])
  
  mat.p.CAS[8,3] <- p.other_trt.tkr_minor
  mat.p.CAS[8,4] <- p.other_trt.tkr_serious.CAS
  mat.p.CAS[8,7] <- p.other_trt.revision_simple
  #mat.p.CAS[8,8] <- p.other_trt.revision_simple
  mat.p.CAS[8,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[8,5] <- 1-sum(mat.p.CAS[8, -c(5)])
  
  #check matrix - rowsum == 1
  # mat <- rowSums(mat.p.CAS)
  # print(mat)
  # 
  # 
  #utilities
  utility <- c(u.tkr,
               u.tkr_normal,
               u.tkr_minor,
               u.tkr_serious,
               u.tkr_normal_after_revision,
               u.revision_simple,
               u.revision_complex,
               u.other_trt,
               0)
  utility <- as.matrix(utility)
  utility <-t(matrix(rep(utility,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  cost <- c(c.tkr,
            c.tkr_normal,
            c.tkr_minor,
            c.tkr_serious,
            c.tkr_normal_after_revision,
            c.revision_simple,
            c.revision_complex,
            c.other_trt,
            c.death)
  cost <- as.matrix(cost)
  cost <-t(matrix(rep(cost,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  
  cost.CAS <- c(c.tkr+c.CAS,
                c.tkr_normal,
                c.tkr_minor,
                c.tkr_serious,
                c.tkr_normal_after_revision,
                c.revision_simple+c.CAS,
                c.revision_complex+c.CAS,
                c.other_trt,
                c.death)
  cost.CAS <- as.matrix(cost.CAS)
  cost.CAS <-t(matrix(rep(cost.CAS,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  
  #discount weight
  disc.weight <- lapply(c(0:time.horizon),discount,disc=disc)
  disc.weight <- as.matrix(unlist(disc.weight))
  disc.weight <-matrix(rep(disc.weight,num.states),nrow=time.horizon+1,ncol=num.states)
  
  #running
  cohort <- matrix(rep(NA,num.states*(time.horizon+1)),nrow=time.horizon+1, ncol=num.states)
  cohort.CAS <- matrix(rep(NA,num.states*(time.horizon+1)),nrow=time.horizon+1, ncol=num.states)
  
  cohort[1,] <- init 
  cohort.CAS[1,] <- init 
  
  for(i in 2:(time.horizon+1)){
    cohort[i,] <-cohort[i-1,]%*%mat.p
    cohort.CAS[i,] <- cohort.CAS[i-1,]%*%mat.p.CAS
  }
  
  #calculate QoI
  LE <- sum(cohort[,1:8])/n
  mean.QALY <- sum((cohort*utility)*disc.weight)/n
  mean.cost <- sum((cohort*cost)*disc.weight)/n
  # QALY <- cohort%*%utility*disc.weight
  # cost <- cohort%*%cost*disc.weight
  QALY <- cohort*utility*disc.weight
  cost <- cohort*cost*disc.weight
  
  LE.CAS <- sum(cohort.CAS[,1:8])/n
  mean.QALY.CAS <- sum((cohort.CAS*utility)*disc.weight)/n
  mean.cost.CAS <- sum((cohort.CAS*cost.CAS)*disc.weight)/n
  QALY.CAS <- cohort.CAS*utility*disc.weight
  cost.CAS <- cohort.CAS*cost.CAS*disc.weight
  
  NMB <- mean.QALY*WTP - mean.cost
  NMB.CAS <- mean.QALY.CAS*WTP - mean.cost.CAS
  INMB <- NMB.CAS-NMB
  
  if(objective.index==1){
    return(INMB)
  }
  else{
    return(-INMB)
  }
  
  # output <- list(INMB, LE, cost, mean.cost, QALY, mean.QALY, LE.CAS, cost.CAS, mean.cost.CAS, QALY.CAS, mean.QALY.CAS)
  # print(cost)
  # return(output)
}


TKR_model_MC <- function(list_params,
                      #objective.index,
                      #rr.CAS,
                      c.CAS,
                      p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                      p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                      p.tkr_minor.tkr_normal,
                      #p.tkr_minor.tkr_serious, 
                      p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                      #p.tkr_serious.tkr_minor,
                      p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                      p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                      p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,
                      #p.revision_simple.other_trt,
                      p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                      #p.other_trt.tkr_minor,p.other_trt.tkr_serious,
                      p.other_trt.tkr_normal_after_revision,
                      #p.other_trt.revision_simple,
                      p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                      c.tkr,
                      c.tkr_normal,
                      c.tkr_minor,
                      c.tkr_serious,
                      c.tkr_normal_after_revision,
                      c.revision_simple,
                      c.revision_complex,
                      c.other_trt,
                      c.death,
                      u.tkr,
                      u.tkr_normal,
                      u.tkr_minor,
                      u.tkr_serious,
                      u.tkr_normal_after_revision,
                      u.revision_simple,
                      u.revision_complex,
                      u.other_trt,
                      disc ,
                      num.states,
                      time.horizon,
                      init,
                      n,
                      WTP)
{
  
  #1 tkr_init
  #2 tkr_normal
  #3 tkr_minor
  #4 tkr_serious
  #5 tkr_normal_after_revision
  #6 revision_simple
  #7 revision_complex
  #8 other_trt
  #9 death
  
  rr.CAS <- list_params[[1]]
  p.tkr_minor.tkr_serious <- list_params[[2]]
  p.tkr_serious.tkr_minor <- list_params[[3]]
  p.revision_simple.other_trt <- list_params[[4]]
  p.other_trt.tkr_minor <- list_params[[5]]
  p.other_trt.tkr_serious <- list_params[[6]]
  p.other_trt.revision_simple <- list_params[[7]]
  
  
  mat.p <- matrix(rep(0,num.states*num.states),nrow=num.states,ncol=num.states)
  
  mat.p.CAS <- matrix(rep(0,num.states*num.states),nrow=num.states,ncol=num.states)
  
  p.tkr_init.tkr_serious.CAS <- (1-rr.CAS)*p.tkr_init.tkr_serious
  p.tkr_normal.tkr_serious.CAS <- (1-rr.CAS)*p.tkr_normal.tkr_serious
  p.tkr_minor.tkr_serious.CAS <- (1-rr.CAS)*p.tkr_minor.tkr_serious
  p.revision_simple.tkr_serious.CAS <- (1-rr.CAS)*p.revision_simple.tkr_serious
  p.revision_complex.tkr_serious.CAS <- (1-rr.CAS)*p.revision_complex.tkr_serious
  p.other_trt.tkr_serious.CAS <- (1-rr.CAS)*p.other_trt.tkr_serious
  
  mat.p[1,2] <- p.tkr_init.tkr_normal
  mat.p[1,3] <- p.tkr_init.tkr_minor
  mat.p[1,4] <- p.tkr_init.tkr_serious
  
  mat.p[2,3] <- p.tkr_normal.tkr_minor
  mat.p[2,4] <- p.tkr_normal.tkr_serious
  mat.p[2,9] <- p.death + p.tkr_normal.death
  mat.p[2,2] <- 1-sum(mat.p[2, -c(2)]) #p.tkr_normal.tkr_normal
  
  
  mat.p[3,2] <- p.tkr_minor.tkr_normal
  mat.p[3,4] <- p.tkr_minor.tkr_serious
  mat.p[3,6] <- p.tkr_minor.revision_simple
  mat.p[3,8] <- p.tkr_minor.other_trt
  mat.p[3,9] <- p.death + p.tkr_normal.death
  mat.p[3,3] <- 1-sum(mat.p[3, -c(3)]) #p.tkr_minor.tkr_minor
  
  mat.p[4,3] <- p.tkr_serious.tkr_minor
  mat.p[4,4] <- p.tkr_serious.other_trt
  mat.p[4,6] <- p.tkr_serious.revision_simple
  mat.p[4,7] <- p.tkr_serious.revision_complex
  # mat.p[4,8] <- rowSums(mat.p[4, -4])
  mat.p[4,9] <- p.death + p.tkr_normal.death
  
  mat.p[5,6] <- p.tkr_normal_after_revision.revision_simple
  mat.p[5,7] <- p.tkr_normal_after_revision.revision_complex
  mat.p[5,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[5,5] <- 1-sum(mat.p[5, -c(5)])
  
  
  mat.p[6,3] <- p.revision_simple.tkr_minor
  mat.p[6,4] <- p.revision_simple.tkr_serious
  # mat.p[6,6] <- rowSums(mat.p[6, -6])
  mat.p[6,8] <- p.revision_simple.other_trt
  mat.p[6,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[6,5] <- 1-sum(mat.p[6, -c(5)])
  
  mat.p[7,4] <- p.revision_complex.tkr_serious
  #mat.p[7,7] <- rowSums(mat.p[7, -7])
  mat.p[7,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p[7,5] <- 1-sum(mat.p[7, -c(5)])
  
  mat.p[8,3] <- p.other_trt.tkr_minor
  mat.p[8,4] <- p.other_trt.tkr_serious
  mat.p[8,7] <- p.other_trt.revision_simple
  #mat.p[8,8] <- p.other_trt.revision_simple
  mat.p[8,9] <- p.death + p.tkr_normal.death
  mat.p[8,5] <- 1-sum(mat.p[8, -c(5)])
  
  mat.p.CAS[1,4] <- p.tkr_init.tkr_serious.CAS
  p.residual <- 1-(p.tkr_init.tkr_serious.CAS+p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)
  mat.p.CAS[1,2] <- p.tkr_init.tkr_normal + p.tkr_init.tkr_normal/(p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)*p.residual
  mat.p.CAS[1,3] <- p.tkr_init.tkr_minor + p.tkr_init.tkr_minor/(p.tkr_init.tkr_normal+p.tkr_init.tkr_minor)*p.residual
  
  mat.p.CAS[2,3] <- p.tkr_normal.tkr_minor
  mat.p.CAS[2,4] <- p.tkr_normal.tkr_serious.CAS
  mat.p.CAS[2,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[2,2] <- 1-sum(mat.p.CAS[2, -c(2)]) #p.tkr_normal.tkr_normal
  
  
  mat.p.CAS[3,2] <- p.tkr_minor.tkr_normal
  mat.p.CAS[3,4] <- p.tkr_minor.tkr_serious.CAS
  mat.p.CAS[3,6] <- p.tkr_minor.revision_simple
  mat.p.CAS[3,8] <- p.tkr_minor.other_trt
  mat.p.CAS[3,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[3,3] <- 1-sum(mat.p.CAS[3, -c(3)]) #p.tkr_minor.tkr_minor
  
  mat.p.CAS[4,3] <- p.tkr_serious.tkr_minor
  mat.p.CAS[4,4] <- p.tkr_serious.other_trt
  mat.p.CAS[4,6] <- p.tkr_serious.revision_simple
  mat.p.CAS[4,7] <- p.tkr_serious.revision_complex
  # mat.p.CAS[4,8] <- rowSums(mat.p.CAS[4, -4])
  mat.p.CAS[4,9] <- p.death + p.tkr_normal.death
  
  mat.p.CAS[5,6] <- p.tkr_normal_after_revision.revision_simple
  mat.p.CAS[5,7] <- p.tkr_normal_after_revision.revision_complex
  mat.p.CAS[5,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[5,5] <- 1-sum(mat.p.CAS[5, -c(5)])
  
  mat.p.CAS[6,3] <- p.revision_simple.tkr_minor
  mat.p.CAS[6,4] <- p.revision_simple.tkr_serious.CAS
  # mat.p.CAS[6,6] <- rowSums(mat.p.CAS[6, -6])
  mat.p.CAS[6,8] <- p.revision_simple.other_trt
  mat.p.CAS[6,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[6,5] <- 1-sum(mat.p.CAS[6, -c(5)])
  
  mat.p.CAS[7,4] <- p.revision_complex.tkr_serious.CAS
  #mat.p.CAS[7,7] <- rowSums(mat.p.CAS[7, -7])
  mat.p.CAS[7,9] <- p.death + p.tkr_normal_after_revision.death
  mat.p.CAS[7,5] <- 1-sum(mat.p.CAS[7,-c(5)])
  
  mat.p.CAS[8,3] <- p.other_trt.tkr_minor
  mat.p.CAS[8,4] <- p.other_trt.tkr_serious.CAS
  mat.p.CAS[8,7] <- p.other_trt.revision_simple
  #mat.p.CAS[8,8] <- p.other_trt.revision_simple
  mat.p.CAS[8,9] <- p.death + p.tkr_normal.death
  mat.p.CAS[8,5] <- 1-sum(mat.p.CAS[8, -c(5)])
  
  #check matrix - rowsum == 1
  # mat <- rowSums(mat.p.CAS)
  # print(mat)
  # 
  # 
  #utilities
  utility <- c(u.tkr,
               u.tkr_normal,
               u.tkr_minor,
               u.tkr_serious,
               u.tkr_normal_after_revision,
               u.revision_simple,
               u.revision_complex,
               u.other_trt,
               0)
  utility <- as.matrix(utility)
  utility <-t(matrix(rep(utility,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  cost <- c(c.tkr,
            c.tkr_normal,
            c.tkr_minor,
            c.tkr_serious,
            c.tkr_normal_after_revision,
            c.revision_simple,
            c.revision_complex,
            c.other_trt,
            c.death)
  cost <- as.matrix(cost)
  cost <-t(matrix(rep(cost,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  
  cost.CAS <- c(c.tkr+c.CAS,
                c.tkr_normal,
                c.tkr_minor,
                c.tkr_serious,
                c.tkr_normal_after_revision,
                c.revision_simple+c.CAS,
                c.revision_complex+c.CAS,
                c.other_trt,
                c.death)
  cost.CAS <- as.matrix(cost.CAS)
  cost.CAS <-t(matrix(rep(cost.CAS,time.horizon+1),nrow=num.states,ncol=time.horizon+1))
  
  #discount weight
  disc.weight <- lapply(c(0:time.horizon),discount,disc=disc)
  disc.weight <- as.matrix(unlist(disc.weight))
  disc.weight <-matrix(rep(disc.weight,num.states),nrow=time.horizon+1,ncol=num.states)
  
  #running
  cohort <- matrix(rep(NA,num.states*(time.horizon+1)),nrow=time.horizon+1, ncol=num.states)
  cohort.CAS <- matrix(rep(NA,num.states*(time.horizon+1)),nrow=time.horizon+1, ncol=num.states)
  
  cohort[1,] <- init 
  cohort.CAS[1,] <- init 
  
  for(i in 2:(time.horizon+1)){
    cohort[i,] <-cohort[i-1,]%*%mat.p
    cohort.CAS[i,] <- cohort.CAS[i-1,]%*%mat.p.CAS
  }
  
  #calculate QoI
  LE <- sum(cohort[,1:8])/n
  mean.QALY <- sum((cohort*utility)*disc.weight)/n
  mean.cost <- sum((cohort*cost)*disc.weight)/n
  # QALY <- cohort%*%utility*disc.weight
  # cost <- cohort%*%cost*disc.weight
  QALY <- cohort*utility*disc.weight
  cost <- cohort*cost*disc.weight
  
  LE.CAS <- sum(cohort.CAS[,1:8])/n
  mean.QALY.CAS <- sum((cohort.CAS*utility)*disc.weight)/n
  mean.cost.CAS <- sum((cohort.CAS*cost.CAS)*disc.weight)/n
  QALY.CAS <- cohort.CAS*utility*disc.weight
  cost.CAS <- cohort.CAS*cost.CAS*disc.weight
  
  NMB <- mean.QALY*WTP - mean.cost
  NMB.CAS <- mean.QALY.CAS*WTP - mean.cost.CAS
  INMB <- NMB.CAS-NMB
  return(INMB)
}

fun_pba <- function(list.params, #list of parameters that are varied - PSA
                                              #list.params.other, # list ofparams - not varied
                                              disc ,
                                              num.states,
                                              time.horizon,
                                              init,
                                              n,
                                              WTP,
                                              index_slice, pbox.int, index.array){
  # trt.index, #which strategy: 1: novel TAH, 2:syncardia
  # #outcome.index=1, #which outcome: 1:cost, 2: effect
  
  #########################################
  #rr.CAS,
  c.CAS <-  unlist(list.params[1])
  p.tkr_init.tkr_normal <- unlist( list.params[2])
  p.tkr_init.tkr_minor <- unlist( list.params[3])
  p.tkr_init.tkr_serious<- unlist( list.params[4])
  p.tkr_normal.tkr_minor <- unlist( list.params[5])
  p.tkr_normal.tkr_serious <- unlist( list.params[6])
  p.tkr_normal.tkr_normal <-  unlist(list.params[7])
  p.tkr_minor.tkr_normal <- unlist( list.params[8])
  #p.tkr_minor.tkr_serious, 
  p.tkr_minor.revision_simple <- unlist( list.params[9])
  p.tkr_minor.other_trt <- unlist( list.params[10])
  p.tkr_minor.tkr_minor <-  unlist(list.params[11])
  #p.tkr_serious.tkr_minor,
  p.tkr_serious.revision_simple <-  unlist(list.params[12])
  p.tkr_serious.revision_complex <- unlist( list.params[13])
  p.tkr_serious.other_trt <-  unlist(list.params[14])
  p.tkr_normal_after_revision.revision_simple <- unlist( list.params[15])
  p.tkr_normal_after_revision.revision_complex <-  unlist(list.params[16])
  p.tkr_normal_after_revision.tkr_normal_after_revision <-  unlist( list.params[17])
  p.revision_simple.tkr_minor <- unlist( list.params[18])
  p.revision_simple.tkr_serious <- unlist( list.params[19])
  p.revision_simple.tkr_normal_after_revision <- unlist(  list.params[20])
  #p.revision_simple.other_trt,
  p.revision_complex.tkr_serious <-  unlist(list.params[21])
  p.revision_complex.tkr_normal_after_revision <-  unlist(list.params[22])
  #p.other_trt.tkr_minor,p.other_trt.tkr_serious <- list.params[1]
  p.other_trt.tkr_normal_after_revision <- unlist(  list.params[23])
  #p.other_trt.revision_simple,
  p.tkr_normal.death <- unlist( list.params[24])
  p.tkr_normal_after_revision.death <-  unlist(list.params[25])
  p.death <-  unlist(list.params[26])
  c.tkr <-  unlist(list.params[27])
  c.tkr_normal <-  unlist(list.params[28])
  c.tkr_minor <-  unlist(list.params[29])
  c.tkr_serious <- unlist( list.params[30])
  c.tkr_normal_after_revision <-  unlist( list.params[31])
  c.revision_simple <- unlist(  list.params[32])
  c.revision_complex <-  unlist( list.params[33])
  c.other_trt <- unlist(  list.params[34])
  c.death <-  unlist( list.params[35])
  u.tkr<- unlist(  list.params[36])
  u.tkr_normal <- unlist(  list.params[37])
  u.tkr_minor <-  unlist( list.params[38])
  u.tkr_serious <-  unlist( list.params[39])
  u.tkr_normal_after_revision <- unlist(  list.params[40])
  u.revision_simple <-  unlist( list.params[41])
  u.revision_complex <-  unlist( list.params[42])
  u.other_trt <-  unlist( list.params[43])
  
 
  #########################################
  dim_slice <- length(index_slice) #size of set K (equation 12) or the number of possible partitions
  result_list <- c(rep(0,(dim_slice*3)))
  dim(result_list) <- c(dim_slice,3)
  #pbox.int <- pbox.int[ , , ,trt.index]
  #pbox.int <- pbox.int[ , , ,trt.index]
  num.pbox.param <- dim(pbox.int)[1] #get the number of pbox params
  for (i in 1:dim_slice){
    k <- index_slice[i] #get current index from the set of indices (K in equation 12)
    #getting lower and upper bounds of each parameter
    # lbound <-rep(0,num.pbox.param)
    # ubound <-rep(0,num.pbox.param)
    lbound <- matrix(rep(0,num.pbox.param*2),nrow=num.pbox.param,ncol=2)
    ubound <- matrix(rep(0,num.pbox.param*2),nrow=num.pbox.param,ncol=2)
    
    for (j in 1:num.pbox.param){
      # lbound <- c(pbox.int[3,index.array[k,1],1],pbox.int[5,index.array[k,2],1])
      # ubound <- c(pbox.int[3,index.array[k,1],2],pbox.int[5,index.array[k,2],2])
      for (jj in 1:2){
        lbound[j,jj] <- pbox.int[j,index.array[k,j],1,jj]
        ubound[j,jj] <- pbox.int[j,index.array[k,j],2,jj]
      }
    }
    
    #using auglag
    objective.index=2
    x.0 <-(lbound[,1]+ubound[,1])/2
    # print(lbound[,1])
    # print(x.0)
    # print(ubound[,1])
    # print("---------------------")
    #temporary fix due to lbound is higher than ubound (vice versa) - not yet fixed as of 15 May 2021
    bound.compare <- (lbound[,1] < ubound[,1])
    lbound_2 <- lbound
    ubound_2 <- ubound
    
    for (i.bound in 1:length(ubound[,1])){
      if(bound.compare[i.bound]==FALSE){
        lbound_2[i.bound,1] <- ubound[i.bound,1]
        ubound_2[i.bound,1] <- lbound[i.bound,1]
      }
    }
    ubound <- ubound_2
    lbound <- lbound_2
    
    optim.max <- bobyqa(lbound[,1],TKR_model,lbound[,1],ubound[,1],
                        nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000),
                        objective.index,
                        #rr.CAS,
                        c.CAS,
                        p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                        p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                        p.tkr_minor.tkr_normal,
                        #p.tkr_minor.tkr_serious, 
                        p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                        #p.tkr_serious.tkr_minor,
                        p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                        p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                        p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,
                        #p.revision_simple.other_trt,
                        p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                        #p.other_trt.tkr_minor,p.other_trt.tkr_serious,
                        p.other_trt.tkr_normal_after_revision,
                        #p.other_trt.revision_simple,
                        p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                        c.tkr,
                        c.tkr_normal,
                        c.tkr_minor,
                        c.tkr_serious,
                        c.tkr_normal_after_revision,
                        c.revision_simple,
                        c.revision_complex,
                        c.other_trt,
                        c.death,
                        u.tkr,
                        u.tkr_normal,
                        u.tkr_minor,
                        u.tkr_serious,
                        u.tkr_normal_after_revision,
                        u.revision_simple,
                        u.revision_complex,
                        u.other_trt,
                        disc ,
                        num.states,
                        time.horizon,
                        init,
                        n,
                        WTP)
    
    #min
    objective.index=1
    optim.min <- bobyqa(lbound[,1],TKR_model,lbound[,1],ubound[,1],
                        nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000),
                        objective.index,
                        #rr.CAS,
                        c.CAS,
                        p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                        p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                        p.tkr_minor.tkr_normal,
                        #p.tkr_minor.tkr_serious, 
                        p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                        #p.tkr_serious.tkr_minor,
                        p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                        p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                        p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,
                        #p.revision_simple.other_trt,
                        p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                        #p.other_trt.tkr_minor,p.other_trt.tkr_serious,
                        p.other_trt.tkr_normal_after_revision,
                        #p.other_trt.revision_simple,
                        p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                        c.tkr,
                        c.tkr_normal,
                        c.tkr_minor,
                        c.tkr_serious,
                        c.tkr_normal_after_revision,
                        c.revision_simple,
                        c.revision_complex,
                        c.other_trt,
                        c.death,
                        u.tkr,
                        u.tkr_normal,
                        u.tkr_minor,
                        u.tkr_serious,
                        u.tkr_normal_after_revision,
                        u.revision_simple,
                        u.revision_complex,
                        u.other_trt,
                        disc ,
                        num.states,
                        time.horizon,
                        init,
                        n,
                        WTP)
    # 
    # optim.min <- directL(TAH.model.single.run_par_inmb,lbound[,1],ubound[,1],
    #                     nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000), 
    #                     objective.index=objective.index, #change sign objective function: 1: positive, 2: negative
    #                     trt.index=1, #which strategy: 1: novel TAH, 2:syncardia
    #                     tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
    #                     #short.term.mortality.crmt=short.term.mortality.crmt, 
    #                     long.term.mortality.crmt=long.term.mortality.crmt,
    #                     short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
    #                     tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
    #                     p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
    #                     # pre-transplantation
    #                     LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
    #                     p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
    #                     p.dis=p.dis,
    #                     c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
    #                     c.m1=c.m1,
    #                     # post-transplantation
    #                     c.HTx=c.HTx,c_pc.post=c_pc.post,
    #                     e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,    
    #                     e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm, 
    #                     # post-transplantation
    #                     e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
    #                     TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc, de=de,
    #                     WTP=WTP,
    #                     randomized = TRUE)
    # optim.min_2 <- bobyqa(lbound[,1],TAH.model.single.run_par,lbound[,1],ubound[,1],
    #                       nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000), 
    #                       objective.index=objective.index, #change sign objective function: 1: positive, 2: negative
    #                       trt.index=2, #which strategy: 1: novel TAH, 2:syncardia
    #                       tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
    #                       #short.term.mortality.crmt=short.term.mortality.crmt, 
    #                       long.term.mortality.crmt=long.term.mortality.crmt,
    #                       short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
    #                       tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
    #                       p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
    #                       # pre-transplantation
    #                       LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
    #                       p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
    #                       p.dis=p.dis,
    #                       c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
    #                       c.m1=c.m1,
    #                       # post-transplantation
    #                       c.HTx=c.HTx,c_pc.post=c_pc.post,
    #                       e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,    
    #                       e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm, 
    #                       # post-transplantation
    #                       e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
    #                       TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc, de=de,
    #                       WTP)
    
    #optim.results[k,2] <- optim.min$value
    #calculate weights (assume independence)
    w=1
    for (j in 1:num.pbox.param){
      #w1 <- pbox.int[3,index.array[k,1],3]
      #w2 <- pbox.int[5,index.array[k,2],3]
      w <- w*pbox.int[j,index.array[k,j],3,1]
    }
    #optim.results[k,3] <- w1*w2
    result_list[i,] <- c(optim.min$value,optim.max$value, w)
  }
  return(result_list)
}

#parallel over set K (set of multi-indices) 5 august 2021
fun_pba_par <- function(list.index, list.params, #list of parameters that are varied - PSA
                    #list.params.other, # list ofparams - not varied
                    disc ,
                    num.states,
                    time.horizon,
                    init,
                    n,
                    WTP,
                    index_slice, pbox.int, index.array){
  # trt.index, #which strategy: 1: novel TAH, 2:syncardia
  # #outcome.index=1, #which outcome: 1:cost, 2: effect
  
  #########################################
  #rr.CAS,
  c.CAS <-  unlist(list.params[1])
  p.tkr_init.tkr_normal <- unlist( list.params[2])
  p.tkr_init.tkr_minor <- unlist( list.params[3])
  p.tkr_init.tkr_serious<- unlist( list.params[4])
  p.tkr_normal.tkr_minor <- unlist( list.params[5])
  p.tkr_normal.tkr_serious <- unlist( list.params[6])
  p.tkr_normal.tkr_normal <-  unlist(list.params[7])
  p.tkr_minor.tkr_normal <- unlist( list.params[8])
  #p.tkr_minor.tkr_serious, 
  p.tkr_minor.revision_simple <- unlist( list.params[9])
  p.tkr_minor.other_trt <- unlist( list.params[10])
  p.tkr_minor.tkr_minor <-  unlist(list.params[11])
  #p.tkr_serious.tkr_minor,
  p.tkr_serious.revision_simple <-  unlist(list.params[12])
  p.tkr_serious.revision_complex <- unlist( list.params[13])
  p.tkr_serious.other_trt <-  unlist(list.params[14])
  p.tkr_normal_after_revision.revision_simple <- unlist( list.params[15])
  p.tkr_normal_after_revision.revision_complex <-  unlist(list.params[16])
  p.tkr_normal_after_revision.tkr_normal_after_revision <-  unlist( list.params[17])
  p.revision_simple.tkr_minor <- unlist( list.params[18])
  p.revision_simple.tkr_serious <- unlist( list.params[19])
  p.revision_simple.tkr_normal_after_revision <- unlist(  list.params[20])
  #p.revision_simple.other_trt,
  p.revision_complex.tkr_serious <-  unlist(list.params[21])
  p.revision_complex.tkr_normal_after_revision <-  unlist(list.params[22])
  #p.other_trt.tkr_minor,p.other_trt.tkr_serious <- list.params[1]
  p.other_trt.tkr_normal_after_revision <- unlist(  list.params[23])
  #p.other_trt.revision_simple,
  p.tkr_normal.death <- unlist( list.params[24])
  p.tkr_normal_after_revision.death <-  unlist(list.params[25])
  p.death <-  unlist(list.params[26])
  c.tkr <-  unlist(list.params[27])
  c.tkr_normal <-  unlist(list.params[28])
  c.tkr_minor <-  unlist(list.params[29])
  c.tkr_serious <- unlist( list.params[30])
  c.tkr_normal_after_revision <-  unlist( list.params[31])
  c.revision_simple <- unlist(  list.params[32])
  c.revision_complex <-  unlist( list.params[33])
  c.other_trt <- unlist(  list.params[34])
  c.death <-  unlist( list.params[35])
  u.tkr<- unlist(  list.params[36])
  u.tkr_normal <- unlist(  list.params[37])
  u.tkr_minor <-  unlist( list.params[38])
  u.tkr_serious <-  unlist( list.params[39])
  u.tkr_normal_after_revision <- unlist(  list.params[40])
  u.revision_simple <-  unlist( list.params[41])
  u.revision_complex <-  unlist( list.params[42])
  u.other_trt <-  unlist( list.params[43])
  
  
  #########################################
  #dim_slice <- length(index_slice) #size of set K (equation 12) or the number of possible partitions
  dim_slice <- 1
  result_list <- c(rep(0,(dim_slice*3)))
  dim(result_list) <- c(dim_slice,3)
  # 
  #result_list <- list()
  #dim(result_list) <- c(dim_slice,3)
  #pbox.int <- pbox.int[ , , ,trt.index]
  #pbox.int <- pbox.int[ , , ,trt.index]
  num.pbox.param <- dim(pbox.int)[1] #get the number of pbox params
  for (i in 1:dim_slice){
    k <- list.index[[8]]
    print("---------------------------")
    print(k)
    print("---------------------------")
    #k <- index_slice[i] #get current index from the set of indices (K in equation 12)
    #getting lower and upper bounds of each parameter
    # lbound <-rep(0,num.pbox.param)
    # ubound <-rep(0,num.pbox.param)
    lbound <- matrix(rep(0,num.pbox.param*2),nrow=num.pbox.param,ncol=2)
    ubound <- matrix(rep(0,num.pbox.param*2),nrow=num.pbox.param,ncol=2)
    
    for (j in 1:num.pbox.param){
      # lbound <- c(pbox.int[3,index.array[k,1],1],pbox.int[5,index.array[k,2],1])
      # ubound <- c(pbox.int[3,index.array[k,1],2],pbox.int[5,index.array[k,2],2])
      for (jj in 1:2){
        print(pbox.int[j,index.array[k,j],1,jj])
        lbound[j,jj] <- pbox.int[j,index.array[k,j],1,jj]
        ubound[j,jj] <- pbox.int[j,index.array[k,j],2,jj]
      }
    }
    
    #using auglag
    objective.index=2
    x.0 <-(lbound[,1]+ubound[,1])/2
    # print(lbound[,1])
    # print(x.0)
    # print(ubound[,1])
    # print("---------------------")
    #temporary fix due to lbound is higher than ubound (vice versa) - not yet fixed as of 15 May 2021
    bound.compare <- (lbound[,1] < ubound[,1])
    lbound_2 <- lbound
    ubound_2 <- ubound
    
    for (i.bound in 1:length(ubound[,1])){
      if(bound.compare[i.bound]==FALSE){
        lbound_2[i.bound,1] <- ubound[i.bound,1]
        ubound_2[i.bound,1] <- lbound[i.bound,1]
      }
    }
    ubound <- ubound_2
    lbound <- lbound_2
    
    optim.max <- bobyqa(lbound[,1],TKR_model,lbound[,1],ubound[,1],
                        nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000),
                        objective.index,
                        #rr.CAS,
                        c.CAS,
                        p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                        p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                        p.tkr_minor.tkr_normal,
                        #p.tkr_minor.tkr_serious, 
                        p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                        #p.tkr_serious.tkr_minor,
                        p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                        p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                        p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,
                        #p.revision_simple.other_trt,
                        p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                        #p.other_trt.tkr_minor,p.other_trt.tkr_serious,
                        p.other_trt.tkr_normal_after_revision,
                        #p.other_trt.revision_simple,
                        p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                        c.tkr,
                        c.tkr_normal,
                        c.tkr_minor,
                        c.tkr_serious,
                        c.tkr_normal_after_revision,
                        c.revision_simple,
                        c.revision_complex,
                        c.other_trt,
                        c.death,
                        u.tkr,
                        u.tkr_normal,
                        u.tkr_minor,
                        u.tkr_serious,
                        u.tkr_normal_after_revision,
                        u.revision_simple,
                        u.revision_complex,
                        u.other_trt,
                        disc ,
                        num.states,
                        time.horizon,
                        init,
                        n,
                        WTP)
    
    #min
    objective.index=1
    optim.min <- bobyqa(lbound[,1],TKR_model,lbound[,1],ubound[,1],
                        nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000),
                        objective.index,
                        #rr.CAS,
                        c.CAS,
                        p.tkr_init.tkr_normal, p.tkr_init.tkr_minor,p.tkr_init.tkr_serious,
                        p.tkr_normal.tkr_minor, p.tkr_normal.tkr_serious,p.tkr_normal.tkr_normal,
                        p.tkr_minor.tkr_normal,
                        #p.tkr_minor.tkr_serious, 
                        p.tkr_minor.revision_simple,p.tkr_minor.other_trt,p.tkr_minor.tkr_minor,
                        #p.tkr_serious.tkr_minor,
                        p.tkr_serious.revision_simple,p.tkr_serious.revision_complex,p.tkr_serious.other_trt,
                        p.tkr_normal_after_revision.revision_simple,p.tkr_normal_after_revision.revision_complex,p.tkr_normal_after_revision.tkr_normal_after_revision,
                        p.revision_simple.tkr_minor,p.revision_simple.tkr_serious,p.revision_simple.tkr_normal_after_revision,
                        #p.revision_simple.other_trt,
                        p.revision_complex.tkr_serious,p.revision_complex.tkr_normal_after_revision,
                        #p.other_trt.tkr_minor,p.other_trt.tkr_serious,
                        p.other_trt.tkr_normal_after_revision,
                        #p.other_trt.revision_simple,
                        p.tkr_normal.death,p.tkr_normal_after_revision.death,p.death,
                        c.tkr,
                        c.tkr_normal,
                        c.tkr_minor,
                        c.tkr_serious,
                        c.tkr_normal_after_revision,
                        c.revision_simple,
                        c.revision_complex,
                        c.other_trt,
                        c.death,
                        u.tkr,
                        u.tkr_normal,
                        u.tkr_minor,
                        u.tkr_serious,
                        u.tkr_normal_after_revision,
                        u.revision_simple,
                        u.revision_complex,
                        u.other_trt,
                        disc ,
                        num.states,
                        time.horizon,
                        init,
                        n,
                        WTP)
    # 
    # optim.min <- directL(TAH.model.single.run_par_inmb,lbound[,1],ubound[,1],
    #                     nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000), 
    #                     objective.index=objective.index, #change sign objective function: 1: positive, 2: negative
    #                     trt.index=1, #which strategy: 1: novel TAH, 2:syncardia
    #                     tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
    #                     #short.term.mortality.crmt=short.term.mortality.crmt, 
    #                     long.term.mortality.crmt=long.term.mortality.crmt,
    #                     short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
    #                     tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
    #                     p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
    #                     # pre-transplantation
    #                     LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
    #                     p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
    #                     p.dis=p.dis,
    #                     c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
    #                     c.m1=c.m1,
    #                     # post-transplantation
    #                     c.HTx=c.HTx,c_pc.post=c_pc.post,
    #                     e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,    
    #                     e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm, 
    #                     # post-transplantation
    #                     e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
    #                     TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc, de=de,
    #                     WTP=WTP,
    #                     randomized = TRUE)
    # optim.min_2 <- bobyqa(lbound[,1],TAH.model.single.run_par,lbound[,1],ubound[,1],
    #                       nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000), 
    #                       objective.index=objective.index, #change sign objective function: 1: positive, 2: negative
    #                       trt.index=2, #which strategy: 1: novel TAH, 2:syncardia
    #                       tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
    #                       #short.term.mortality.crmt=short.term.mortality.crmt, 
    #                       long.term.mortality.crmt=long.term.mortality.crmt,
    #                       short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
    #                       tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
    #                       p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
    #                       # pre-transplantation
    #                       LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
    #                       p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
    #                       p.dis=p.dis,
    #                       c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
    #                       c.m1=c.m1,
    #                       # post-transplantation
    #                       c.HTx=c.HTx,c_pc.post=c_pc.post,
    #                       e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,    
    #                       e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm, 
    #                       # post-transplantation
    #                       e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
    #                       TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc, de=de,
    #                       WTP)
    
    #optim.results[k,2] <- optim.min$value
    #calculate weights (assume independence)
    w=1
    for (j in 1:num.pbox.param){
      #w1 <- pbox.int[3,index.array[k,1],3]
      #w2 <- pbox.int[5,index.array[k,2],3]
      w <- w*pbox.int[j,index.array[k,j],3,1]
    }
    #optim.results[k,3] <- w1*w2
    result_list <- c(optim.min$value,optim.max$value, w)
  }
  return(result_list)
}


estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha,beta))
}