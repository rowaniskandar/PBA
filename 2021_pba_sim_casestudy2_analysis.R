#R file for case study 2
#Probability bound analysis: A novel approach for quantifying parameter uncertainty in decision-analytic modeling and cost-effectiveness analysis
#author: Rowan Iskandar
#analysis file
#based on
#Dong, Hengjin, and Martin Buxton. "Early assessment of the likely cost-effectiveness of a new technology: a Markov model with probabilistic sensitivity analysis of computer-assisted total knee replacement." International Journal of Technology Assessment in Health Care 22.2 (2006): 191-202.
#last updated: 13 August 2021
currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  
setwd(currpath)

library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
#library(INLA, lib.loc = "C:/R_libraries")
#library(flexsurv)

library(plyr)
library(ggplot2)
library(nloptr) #optimization
library(doParallel)  
library(future)
library(future.apply)
library(purrr)
library(truncnorm)
# library(gmailr)
library(support)
library(openxlsx)
library(beepr)
library(RColorBrewer)

source("2021_pba_sim_casestudy2_model.R") #make sure you put this model file in the same folder level
###################################################################################
ptm1 <- proc.time()
pbox.parallel <- 1 #1 = parallel mode over the sub-intervals

#estimates from the Dong 2006
p.tkr_init.tkr_normal <- 0.94220
p.tkr_init.tkr_minor <- 0.04285
p.tkr_init.tkr_serious <- 0.01495
#no death

p.tkr_normal.tkr_minor <- .01385
p.tkr_normal.tkr_serious <- .00921
p.tkr_normal.death <- 0.00046 + 0.00341

p.tkr_normal.tkr_normal <- .97307  # seems 1-#

p.tkr_minor.tkr_normal <- .94236
p.tkr_minor.tkr_serious <- .00921
p.tkr_minor.revision_simple <- .00250
p.tkr_minor.other_trt <- .01701
p.tkr_minor.death <- 0.00046 + 0.00341

p.tkr_minor.tkr_minor <- .02505  #seems 1-#

p.tkr_serious.tkr_minor <- .01385
p.tkr_serious.revision_simple <- .00523
p.tkr_serious.revision_complex <- .02469
p.tkr_serious.other_trt <- .95236 # seems 1-#
p.tkr_serious.death <- 0.00341 + 0.00046
#cannot stay

p.tkr_normal_after_revision.revision_simple <- .01038
p.tkr_normal_after_revision.revision_complex <- .02003
p.tkr_normal_after_revision.death <- .00151+ 0.00341

p.tkr_normal_after_revision.tkr_normal_after_revision <- .96468 #seems  1-#

p.revision_simple.tkr_minor <- .00816
p.revision_simple.tkr_serious <- .01590
p.revision_simple.tkr_normal_after_revision <- .95400 #seems  1-#
p.revision_simple.other_trt <- .01701
p.revision_simple.death <- .00151+.00341
#cannot stay

p.revision_complex.tkr_serious <- .02545
p.revision_complex.tkr_normal_after_revision <- .96963 #seems  1-#
p.revision_complex.death <- 0.00151+.00341
#cannot stay

p.other_trt.tkr_minor <- .01385
p.other_trt.tkr_serious <- .00921
p.other_trt.tkr_normal_after_revision <- .97057 # seems 1-#
p.other_trt.revision_simple <- .00250
p.other_trt.death <- 0.00341 + 0.00046

p.death <- 0.00341
p.tkr_normal.death <- 0.00046
p.tkr_normal_after_revision.death <- .00151

rr.CAS <- 0.66
c.CAS <- 235
  
c.tkr <- 5197
c.tkr_normal <- 0
c.tkr_minor <- 0
c.tkr_serious <-0
c.tkr_normal_after_revision <-0
c.revision_simple <- 6234
c.revision_complex <- 7326
c.other_trt <- 2844
c.death <- 0
 
u.tkr <- 0.72
u.tkr_normal <- 0.78
u.tkr_minor <- 0.66
u.tkr_serious <- 0.35
u.tkr_normal_after_revision <- 0.68
u.revision_simple <- 0.66
u.revision_complex <- 0.51
u.other_trt <- 0.72

disc  <- 0.035
num.states  <- 9
time.horizon <- 10*12
WTP <- 30000
n <- 1000
init <- c(n,0,0,0,0,0,0,0,0)

min.x <- -1200000 #min(df.nmb.upper,df.nmb.lower)
max.x <- 700000#max(df.nmb.upper,df.nmb.lower)
n.x <-200000
x.values <- seq(min.x,max.x,length.out=n.x)
################################################################################
#start of PBA
##############################################################################
set.seed(1234)

pbox.num <-5 #number of discretization
n.sp <- 25
n.MC <- 1000
n.sims <- n.sp
sp.yes <- 1 #1: do monte carlo simulation 0: no monte carlo

sp = 0 #sp =0 uses monte carlo approach sp=1 uses support point
if(sp==0){n.sims <- n.MC}
MC.result <- vector(mode = "list", length = n.sims)

#1. uncertainty representation
num.pbox.params <- 7
# p.tkr_minor.tkr_serious
# p.tkr_serious.tkr_minor
# p.revision_simple.other_trt
# p.other_trt.tkr_minor
# p.other_trt.tkr_serious
# p.other_trt.revision_simple

#PBA parameters to vary (#Carlo# change these values for the summary statistics on PBA parameter)
#default: summary statistics calculated from PSA data
#pba params:
rr.CAS.mean <- rr.CAS
rr.CAS.min <- 0
rr.CAS.max <- 1

p.tkr_minor.tkr_serious.mean <- p.tkr_minor.tkr_serious
p.tkr_minor.tkr_serious.min <- 0.00327
p.tkr_minor.tkr_serious.max <- 0.02704

p.tkr_serious.tkr_minor.mean <- p.tkr_serious.tkr_minor
p.tkr_serious.tkr_minor.min <- .00428
p.tkr_serious.tkr_minor.max <- 03200

p.revision_simple.other_trt.mean <- p.revision_simple.other_trt
p.revision_simple.other_trt.min <- .00621
p.revision_simple.other_trt.max <- .03753

p.other_trt.tkr_minor.mean <- p.other_trt.tkr_minor
p.other_trt.tkr_minor.min <- .00633
p.other_trt.tkr_minor.max <- .03213

p.other_trt.tkr_serious.mean <- p.other_trt.tkr_serious
p.other_trt.tkr_serious.min <- .00232
p.other_trt.tkr_serious.max <- .02643

p.other_trt.revision_simple.mean <- p.other_trt.revision_simple
p.other_trt.revision_simple.min <- .00007
p.other_trt.revision_simple.max <- .01361

#bound scenarios (variable: bound):
#"extreme": [0,1]
#"PSA": [min, max]
#"user" supplied

#summary statistics supplied, choose one of the following (variable: statistics)
#minmax: min, max
#"mean": min, max, mean
#"std": min, max, mean, std
#"median", min, max, median
#"median_std": min, max, mean, std, median
#num.pbox.int <-12

#pbox.int.loop <- c(1,2,3,4,5,6,7,8)
#pbox.int.loop <- c(2,3,4,5,12)
pbox.int.loop <-c(pbox.num)

# statistics.loop <-c("minmax","median","mean","std","median_std")
statistics.loop <-c("mean") #we have data only on min, max, mean

bound <- "PBA" #use published ranges as in Song 2006

#statistics.loop <-c("minmax","mean","std")

#preparing matrices for results by combinations of pbox intervals and summary statistics
#columns: number of x, bounds: lower upper, statistics, pbox int
cdf.all <- array(c(rep(0,n.x*2*length(statistics.loop),length(pbox.int.loop))), dim=c(n.x,2,length(statistics.loop),length(pbox.int.loop)))
################################################################################
#start of loop over pbox intervals and over which statistics
for (i.int in 1:length(pbox.int.loop)){
  for (i.stat in 1:length(statistics.loop)){
    
    statistics <- statistics.loop[i.stat]
    num.pbox.int <- pbox.int.loop[i.int] #number of p-box discretization for UP
    #combination.med.std <- "median_std"
    #################################################
    if(bound=="extreme"){
      
      rr.CAS.max <- 1
      rr.CAS.min <- 0
      
      p.tkr_minor.tkr_serious.min <- 0.0001
      p.tkr_minor.tkr_serious.max <- 0.2
      
      p.tkr_serious.tkr_minor.min <- 0.0001
      p.tkr_serious.tkr_minor.max <- 0.2
      
      p.revision_simple.other_trt.min <- 0.0001
      p.revision_simple.other_trt.max <- 0.2
      
      p.other_trt.tkr_minor.min <- 0.0001
      p.other_trt.tkr_minor.max <- 0.2
      
      p.other_trt.tkr_serious.min <- 0.0001
      p.other_trt.tkr_serious.max <- 0.2
      
      p.other_trt.revision_simple.min <- 0.00001
      p.other_trt.revision_simple.max <- 0.2
    } else if(bound=="PSA"){
      
      rr.CAS.max <- 1
      rr.CAS.min <- 0.3
      
      p.tkr_minor.tkr_serious.min <- 0.00327
      p.tkr_minor.tkr_serious.max <- 0.02704
      
      p.tkr_serious.tkr_minor.min <- 0.00428
      p.tkr_serious.tkr_minor.max <- 0.03200
      
      p.revision_simple.other_trt.min <- 0.00621
      p.revision_simple.other_trt.max <- 0.03753
      
      p.other_trt.tkr_minor.min <- 0.00633
      p.other_trt.tkr_minor.max <- 0.03213
      
      p.other_trt.tkr_serious.min <- 0.00232
      p.other_trt.tkr_serious.max <- 0.02643
      
      p.other_trt.revision_simple.min <- 0.00007
      p.other_trt.revision_simple.max <- 0.01361
    } 

    
    ##################################################################
    #table of parameters which are represented by pbox
    #rows: parameter
    #columns: min, max, median, mean, std
    #3rd dim: treatment
    #params.table <- matrix(rep(0,num.pbox.params*5),nrow=num.pbox.params,ncol=5)
    params.table <- array(rep(0,num.pbox.params*5*2),dim=c(num.pbox.params, 5, 2))
    
    mean.std <- 0.2 #std as a fraction of mean
    #crmt params for the PBA params
    # params.table[1, ,1] = c(0,1,9999,short.term.mortality.crmt.1.mean, short.term.mortality.crmt.1.std )
    # params.table[2, ,1] = c(0,1,9999,short.term.mortality.crmt.2.mean, short.term.mortality.crmt.2.std )
    # params.table[3, ,1]= c(0,1,9999,short.term.mortality.crmt.3.mean, short.term.mortality.crmt.3.std )
    # params.table[4, ,1]= c(0,1,9999,p.majdm.crmt.mean,p.majdm.crmt.std)
    
    params.table[1, ,1] = c(rr.CAS.min,rr.CAS.max,0,
                            rr.CAS.mean, rr.CAS.mean*mean.std )
    params.table[2, ,1] = c(p.tkr_minor.tkr_serious.min,p.tkr_minor.tkr_serious.max,0,
                            p.tkr_minor.tkr_serious.mean, p.tkr_minor.tkr_serious.mean*mean.std )
    params.table[3, ,1] = c(p.tkr_serious.tkr_minor.min,p.tkr_serious.tkr_minor.max,0,
                            p.tkr_serious.tkr_minor.mean, p.tkr_serious.tkr_minor.mean*mean.std )
    params.table[4, ,1]= c(p.revision_simple.other_trt.min,p.revision_simple.other_trt.max,0,
                           p.revision_simple.other_trt.mean, p.revision_simple.other_trt.mean*mean.std )
    params.table[5, ,1]= c(p.other_trt.tkr_minor.min,p.other_trt.tkr_minor.max,0,
                           p.other_trt.tkr_minor.mean,p.other_trt.tkr_minor.mean*mean.std)
    params.table[6, ,1]= c(p.other_trt.tkr_serious.min,p.other_trt.tkr_serious.max,0,
                           p.other_trt.tkr_serious.mean,p.other_trt.tkr_serious.mean*mean.std)
    params.table[7, ,1]= c(p.other_trt.revision_simple.min,p.other_trt.revision_simple.max,0,
                           p.other_trt.revision_simple.mean,p.other_trt.revision_simple.mean*mean.std)

    #discretization over the cdf values (u) (NOTE: this is range of p-box)
    #partitioning u domain
    #columns: (1) lower bound - interval (2) upper bound - interval (3) width of interval
    inv.int <- matrix(rep(0,(num.pbox.int-1)*3),nrow=num.pbox.int-1,ncol=3) 
    dummy <- seq(0,1,length.out=num.pbox.int)
    for (i in 1:(num.pbox.int-1)){
      inv.int[i,1]=dummy[i]
      inv.int[i,2]=dummy[i+1]
      inv.int[i,3]=dummy[i+1]-dummy[i]
    }
    #over x values, using inverse pbox functions
    #dimensions: (1) variable, (2) intervals, (3) lower or upper bounds or weight
    #pbox.int <-array(rep(0,num.pbox.params*(num.pbox.int-1)*3), dim=c(num.pbox.params, num.pbox.int-1, 3))
    pbox.int <-array(rep(0,num.pbox.params*(num.pbox.int-1)*3*2), dim=c(num.pbox.params, num.pbox.int-1, 3,2))
    for (j in 1:num.pbox.params){
      for (i in 1:(num.pbox.int-1)){
        for (jj in 1:1){
          if(statistics=="minmax"){
            pbox.int[j,i,1,jj]=pbox_minmax_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj])
            pbox.int[j,i,2,jj]=pbox_minmax_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj])
          } else if (statistics=="median"){
            pbox.int[j,i,1,jj]=pbox_median_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj])
            pbox.int[j,i,2,jj]=pbox_median_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj])
          } else if (statistics=="mean"){
            pbox.int[j,i,1,jj]=pbox_mean_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj])
            pbox.int[j,i,2,jj]=pbox_mean_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj])
          } else if (statistics=="std"){
            pbox.int[j,i,1,jj]=pbox_std_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj],params.table[j,5,jj])
            pbox.int[j,i,2,jj]=pbox_std_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj],params.table[j,5,jj])
          } else if (statistics=="median_std"){
            pbox.int[j,i,1,jj]=pbox_median_std_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj],params.table[j,4,jj],params.table[j,5,jj])
            pbox.int[j,i,2,jj]=pbox_median_std_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj],params.table[j,4,jj],params.table[j,5,jj])
          }  
          pbox.int[j,i,3,jj]=inv.int[i,3]
        }
      }
    }
    
    #for PSA params
    if (sp.yes==1){
      if(i.int==1 & i.stat==1){
        if(sp==1){
          dist.param <- vector("list",7)
          dist.param[[1]] <- c(params.table[1,4,1]^2/params.table[1,5,1]^2,params.table[1,5,1]^2/(params.table[1,4,1]))
          
          dist.param[[2]] <- estBetaParams(params.table[2,4,1],params.table[2,5,1]^2)
          dist.param[[3]] <- estBetaParams(params.table[3,4,1],params.table[3,5,1]^2)
          dist.param[[4]] <- estBetaParams(params.table[4,4,1],params.table[4,5,1]^2)
          dist.param[[5]] <- estBetaParams(params.table[5,4,1],params.table[5,5,1]^2)
          dist.param[[6]] <- estBetaParams(params.table[6,4,1],params.table[6,5,1]^2)
          dist.param[[7]] <- estBetaParams(params.table[7,4,1],params.table[7,5,1]^2)
          
          D <- sp(n.sp,7,dist.str=c("gamma",rep("beta",6)),dist.param=dist.param)
          
          x1.sp <- D$sp[,1]
          x2.sp <- D$sp[,2]
          x3.sp <- D$sp[,3]
          x4.sp <- D$sp[,4]    
          x5.sp <- D$sp[,5]
          x6.sp <- D$sp[,6]
          x7.sp <- D$sp[,7]
          
          list_test.sp <- list()
          for (i in 1:n.sp){
            list_test.sp[[i]] <- list(x1.sp[i],x2.sp[i],x3.sp[i],x4.sp[i],x5.sp[i],x6.sp[i],x7.sp[i])
          }
        }
        else {
          x1.sp <- rgamma(n.MC, shape=(params.table[1,4,1])^2/params.table[1,5,1]^2, scale=params.table[1,5,1]^2/(params.table[1,4,1]))
          x2.sp <- rbeta(n.MC, estBetaParams(params.table[2,4,1],params.table[2,5,1]^2)[1], estBetaParams(params.table[2,4,1],params.table[2,5,1]^2)[2])
          x3.sp <- rbeta(n.MC, estBetaParams(params.table[3,4,1],params.table[3,5,1]^2)[1], estBetaParams(params.table[3,4,1],params.table[3,5,1]^2)[2])
          x4.sp <- rbeta(n.MC, estBetaParams(params.table[4,4,1],params.table[4,5,1]^2)[1], estBetaParams(params.table[4,4,1],params.table[4,5,1]^2)[2])
          x5.sp <- rbeta(n.MC, estBetaParams(params.table[5,4,1],params.table[5,5,1]^2)[1], estBetaParams(params.table[5,4,1],params.table[5,5,1]^2)[2])
          x6.sp <- rbeta(n.MC, estBetaParams(params.table[6,4,1],params.table[6,5,1]^2)[1], estBetaParams(params.table[6,4,1],params.table[6,5,1]^2)[2])
          x7.sp <- rbeta(n.MC, estBetaParams(params.table[7,4,1],params.table[7,5,1]^2)[1], estBetaParams(params.table[7,4,1],params.table[7,5,1]^2)[2])
          
          # x1.sp <- runif(n.MC, params.table[1,1,1], params.table[1,2,1])
          # x2.sp <- runif(n.MC,  params.table[2,1,1], params.table[2,2,1])
          # x3.sp <- runif(n.MC,  params.table[3,1,1], params.table[3,2,1])
          # x4.sp <- runif(n.MC, params.table[4,1,1], params.table[4,2,1])  
          # x5.sp <- runif(n.MC,  params.table[5,1,1], params.table[5,2,1])
          # x6.sp <- runif(n.MC,  params.table[6,1,1], params.table[6,2,1])
          # x7.sp <- runif(n.MC,  params.table[7,1,1], params.table[7,2,1])
          list_test.sp <- list()
          for (i in 1:n.MC){
            list_test.sp[[i]] <- list(x1.sp[i],x2.sp[i],x3.sp[i],x4.sp[i],x5.sp[i],x6.sp[i],x7.sp[i])
          }
        }
      }
    }
    

    
    params_MC <- list_test.sp
    ################################################################################
    #2. uncertainty propagation
    ################################################################################
    #create multi-indices (for each p-box parameter)
    num.param.pba <- num.pbox.params #number of model parameters with p-box
    n.indices <- (num.pbox.int-1)**num.param.pba #number of unique indices - for each parameter, there will be num.pbox.int index
    input_data <- 1:n.indices

    vec <- c(1:(num.pbox.int-1))
    lst <- lapply(numeric(num.param.pba), function(x) vec)
    index.array <- as.matrix(expand.grid(lst))
    #print(dim(index.array)[1])
    index.array.last <- rep(0,dim(index.array)[1])
    index.array <- cbind(index.array,index.array.last)
    #print(dim(index.array)[1])
    
    list.index <- list()
    index_in_index.array <- seq(1,n.indices)
    for (zz in 1:n.indices){
     list.index[[zz]] <- list(index.array[zz,1],index.array[zz,2],index.array[zz,3],index.array[zz,4],index.array[zz,5],index.array[zz,6],index.array[zz,7],index_in_index.array[zz]) 
    }
    
    ##############################################################################
    #par version
    list_test.pbox <- list()
    #fixed PSA parameters at their average values
list_test.pbox <- list(  c.CAS ,
                         p.tkr_init.tkr_normal ,
                         p.tkr_init.tkr_minor ,
                         p.tkr_init.tkr_serious,
                         p.tkr_normal.tkr_minor,
                         p.tkr_normal.tkr_serious ,
                         p.tkr_normal.tkr_normal ,
                         p.tkr_minor.tkr_normal ,
                         #p.tkr_minor.tkr_serious, 
                         p.tkr_minor.revision_simple ,
                         p.tkr_minor.other_trt ,
                         p.tkr_minor.tkr_minor,
                         #p.tkr_serious.tkr_minor,
                         p.tkr_serious.revision_simple,
                         p.tkr_serious.revision_complex ,
                         p.tkr_serious.other_trt ,
                         p.tkr_normal_after_revision.revision_simple,
                         p.tkr_normal_after_revision.revision_complex ,
                         p.tkr_normal_after_revision.tkr_normal_after_revision ,
                         p.revision_simple.tkr_minor ,
                         p.revision_simple.tkr_serious ,
                         p.revision_simple.tkr_normal_after_revision ,
                         #p.revision_simple.other_trt,
                         p.revision_complex.tkr_serious ,
                         p.revision_complex.tkr_normal_after_revision ,
                         #p.other_trt.tkr_minor,p.other_trt.tkr_serious <- list.params[1]
                         p.other_trt.tkr_normal_after_revision ,
                         #p.other_trt.revision_simple,
                         p.tkr_normal.death,
                         p.tkr_normal_after_revision.death,
                         p.death,
                         c.tkr ,
                         c.tkr_normal ,
                         c.tkr_minor,
                         c.tkr_serious ,
                         c.tkr_normal_after_revision ,
                         c.revision_simple ,
                         c.revision_complex ,
                         c.other_trt,
                         c.death ,
                         u.tkr ,
                         u.tkr_normal ,
                         u.tkr_minor ,
                         u.tkr_serious ,
                         u.tkr_normal_after_revision ,
                         u.revision_simple ,
                         u.revision_complex ,
                         u.other_trt
    )
    
    params_pdf.pbox <- list_test.pbox
    # params_pdf.MC <- list_test.MC
    ##############################################################################
    #running PBA
    ################################################################################
    #source("TAH_functions_new_model_v4.R") #load 
    

    if(pbox.parallel==0){ #no parallel over sub-intervals (K set)
    #running the optimization
    result_list_outer.sp.pbox <- fn_bbmod1_simp_average_inmb_bobyqa(params_pdf.pbox,
                                                                    disc ,
                                                                    num.states,
                                                                    time.horizon,
                                                                    init,
                                                                    n,
                                                                    WTP,
                                                                    index_slice=input_data,pbox.int=pbox.int,index.array=index.array)

    } else {
    
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
    #running the optimization
    result_list_outer.sp.pbox <- future.apply::future_lapply(list.index, fun_pba_par ,params_pdf.pbox,
                                                                    disc ,
                                                                    num.states,
                                                                    time.horizon,
                                                                    init,
                                                                    n,
                                                                    WTP,
                                                                    index_slice=input_data,pbox.int=pbox.int,index.array=index.array)
    }

    
    #approach 2 MC only + average values (PSA)
    if(sp.yes==1){
      if(i.int==1 & i.stat==1){
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
        MC.result <- future.apply::future_lapply(params_MC, TKR_model_MC, 
                                                 c.CAS ,
                                                 p.tkr_init.tkr_normal ,
                                                 p.tkr_init.tkr_minor ,
                                                 p.tkr_init.tkr_serious,
                                                 p.tkr_normal.tkr_minor,
                                                 p.tkr_normal.tkr_serious ,
                                                 p.tkr_normal.tkr_normal ,
                                                 p.tkr_minor.tkr_normal ,
                                                 #p.tkr_minor.tkr_serious, 
                                                 p.tkr_minor.revision_simple ,
                                                 p.tkr_minor.other_trt ,
                                                 p.tkr_minor.tkr_minor,
                                                 #p.tkr_serious.tkr_minor,
                                                 p.tkr_serious.revision_simple,
                                                 p.tkr_serious.revision_complex ,
                                                 p.tkr_serious.other_trt ,
                                                 p.tkr_normal_after_revision.revision_simple,
                                                 p.tkr_normal_after_revision.revision_complex ,
                                                 p.tkr_normal_after_revision.tkr_normal_after_revision ,
                                                 p.revision_simple.tkr_minor ,
                                                 p.revision_simple.tkr_serious ,
                                                 p.revision_simple.tkr_normal_after_revision ,
                                                 #p.revision_simple.other_trt,
                                                 p.revision_complex.tkr_serious ,
                                                 p.revision_complex.tkr_normal_after_revision ,
                                                 #p.other_trt.tkr_minor,p.other_trt.tkr_serious <- list.params[1]
                                                 p.other_trt.tkr_normal_after_revision ,
                                                 #p.other_trt.revision_simple,
                                                 p.tkr_normal.death,
                                                 p.tkr_normal_after_revision.death,
                                                 p.death,
                                                 c.tkr ,
                                                 c.tkr_normal ,
                                                 c.tkr_minor,
                                                 c.tkr_serious ,
                                                 c.tkr_normal_after_revision ,
                                                 c.revision_simple ,
                                                 c.revision_complex ,
                                                 c.other_trt,
                                                 c.death ,
                                                 u.tkr ,
                                                 u.tkr_normal ,
                                                 u.tkr_minor ,
                                                 u.tkr_serious ,
                                                 u.tkr_normal_after_revision ,
                                                 u.revision_simple ,
                                                 u.revision_complex ,
                                                 u.other_trt,
                                                 disc ,
                                                 num.states,
                                                 time.horizon,
                                                 init,
                                                 n,
                                                 WTP)
      }
    }

    #save.image(file=paste(currpath,"/Rdata/Rdata_TAH_PBA_statistics_15",statistics,"_bound_" ,bound,"_npbox_",num.pbox.int ,"_" , date(), "..RData",sep=""))
    
    ################################################################################
    #processing PBA outputs
    nmb <- sapply (result_list_outer.sp.pbox, function (x) {length (x) <- 3; return (x)})
    nmb.t <- t(nmb)
    #INMB approach (as output) 16 May 2021
    df.nmb.lower <- data.frame(nmb.t[,1])
    df.nmb.upper <- data.frame(abs(nmb.t[,2]))
    df.weight <- data.frame(nmb.t[,3]) #weight
    
    cdf.x.upper.sp.pbox <- array(seq(0,n.x), c(n.x))
    cdf.x.lower.sp.pbox <- array(seq(0,n.x), c(n.x))
    
    #calculate empirical cdf
    for (j in 1:n.x){
      #cdf.x.lower.sp[j]<- sum(df.x$V2[df.x$V1<=x.values[j]]) #note reverse, see Ferson 2015
      cdf.x.lower.sp.pbox[j]<- sum(df.weight[df.nmb.upper<=x.values[j]]) #note reverse, see Ferson 2015
      cdf.x.upper.sp.pbox[j]<- sum(df.weight[df.nmb.lower<=x.values[j]])
    }
    
    x.mean.lower.sp.df.pbox <- as.data.frame(cbind(x.values,cdf.x.lower.sp.pbox))
    colnames(x.mean.lower.sp.df.pbox) <- c("x","y")
    
    x.mean.upper.sp.df.pbox <- as.data.frame(cbind(x.values,cdf.x.upper.sp.pbox))
    colnames(x.mean.upper.sp.df.pbox) <- c("x","y")
    
    #saving results

    cdf.all[,1,i.stat,i.int] <-cdf.x.lower.sp.pbox
    cdf.all[,2,i.stat,i.int] <-cdf.x.upper.sp.pbox
    
  }
}

########################################################################
save.image(file=paste(currpath,"/rdata/manuscript/rdata_tkr_pbox_",pbox.num, "_",statistics ,"_", date(),"_rrCASminmax_", rr.CAS.min ,"_"  , rr.CAS.max, ".RData",sep=""))

########################################################################
#vary i.stat only, not i.int
# pbox1.lower.psa <- as.data.frame(cbind(x.values,cdf.all[,1,1,1]))
# pbox1.upper.psa <- as.data.frame(cbind(x.values,cdf.all[,2,1,1]))

pbox1.lower <- as.data.frame(cbind(x.values,cdf.all[,1,1,1]))
pbox1.upper <- as.data.frame(cbind(x.values,cdf.all[,2,1,1]))

c.color_2 <- c(c0_precise="#ce1256",
               c1_pbox1="#7E7E7E")

x.breaks <- round_any(seq(min.x,900000,200000),100000)

wt <- c(rep(1/n.sims,n.sims))
df.x <- as.data.frame(cbind(unlist(MC.result),wt)) 
colnames(df.x) <- c("V1","V2")
cdf.x.lower.basecase <- array(seq(0,n.x), c(n.x))

for (j in 1:n.x){
  cdf.x.lower.basecase[j]<- sum(df.x$V2[df.x$V1<=x.values[j]]) 
}
x.mean.lower.basecase <- as.data.frame(cbind(x.values,cdf.x.lower.basecase))
colnames(x.mean.lower.basecase) <- c("x.values","V2")
colnames(pbox1.lower.mean) <- c("x.values","V2")
colnames(pbox1.upper.mean) <- c("x.values","V2")


ltype <- "solid"
lsize <- 1
#ltype <- "twodash"
########################################################################
#plot
plot1 <- ggplot()+  geom_line(data = x.mean.lower.basecase, aes(x=x.values, y=V2, color = "c0_precise"), alpha=1, size=lsize, linetype="solid") +
  geom_line(data = pbox1.upper, aes(x=x.values, y=V2, color = "c1_pbox1"), alpha=1, size=lsize,linetype=ltype) +
  geom_line(data = pbox1.lower, aes(x=x.values, y=V2, color = "c1_pbox1"), alpha=1, size=lsize,linetype=ltype) +
   geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3)+
  scale_x_continuous(name = "Incremental net monetary benefit",breaks=x.breaks,expand=c(0,0),limits=c(min.x,900000),labels=function(x) format(x, scientific = TRUE)) +
  scale_y_continuous(name = "CDF",expand=c(0,0),limits=c(0, 1.01)) +
  scale_color_manual(values = c.color_2, labels=c("Precise","P-box (mean)","P-box (min+max)" ),
                     name = "CDF Types")+
  # scale_color_manual(values = c.color_6, labels=c("Precise","P-box (3)","P-box (4)","P-box (5)","P-box (10)","P-box (50)"),
  #                    name = "CDF Types")+
  #ggtitle("P-box of quantity of interest (QoI)" ) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(color = "black", size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(color = "black", size = 8), legend.position = "bottom") +
  guides(color = guide_legend(title="CDF types"), linetype = guide_legend(title="CDF types"))
ggsave(paste(currpath,"/plot/manuscript/tkr_pbox_", pbox.num, "_1000_",statistics, "_",date(), "_rrCASminmax_", rr.CAS.min ,"_"  , rr.CAS.max,".tiff",sep="") , units="in", width=7, height=4, dpi=200, compression = 'lzw')


time1 <- (proc.time()-ptm1)
ela1 <-time1[3]/60
print(paste("time to run PBA:  ",ela1 ," minutes"))
########################################################################
#end of code
########################################################################