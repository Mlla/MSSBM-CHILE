#-------------------------------------------------------------------------------
#
# Final calibration file (all functions included) for
#   Blanchard et al. JAppliedEcol
#
#-------------------------------------------------------------------------------

rm(list=ls())
library(minpack.lm)

# - Set working directory of Dropbox folder
# setwd("/Users/jlblanch/Dropbox/sizeBasedStuff/NorthSea")
# setwd("D:/Repository/Dropbox/sizeBasedStuff/NorthSea")
# setwd("/Users/Mlla/Dropbox/Size Spectrum Model Examples/3. Multispecies size spectrum model/multispecies_sizebasedmodels/NorthSea/R")

# Mlla: The files for setup and project the model are all in the same folder now
source("paramNorthSeaModel.R")
source("calibration_funcs.r")
source("SizeBasedModel.r")
source("SelectivityFuncs.r")
source("plots.r")

#----------------------------------------------
# 1) Load the parameter input data
#----------------------------------------------

#- interaction matrix
load("interactionmatrix_Schoener_twostage2D.RData")

#- parameter file
# paramfile <- "IMAGE_paper_params_new.txt" # Mlla: Here I did not find this particular file within your folder, so 
# I tried to run the code with a similar: ""IMAGE_paper_params_newNH.txt"
paramfile <- "IMAGE_paper_params_newNH.txt"

#----------------------------------------------
# 2) Setup the parameter file
#----------------------------------------------

param <-  paramNSModel(fileSpecies  = paramfile,
                       theta        = theta,
                       kap          = 1e11)
param$species$eRepro    <- 1
param$tmax              <- 1000
param$dt                <- 1

#- Get Rmax values
# source("R/getRmax.R") #Mlla: I commented this function becuase I did not find it. So instead I used a vector or Rmax from 
# North sea data in MIZER examples. No quiet sure if the type of object is the same or will cause some problems later.
# newRmax are in log scale
# param$species$R0        <- Rmax$RmaxAB

param$species$R0 <- c(26.52, 25.92, 24.22, 24.77, 22.41, 21.60, 19.41, 20.34, 22.45, 21.02, 17.73, 19.64)

#----------------------------------------------
# 3) Setup initial model
#----------------------------------------------

Initialcomm             <- Setup(param)
Initialcomm             <- Project(Initialcomm)
plotResults(Initialcomm)

#- Check that initial run is at equib and has SSB less than SSBmax
Nhat                    <- Initialcomm$N[dim(Initialcomm$N)[1],,]
SSBhat                  <- apply(sweep(Initialcomm$psi * Nhat,2,Initialcomm$w * Initialcomm$dw,"*"),1,sum) /1e6
if(any(SSBhat > (Initialcomm$param$species$maxSSB))) stop("SSB greater than SSBmax")

#----------------------------------------------
# 4) Setup optimization criteria
#----------------------------------------------

meantsteps              <- 100
maxSSBpenalty           <- FALSE
logerror                <- TRUE
Initialcomm$param$tmax  <- 300
lower                   <- rep(1,13)
upper                   <- rep(100,13)
logParams               <- c(log(Initialcomm$param$species$R0),log(param$kap))

#----------------------------------------------
# 5) Optimization call
#----------------------------------------------
count                 <- 0 #reset counting
# We modified the following lines becuase Julia used optim in the North Sea work.

# Mlla: the function that is call by optim function "calibrateR0_yield_SSB_kappa_optim", it is not in the folder calibration_func
# but in the line 116. 

optim_yield_ssb_kappa <- optim(par=logParams,fn = calibrateR0_yield_SSB_kappa_optim,
                                method="L-BFGS-B",hessian=F,
                                lower=c(rep(15,12),20),upper=c(rep(35,12),30),
                                Initialcomm = Initialcomm, meantsteps = meantsteps,
                                maxSSBpenalty=maxSSBpenalty, logerror = logerror,
                                paramfile=paramfile)

# Mlla: I am having an error in the optimization, and I am not able to create the output from optim function.

# nlslm_yield_ssb_kappa <- nls.lm(par=logParams,fn = calibrateR0_yield_SSB_kappa_nls,
#                                 lower=c(rep(15,12),20),upper=c(rep(35,12),30),
#                                 Initialcomm = Initialcomm, meantsteps = meantsteps,
#                                 maxSSBpenalty=maxSSBpenalty, logerror = logerror,
#                                 paramfile=paramfile)

# Nonlinear regression via the Levenberg-Marquardt algorithm
# parameter estimates: 27.5837024746905, 27.1561754623063, 30.7349325423871, 27.3123297759941, 25.1961943832892, 25.7774081810535, 23.0108970612388, 31.1339477152861, 28.4660352897801, 26.2304390655457, 23.4165712287494, 27.8576277457183, 25.899867253659 
# residual sum-of-squares: 4.058
# reason terminated: Relative error between `par' and the solution is at most `ptol'.


#----------------------------------------------
# 6) Function calls
#----------------------------------------------

#---------------------
#- Calibration routine optim
#---------------------

calibrateR0_yield_SSB_kappa_optim <- function(logR0hat, Initialcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE,paramfile="input/IMAGE_paper_params_new.txt")
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)
  
  kap <- exp(logR0hat[length(logR0hat)])
  cat("count: ", count, ". logR0", signif(logR0hat,5), "\n")
  
  Newparam<-paramNSModel(fileSpecies=paramfile,
                         theta = theta,
                         kap=kap)
  
  Newparam$species$R0 <- exp(logR0hat[-length(logR0hat)])
  
  Newparam$theta      <- Initialcomm$param$theta
  Newparam$tmax       <- Initialcomm$param$tmax
  Newparam$dt         <- Initialcomm$param$dt
  Newparam$wPPcut     <- Initialcomm$param$wPPcut
  
  # Project model with the new R0, continuing from where the initial run left off
  Newcomm             <- Setup(Newparam, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm             <- Project(Newcomm)
  ye                  <- Yielderror(Newcomm,meantsteps=meantsteps, logerror = logerror)
  se                  <- SSBerror(Newcomm,meantsteps=meantsteps, maxSSBpenalty=maxSSBpenalty, logerror = logerror)
  error <- c(se,ye)
  # error <- sum(sqrt(error[!is.na(error)]^2)) #We need to sum here already, because else we add 1e2 to all the 21 parameters
  
  
  # Do we need to apply an extinction penalty?
  extinct             <- Extinct_test(Newcomm,extinct_nsteps = 10)
  if(extinct)
  {
    error             <- error + 1e1
    cat("Extinction!\n")
  }
  
  # Project model without fishing and apply an extinction penalty using new parameters and where the new model left off
  Newparam$Q[,]       <- 0
  Newcomm             <- Setup(Newparam, ContinueCalculation=T, initialcommunity=Newcomm)
  Newcomm             <- Project(Newcomm)
  
  extinct2            <- Extinct_test(Newcomm,extinct_nsteps = 10)
  if(extinct2)
  {
    error             <- error + 1e1
    cat("Extinction!\n")
  }
  
  cat("error: ", error, "\n\n")
  return(error)
}

#---------------------
#- Calibration routine optim #Mlla: I commented this function because is not used eaither
#---------------------

# calibrateR0_yield_SSB_kappa_nls <- function(logR0hat=logParams, Initialcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE,paramfile="input/IMAGE_paper_params_new.txt")
# {
#   #count <- count + 1
#   assign("count", count+1, pos = .GlobalEnv)
#   
#   kap <- exp(logR0hat[length(logR0hat)])
#   cat("count: ", count, ". logR0", signif(logR0hat,5), "\n")
#   
#   Newparam<-paramNSModel(fileSpecies=paramfile,
#                          theta = theta,
#                          kap=kap)
#   
#   Newparam$species$R0 <- exp(logR0hat[-length(logR0hat)])
#   
#   Newparam$theta      <- Initialcomm$param$theta
#   Newparam$tmax       <- Initialcomm$param$tmax
#   Newparam$dt         <- Initialcomm$param$dt
#   Newparam$wPPcut     <- Initialcomm$param$wPPcut
#   
#   # Project model with the new R0, continuing from where the initial run left off
#   Newcomm             <- Setup(Newparam, ContinueCalculation=T, initialcommunity=Initialcomm)
#   Newcomm             <- Project(Newcomm)
#   ye                  <- Yielderror(Newcomm,meantsteps=meantsteps, logerror = logerror)
#   se                  <- SSBerror(Newcomm,meantsteps=meantsteps, maxSSBpenalty=maxSSBpenalty, logerror = logerror)
#   error <- c(se,ye)
#   error <- error[!is.na(error)]
#   # error <- sqrt(error[!is.na(error)]^2) #Make sure the errors are positive
#   #Don't think we need this, nls.lm takes a vector of residuals and then computes the sum of squared errors 
#   
#   # Do we need to apply an extinction penalty?
#   extinct             <- Extinct_test(Newcomm,extinct_nsteps = 10)
#   if(extinct)
#   {
#     error          <- error * 1e1
#     cat("Extinction!\n")
#   }
#   
#   # Project model without fishing and apply an extinction penalty using new parameters and where the new model left off
#   Newparam$Q[,]       <- 0
#   Newcomm             <- Setup(Newparam, ContinueCalculation=T, initialcommunity=Newcomm)
#   Newcomm             <- Project(Newcomm)
#   
#   extinct2            <- Extinct_test(Newcomm,extinct_nsteps = 10)
#   if(extinct2)
#   {
#     error          <- error * 1e1
#     cat("Extinction!\n")
#   }
#   
#   cat("error: ", sum(error^2), "\n\n")
#   return(error)
# }

#---------------------
#- Extinction routine
#---------------------
Extinct_test <- function(Newcomm, extinct_nsteps=10, extinct_threshold=0.1,ma_nsteps=20)
{
  Biomass <- rowSums(sweep(Newcomm$N,3,Newcomm$w*Newcomm$dw,"*"),dims=2)
  
  ma <- apply(Biomass, 2, function(x) filter(x, rep(1/ma_nsteps,ma_nsteps), sides=1))
  
  # Get biomass relative to last time point
  relma <- sweep(ma, 2, ma[dim(ma)[1],], "/")
  
  # Is biomass nsteps ago some % higher
  extinct <- any(relma[dim(relma)[1]-extinct_nsteps+1,] > (1+extinct_threshold))
  
  return(extinct)
}

#---------------------
#- SSB error
#---------------------
SSBerror <- function(Newcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE)
{
  if (is.na(meantsteps))
  {
    Nhat <- Newcomm$N[dim(Newcomm$N)[1],,]
    SSBhat <- apply(sweep(Newcomm$psi * Nhat,2,Newcomm$w * Newcomm$dw,"*"),1,sum) /1e6
  }
  else
  {
    Nhat    <- apply(Newcomm$N[(dim(Newcomm$N)[1]-meantsteps+1):dim(Newcomm$N)[1],,],c(2,3),mean)
    SSBhat <- apply(sweep(Newcomm$psi * Nhat,2,Newcomm$w * Newcomm$dw,"*"),1,sum) /1e6
  }
  
  SSBobs <- Newcomm$param$species$SSB_8595
  if (!logerror)
    SSBerr <- (SSBhat / SSBobs) - 1
  if (logerror)
    SSBerr <- log(SSBhat / SSBobs)
  
  # Apply penalty if predicted SSB is greater than maximum. Maximum taken from Rochet et al
  if (maxSSBpenalty)
  {
    SSBerr[SSBhat > Newcomm$param$species$maxSSB] <- 1e2
    if (any(SSBhat > Newcomm$param$species$maxSSB))
      cat("SSB maxed out\n")
  }
  
  return(SSBerr)
}

#---------------------
#- Yield error
#---------------------
Yielderror <- function(Newcomm, meantsteps= NA, logerror = FALSE)
{
  # Yield and N - final value or mean of last tsteps?
  if (is.na(meantsteps))
    Yieldhat <- Newcomm$Yield[dim(Newcomm$Yield)[1],]/1e6
  else
    Yieldhat <- apply(Newcomm$Yield[(dim(Newcomm$Yield)[1]-meantsteps+1):dim(Newcomm$Yield)[1],],2,mean)/1e6
  
  # Use mean Catch 1985 - 1995
  Yieldobs <- Newcomm$param$species$Catch_8595
  
  if (!logerror)
    Yielderr <- (Yieldhat / Yieldobs) - 1
  if (logerror)
    Yielderr <- log(Yieldhat) - log(Yieldobs)
  
  return(Yielderr)
}
