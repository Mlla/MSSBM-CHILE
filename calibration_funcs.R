# Mlla's copy from Dr Julia Blanchard 02/10/2017
# ==============================================


# Calibration

# Functions for calibration including:

# Error functions:

# Yielderror - sets up and projects a new model from an existing commmunity
#               - returns a vector of errors based on observed and actual yields

# YieldSSBerror - sets up and projects a new model from an existing commmunity
#               - returns a vector of errors based on observed and actual yields, and
#               - SSB


# Objective functions:
# 
# calibrateR0 - calls YieldSSBerror with new values of R0, returns sum of squared errors
#             - used with optim

# calibrateR0.nlslm - as calibrateR0 but returns a vector of errors. This is because
#                   - nls.lm takes a vector of errors , not the sum of squared errors.


#**************************** Wrapper functions ********************************
count <- 0

calibrateR0_yield <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)
  Initialcomm$param$species$R0 <- exp(logR0hat)
  cat("count: ", count, ". logR0", signif(logR0hat,3), "\n")
  
  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  error <- sum((Yielderror(Newcomm,meantsteps=meantsteps)^2), na.rm=T)
  
  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }
  
  cat("error: ", error, "\n\n")
  return(error)
}

calibrateR0_SSB <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  assign("count", count+1, pos = .GlobalEnv)
  Initialcomm$param$species$R0 <- exp(logR0hat)
  cat("count: ", count, ". logR0", signif(logR0hat,3), "\n")

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  error <- sum((SSBerror(Newcomm,meantsteps=meantsteps)^2), na.rm=T)

  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", error, "\n\n")
  return(error)
}


calibrateR0_survey <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)
  Initialcomm$param$species$R0 <- exp(logR0hat)
  cat("count: ", count, ". logR0", signif(logR0hat,3), "\n")

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  error <- sum((Surveyerror(Newcomm,meantsteps=meantsteps)^2), na.rm=T)
  
  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", error, "\n\n")
  return(error)
}

calibrateR0_yield_SSB <- function(logR0hat, Initialcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE, fixedFeedingLevel = NA, fixedPredationMortality = NA)
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)
  Initialcomm$param$species$R0 <- exp(logR0hat)
  cat("count: ", count, ". logR0", signif(logR0hat,5), "\n")

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm, fixedFeedingLevel = fixedFeedingLevel, fixedPredationMortality = fixedPredationMortality)
  ye <- Yielderror(Newcomm,meantsteps=meantsteps, logerror = logerror)
  se <- SSBerror(Newcomm,meantsteps=meantsteps, maxSSBpenalty=maxSSBpenalty, logerror = logerror)
  error <- sum((c(se,ye)^2), na.rm=T)

  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", error, "\n\n")
  return(error)
}

calibrateR0_eRepro_yield_SSB <- function(logParams, Initialcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE, fixedFeedingLevel = NA, fixedPredationMortality = NA)
{
    #count <- count + 1
    assign("count", count+1, pos = .GlobalEnv)
    Initialcomm$param$species$R0 <- exp(logParams[1:(length(logParams)-1)])
    Initialcomm$param$species$eRepro <- exp(logParams[length(logParams)])
    cat("count: ", count, "\n")
    cat("logR0", signif(logParams[1:(length(logParams)-1)],5), "\n")
    cat("eRepro", signif(exp(logParams[length(logParams)]),5), "\n")

    # Project model with the new R0, continuing from where the initial run left off
    Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
    Newcomm <- Project(Newcomm, fixedFeedingLevel = fixedFeedingLevel, fixedPredationMortality = fixedPredationMortality)
    ye <- Yielderror(Newcomm,meantsteps=meantsteps, logerror = logerror)
    se <- SSBerror(Newcomm,meantsteps=meantsteps, maxSSBpenalty=maxSSBpenalty, logerror = logerror)
    error <- sum((c(se,ye)^2), na.rm=T)

  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", error, "\n\n")
  return(error)
}

calibrateR0_yield_SSB_kappa <- function(logR0hat, Initialcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE)
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)

  kap <- exp(logR0hat[length(logR0hat)])
  cat("count: ", count, ". logR0", signif(logR0hat,5), "\n")

  Newparam<-paramNSModel(fileSpecies="input/IMAGE_paper_params_noSaithe.csv",
                          theta = theta,
                          kap=kap)

  Newparam$species$R0 <- exp(logR0hat[-length(logR0hat)])

  Newparam$tmax <- 75
  Newparam$dt <- 1

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Newparam, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  ye <- Yielderror(Newcomm,meantsteps=meantsteps, logerror = logerror)
  se <- SSBerror(Newcomm,meantsteps=meantsteps, maxSSBpenalty=maxSSBpenalty, logerror = logerror)
  error <- sum((c(se,ye)^2), na.rm=T)

  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", error, "\n\n")
  return(error)
}



calibrateR0_yield_survey <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)
  Initialcomm$param$species$R0 <- exp(logR0hat)
  cat("count: ", count, ". logR0", signif(logR0hat,3), "\n")

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  ye <- Yielderror(Newcomm,meantsteps=meantsteps)
  se <- Surveyerror(Newcomm,meantsteps=meantsteps)
  error <- sum((c(se,ye)^2), na.rm=T)

  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", error, "\n\n")
  return(error)
}


calibrateR0_yield.nlslm <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  Initialcomm$param$species$R0 <- exp(logR0hat)

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  return(Yielderror(Newcomm,meantsteps=meantsteps))
}

calibrateR0eRepro_yield_blim <- function(logparams, Initialcomm, meantsteps=NA)
{
  Initialcomm$param$species$R0 <- exp(logparams[1:Initialcomm$param$nspp])
  Initialcomm$param$species$eRepro <- exp(logparams[(Initialcomm$param$nspp+1):(Initialcomm$param$nspp*2)])  

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  return(sum((YieldBlimerror(Newcomm,meantsteps=meantsteps)^2)))
}

calibrateR0eRepro_yield_blim.nlslm <- function(logparams, Initialcomm, meantsteps=NA)
{
  Initialcomm$param$species$R0 <- exp(logparams[1:Initialcomm$param$nspp])
  Initialcomm$param$species$eRepro <- exp(logparams[(Initialcomm$param$nspp+1):(Initialcomm$param$nspp*2)])  

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  return(YieldBlimerror(Newcomm,meantsteps=meantsteps))
}

#*******************************************************************************
# nls.lm() wrapper functions

calibrateR0_yield_SSB_nlslm <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  #count <- count + 1
  assign("count", count+1, pos = .GlobalEnv)
  Initialcomm$param$species$R0 <- exp(logR0hat)
  cat("count: ", count, ". logR0", signif(logR0hat,3), "\n")

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  ye <- Yielderror(Newcomm,meantsteps=meantsteps)
  se <- SSBerror(Newcomm,meantsteps=meantsteps)
  #error <- sum((c(se,ye)^2), na.rm=T)
  error <- c(se,ye)
  error <- error[!is.na(error)]

  # Do we need to apply an extinction penalty?
  extinct <- Extinct_test(Newcomm)
  if(extinct)
  {
    error <- error + 1e9
    cat("Extinction!\n")
  }

  cat("error: ", sum(error^2, na.rm=T), "\n\n")
  return(error)
}



# ************ error functions *************************************************

# Function to test if anything is going extinct
# Look at mean total biomass of last nsteps
# Compare to mean total biomass of last (nsteps*2):nsteps
# if less assume that species is on the way out
# So (default) if a species loses 50% of its biomass in the last 10 time steps - it's doomed
# Should be wide enough margin to cope with wobbling caused by solver and stocks 'settling down'
Extinct_test <- function(Newcomm, extinct_nsteps=10, extinct_threshold=0.5,ma_nsteps=10)
{
  #browser()

  Biomass <- rowSums(sweep(Newcomm$N,3,Newcomm$w*Newcomm$dw,"*"),dims=2)

  # Get biomass relative to final timestep
  #relB <- sweep(Biomass, 2, Biomass[dim(Biomass)[1],], "/")
  #extinct <- any(relB[dim(relB)[1]-extinct_nsteps+1,] > (1+extinct_threshold))
  # Need to smooth in some way - else crank up extinct_threshold

  # I tried using a moving average to help sort out wrinkles caused by solution
  # But can cause more wrinkles due to periodic affects -
  ma <- apply(Biomass, 2, function(x) filter(x, rep(1/ma_nsteps,ma_nsteps), sides=1))
  # Get biomass relative to last time point
  relma <- sweep(ma, 2, ma[dim(ma)[1],], "/")
  # Is biomass nsteps ago some % higher
  extinct <- any(relma[dim(relma)[1]-extinct_nsteps+1,] > (1+extinct_threshold))

  return(extinct)
}

SSBerror <- function(Newcomm, meantsteps=NA, maxSSBpenalty=FALSE, logerror = FALSE)
{
 #browser()
 if (is.na(meantsteps))
  {
    Nhat <- Newcomm$N[dim(Newcomm$N)[1],,]
    SSBhat <- apply(sweep(Newcomm$psi * Nhat,2,Newcomm$w * Newcomm$dw,"*"),1,sum) /1e6
  }
  else
  {
    Nhat <- apply(Newcomm$N[(dim(Newcomm$N)[1]-meantsteps+1):dim(Newcomm$N)[1],,],c(2,3),mean)
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
    #cat(SSBhat > Newcomm$param$species$maxSSB, "\n")
    SSBerr[SSBhat > Newcomm$param$species$maxSSB] <- 1e9
    if (any(SSBhat > Newcomm$param$species$maxSSB))
      cat("SSB maxed out\n")
  }
  
  return(SSBerr)
}

Yielderror <- function(Newcomm, meantsteps= NA, logerror = FALSE)
{
  # Yield and N - final value or mean of last tsteps?
  # Note the scaling
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
#  Yielderr <- (Yieldhat / Yieldobs) - 1
  #cat("Yielderr: ", Yielderr, "\n")
  return(Yielderr)
}

Surveyerror <- function(Newcomm, meantsteps= NA)
{
  # N - final value or mean of last tsteps?
  if (is.na(meantsteps))
    Nhat <- Newcomm$N[dim(Newcomm$N)[1],,]
  else
    Nhat <- apply(Newcomm$N[(dim(Newcomm$N)[1]-meantsteps+1):dim(Newcomm$N)[1],,],c(2,3),mean)

    # Do one species at a time
    Nhatw <- rep(NA,nrow(Newcomm$param$species))
    #i <- 1
    for (i in 1:nrow(Newcomm$param$species))
    {
      if (is.na(Newcomm$param$species$Nden[i])) Nhatw[i] <- NA
      else
      {
        # get the two end size bins which only make a partial contribution
        startw <- max(which((Newcomm$w - Newcomm$param$species$lower_w[i])<0))
        endw <- min(which((Newcomm$w - Newcomm$param$species$upper_w[i])>0)) - 1
        # Get size bins that should be included
        welem <- rep(0,length(Newcomm$w))
        welem[startw:endw] <- 1
        # Get contribution of each bin
        wcontr <- welem
        # Set proportion contribution of the end bins
        wcontr[startw] <- (Newcomm$w[startw+1] - Newcomm$param$species$lower_w[i]) / Newcomm$dw[startw]
        wcontr[endw] <- (Newcomm$param$species$upper_w[i] - Newcomm$w[endw]) / Newcomm$dw[endw]
        Nhatw[i] <- sum((Nhat[i,] * wcontr) * (Newcomm$dw * welem))
      }
    }

    survey_error <- (Nhatw / Newcomm$param$species$Nden) - 1
    #cat("survey_error: ", survey_error, "\n")
    return(survey_error)
}

YieldBlimerror <- function(Newcomm, meantsteps= NA)
{
  # Yield and N - final value or mean of last tsteps?
  if (is.na(meantsteps))
  {
    Yieldhat <- Newcomm$Yield[dim(Newcomm$Yield)[1],]/1e6
    Nhat <- Newcomm$N[dim(Newcomm$N)[1],,]
    SSBhat <- apply(sweep(Newcomm$psi * Nhat,2,Newcomm$w * Newcomm$dw,"*"),1,sum)
    RDI <- Newcomm$RDI[dim(Newcomm$RDI)[1],]
  }
  else
  {
    Yieldhat <- apply(Newcomm$Yield[(dim(Newcomm$Yield)[1]-meantsteps+1):dim(Newcomm$Yield)[1],],2,mean)/1e6
    Nhat <- apply(Newcomm$N[(dim(Newcomm$N)[1]-meantsteps+1):dim(Newcomm$N)[1],,],c(2,3),mean)
    SSBhat <- apply(sweep(Newcomm$psi * Nhat,2,Newcomm$w * Newcomm$dw,"*"),1,sum)
    RDI <- apply(Newcomm$RDI[(dim(Newcomm$RDI)[1]-meantsteps+1):dim(Newcomm$RDI)[1],],2,mean)
  }
  # Use mean Catch 1985 - 1995
  Yieldobs <- Newcomm$param$species$Catch_8595
  Yielderr <- (Yieldhat / Yieldobs) - 1

  # Blim error
  # Calc SSB
#  SSBm <- Newcomm$psi * e
#  SSBhat <- SSBm * N %*% dw
#  SSBhat <- (Newcomm$psi * Nhat) %*% (Newcomm$w*Newcomm$dw)
  SSBerr <- ((SSBhat * Newcomm$param$species$R0) / (Newcomm$param$species$Blim * 1e6 * RDI)) - 1
  # S.R = Rsave = R = = 0.5*param.rho(i)*SSB/w(1);
  # in R code is RDI, i.e. the density independent recruitment
  error <- c(Yielderr, SSBerr)
  
  return(error)
}



#************************************************************************
# Calibration plots

plotYieldCalib <- function(model, meantsteps = 10, yieldname="Catch_8595",ylim=c(NA, NA))
{
  Yobs <- model$param$species[,yieldname]
  endt <- dim(model$Yield)[1]
  Yhat <- apply(model$Yield[(endt-meantsteps+1):endt,],2,mean)*1e-6

#cat("Yhat / Yobs", Yhat / Yobs, "\n")

  if(is.na(ylim[1])) ylim[1] <- min((Yhat/Yobs)-1)
  if(is.na(ylim[2])) ylim[2] <- max((Yhat/Yobs)-1)
  plot(x=1:12, y=(Yhat / Yobs)-1, type="n", xlab="Species",ylab="(Yhat / Yobs) - 1", ylim=ylim, xaxt="n")
  axis(1,1:12,model$param$species$species)
  points(x=1:12,y=(Yhat / Yobs)-1)
  lines(x=c(1,12),y=c(0,0),lty=2)
}

#plotYieldCalib(NSmodel,yieldname="Catch_8595")

ploteReproCalib <- function(model, meantsteps=10, ylim=c(NA,NA))
{
  Nhat <- apply(model$N[(dim(model$N)[1]-meantsteps+1):dim(model$N)[1],,],c(2,3),mean)
  SSBhat <- apply(sweep(model$psi * Nhat,2,model$w * model$dw,"*"),1,sum)
  RDI <- apply(model$RDI[(dim(model$RDI)[1]-meantsteps+1):dim(model$RDI)[1],],2,mean)
#  SSBerr <- ((SSBhat * model$param$species$R0) / (model$param$species$Blim * 1e6 * RDI))

  SSBerr <- ((SSBhat * model$param$species$R0) / (model$param$species$Blim * 1e6 * RDI)) - 1

  if(is.na(ylim[1])) ylim[1] <- min(SSBerr)
  if(is.na(ylim[2])) ylim[2] <- max(SSBerr)


  plot(x=1:12,y=SSBerr,ylim=ylim,xlab = "Species", ylab="SSB error", xaxt="n")
  axis(1,1:12,model$param$species$species)
  points(x=1:12,y=SSBerr)
  lines(x=c(1,12),y=c(0,0),lty=2)
}

#ploteReproCalib(NSmodel,meantsteps=10)

#*********************************************************************************
# Check calibration routine is OK
# e.g. check why RSS in routine is different to final
# Hessian? BFGS?


# Profile plots
# With so many species we can only do a series of 1D profiles, rather than
# joint profiles (suppose we could do a series of 2D profiles)

# Two parameters per species
profile1D <- function(model,sp,param,error = "YieldSSB", range=1/3, length=11)
{
    # Use species sp
    # Set up a parameter range of +- range, length
    # Get the RSS
    orig_param <- model$param$species[sp,param]
    param_seq <- seq(from=orig_param*(1-range), to = orig_param*(1+range), length=11)
    RSS_store <- rep(NA,length)
    for (i in 1:length)
    {
	#cat("iter: ", i, "\n")
	model$param$species[sp,param] <- param_seq[i]
	if (error == "YieldSSB")
	    RSS_store[i] <- sum(YieldSSBerror(model,meantsteps=10)^2)
	if (error == "Yield")
	    RSS_store[i] <- sum(Yielderror(model,meantsteps=10)^2)
    }
    return(list(param=param_seq,RSS=RSS_store))
}


plotCalib <- function(model,param,prof_dat = NULL,error="YieldSSB",range=1/3,length=11)
{
    nspp <- model$param$nspp
    par(mfrow=c(3,4))
    prof_sp <- list()
    for (sp in 1:nspp)
    {
	if(is.null(prof_dat)){
	    prof <- profile1D(model,sp,param,error,range,length)
	    prof_sp[[sp]] <- prof
	}
	else prof <- prof_dat[[sp]]
	plot(x=prof[["param"]],y=prof[["RSS"]],xlab=param,ylab="RSS",type="l")
	lines(x=rep(model$param$species[sp,param],2),y=c(0,1e9),lty=2)
	title(model$param$species$species[sp])
    }
    if(is.null(prof_dat)) return(prof_sp)
}


#test <- profile1D(Finalcomm2a,1,"R0")
#plot(test[["param"]],test[["RSS"]])


calibrateR0_yield.nlslm <- function(logR0hat, Initialcomm, meantsteps=NA)
{
  Initialcomm$param$species$R0 <- exp(logR0hat)

  # Project model with the new R0, continuing from where the initial run left off
  Newcomm <- Setup(Initialcomm$param, ContinueCalculation=T, initialcommunity=Initialcomm)
  Newcomm <- Project(Newcomm)
  return(Yielderror(Newcomm,meantsteps=meantsteps))
}



