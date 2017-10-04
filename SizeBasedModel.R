# Mlla's copy from Dr Julia Blanchard 02/10/2017
# ==============================================

# Models multiple size spectra with different life-history traits.
# Originally based on Pope et al. approach but with food dependent growth and other processes
# following approach of Ken Andersen, Martin Pedersen et al.

# Two functions:
# Setup - takes the model parameter list and creates the objects used in the projection
# Project - projects forward in time
# CheckParams - function to check model object has required parameters to enable elegent failing

# To Do
# Sort out SRR - add it to parameter list
# Multiple iters
# What are we saving?
#         e.g. Why does R and RDI store same thing (ish)
#         Do we need all these objects in the model object? Like N and Nsave?
#         Lots of duplication




# ************************** CHECK PARAMS **************************************

# Function to check parameters are consistent

CheckParams <- function(param)
{
  # Required list elements
  
  # Don't actually need nspp - can calc that from dims of species element - redundant
  required_elements <- c("species", "theta", "nspp", "ngrid", "ngridPP", "tmax", "dt",
                         "isave", "w0", "wMax", "slope0", "wPPcut", "rPP", "kap",
                         "lambda", "Z0pre", "Z0exp", "n", "p", "q")
  
  # Check that params is a list
  if (!is(param,"list"))
    stop(paste("param should be a list with the following elements: ", paste(required_elements, collapse=", ")))  
  
  if (!all(required_elements %in% names(param)))
  {
    missing_elements <- required_elements[!(required_elements %in% names(param))]
    stop(paste("You are missing these elements from the parameter list: ", paste(missing_elements, collapse=", ")))
  }
  
  # Also need to include and check Fishing and SRR parameters - not done yet
  
  # Check species
  # Number of species OK?
  if(!is(param$species,"data.frame"))
    stop("species element should be a data.frame")
  if(nrow(param$species) != param$nspp)
    stop("Number of species in species element does not equal nspp element")
  
  species_cols <- c("Winf", "Wmat", "h", "alpha", "R0", "beta", "sigma", "eRepro", "gamma",
                    "k", "ks", "N0")
  # What about Wegg? Not used yet.
  # Check species matrix has appropriate columns
  if(!all(species_cols %in% names(param$species)))
    stop(paste("The species element must have columns with these names: ", paste(species_cols, collapse= ", ")))
  
  # Check interaction matrix has correct dims
  if(!is(param$theta,"matrix") || all(dim(param$theta) != param$nspp))
    stop("theta (the interaction matrix) must be a square matrix with dim = nspp")
  
  # Check other elements are of required class etc
  #required_elements[-which(required_elements == c("species", "theta"))]
  
  # integers
  integer_elements <- c("nspp", "ngrid", "ngridPP")
  #  if(any(unlist(lapply(param[integer_elements],is,"integer")) == FALSE))
  #    warning(paste(paste(integer_elements, collapse=", "), "should be integers"))
  
  # Other stuff
  if (param$dt >= param$tmax)
    stop("dt must be less than tmax")
  
  # Check if there is a Wegg column in the species data. If not set Wegg to w0
  #if(!("Wegg" %in% names(param$species))) param$species$Wegg <- param$w0
  
  # Check that Wegg is not smaller than w0
  if (any(param$species$Wegg < param$w0)) stop("Egg size is smaller than minimum size class")
  
  
  
  # Having got this far, everything  must be OK
  return (TRUE)
}

# ************************** SETUP *********************************************

# Do we need all these arguments? Continue calc and initial community?
# Yes, but not for project
# Setup returns a 'model' object (a list of everything)

# ngrid and ngridfull are the same - is this right?

Setup <- function(param, ContinueCalculation=FALSE, initialcommunity=NULL,SRR="BH") {
  
  #  print("Set up time")
  #  print(system.time({
  
  #if(!CheckParams(param)) stop ("Problem with the parameters")
  
  # Pull out some useful parameters - just a shortcut
  sp <- param$species
  ngrid <- param$ngrid
  nspp <- param$nspp
  dt <- param$dt
  
  # Set up grid:
  w <- 10^(seq(from=log10(param$w0),to=log10(param$wMax),length.out=ngrid))
  dw <- diff(w)
  dw[ngrid] <- dw[ngrid-1] # Set final dw as same as one before
  
  # Plankton spectrum
  if(is.na(param$ngridPP))
    param$ngridPP <- round(ngrid)*0.3
  ngridPP <- param$ngridPP
  
  # Combine fish and plankton w
  # Might keep seperate for future...
  wFull <- c(10^seq(from=log10(w[1]/(4*max(sp$beta))), to = log10(w[1]-dw[1]),length.out=ngridPP),w)
  ngridfull<-length(wFull) # just a reference, could be calculated on the fly
  dwFull <- diff(wFull)
  dwFull[ngridfull] <- dwFull[ngridfull-1]
  idxGrid <- (ngridPP+1):ngridfull       # short cut to index just the fish
  
  nsave   <- floor(param$tmax/(param$dt*param$isave)) #num of time slots to save
  
  # Make the model object
  # Just a big list full of stuff
  model <- list(
    param = param,
    
    # Grid parameters
    w = w,
    dw = dw,
    wFull = wFull,
    ngridfull = ngridfull,
    dwFull = dwFull,
    idxGrid = idxGrid,
    
    # Species parameters
    psi = matrix(NA,nrow=nspp,ncol=ngrid),
    IntakeMax = matrix(NA,nrow=nspp,ncol=ngrid),
    SearchVol = matrix(NA,nrow=nspp,ncol=ngrid),
    Activity = matrix(NA,nrow=nspp,ncol=ngrid),
    StdMetab = matrix(NA,nrow=nspp,ncol=ngrid),
    predkernel = array(NA,dim=c(nspp,ngrid,length(wFull))),      
    
    
    Z0 = param$Z0pre*sp$Winf^param$Z0exp,    # background natural mortality
    #Z0= (param$Z0pre*sp$k_vb)               # background mortality of Pope et al. 2006
    # Background spectrum
    rrPP = param$rPP*wFull^(param$n-1), #weight specific plankton growth rate ##
    NinfPP = rep(NA,ngridfull),
    
    # Arrays to save output - not everything is saved as some stuff can be calced on the fly after the simulation
    # Can add more if necessary
    # Not everything needs to be saved. Can calculate a lot of these if necessary
    N = array(0,dim=c(nsave,nspp,ngrid)), # abundance
    F = array(0,dim=c(nsave,nspp,ngrid)), # Fishing mortality
    RDI = matrix(0,nrow=nsave,ncol=nspp), # actual recruitment (i.e. physio minus effects of DD)
    RDD = matrix(0,nrow=nsave,ncol=nspp), # physiological recruitment
    M2 = array(0,dim=c(nsave,nspp,length(w))), # predation mortality on fish species
    M2background = matrix(0,nrow=nsave,ncol=length(wFull)),
    #SSB = matrix(0,nrow=nsave,ncol=nspp),
    #SSBm = array(0,dim=c(nsave,nspp,ngrid)),
    eSpawningPopulation = matrix(0,nrow=nsave,ncol=nspp),
    eSpawning = array(0,dim=c(nsave,nspp,ngrid)),
    f = array(0,dim=c(nsave,nspp,ngrid)),  # feeding level
    gg = array(0,dim=c(nsave,nspp,ngrid)), # growth
    Yield = matrix(0,nrow=nsave,ncol=nspp),
    Biomass = matrix(0,nrow=nsave,ncol=nspp),
    nPP = matrix(0,nrow=nsave,ncol=length(wFull))
  )
  
  # Now calculate the species and background parameters
  
  # Set up plankton spectrum
  model$NinfPP <- param$kap*wFull^(-param$lambda) # the resource carrying capacity - one for each mp and m (130 of them)
  model$NinfPP[wFull>param$wPPcut] <- 0      #set density of sizes < plankton cutoff size
  if (ContinueCalculation==TRUE)
    model$nPP[1,] <- initialcommunity$nPP[dim(initialcommunity$nPP)[1],]
  else
    model$nPP[1,] <- model$NinfPP # plankton numbers - set the initial values at carrying cap
  
  #browser()
  
  # Initial population abundance
  if (ContinueCalculation)
    model$N[1,,] = initialcommunity$N[dim(initialcommunity$N)[1],,]
  else
  {
    model$N[1,,] <- unlist(tapply(w,1:length(w),function(wx,N0,w0,slope0,Winf,n)
      N0 * (wx/w0)^slope0,
      N0=sp$N0,w0=param$w0,slope0=param$slope0,Winf=sp$Winf,n=param$n))
    # set densities at w > Winf to 0
    # This next bit is long winded...
    tempN <- model$N[1,,]
    tempN[unlist(tapply(w,1:length(w),function(wx,Winf)Winf<wx,Winf=sp$Winf))] <- 0
    # Also any densities at w < Wmin set to 0
    tempN[unlist(tapply(w,1:length(w),function(wx,Wmin)Wmin>wx,Wmin=sp$Wmin))] <- 0    
    model$N[1,,] <- tempN
  }
  
  # Allocation to reproduction
  model$psi[] <- unlist(tapply(w,1:length(w),function(wx,Winf,Wmat,n){
    ((1 + (wx/(Wmat))^-10)^-1) * (wx/Winf)^(1-n)},Winf=sp$Winf,Wmat=sp$Wmat,n=param$n))
  # Set w < 10% of Wmat to 0
  model$psi[unlist(tapply(w,1:length(w),function(wx,Wmat)wx<(Wmat*0.1)  ,Wmat=sp$Wmat))] <- 0
  # Set all m > M to 1 # Check this is right...
  model$psi[unlist(tapply(w,1:length(w),function(wx,Winf)(wx/Winf)>1,Winf=sp$Winf))] <- 1
  
  model$IntakeMax[] <- unlist(tapply(w,1:ngrid,function(wx,h,n)h * wx^n,h=sp$h,n=param$n))
  model$SearchVol[] <- unlist(tapply(w,1:ngrid,function(wx,gamma,q)gamma * wx^q,gamma=param$species$gamma,q=param$q))
  model$Activity[] <- unlist(tapply(w,1:ngrid,function(wx,k)k * wx,k=param$species$k))
  model$StdMetab[] <- unlist(tapply(w,1:ngrid,function(wx,ks,p)ks * wx^p,ks=param$species$ks,p=param$p))
  
  # Could maybe improve this. Pretty ugly at the moment - but only called once I guess
  model$predkernel[] <- sp$beta
  model$predkernel <- exp(-0.5*sweep(log(sweep(sweep(model$predkernel,3,wFull,"*")^-1,2,w,"*")),1,sp$sigma,"/")^2)
  model$predkernel <- sweep(model$predkernel,c(2,3),combn(wFull,1,function(x,w)x<w,w=w),"*") # find out the untrues and then multiply
  
  #********************* Set up fishing stuff **********************************
  # Fishing stuff
  # Fmort = Sel * Effort * Q
  # Sel is determined by w, function and parameters
  # Selectivity is a 3D array: species, w and gear (for multiple gears)
  
  gears <- unique(model$param$sel_params$gear)
  model$selectivity <- array(0,dim=c(model$param$nspp,length(model$w),length(gears)),
                             dimnames=list(species=sp$species,w=as.integer(w),gear=gears))
  # This is pretty rank code. Might need to rethink structure of selectivity params
  # Idea is to make a list of the parameters, pulled from the selectivity parameters data frame, then
  # call the selectivity function, using do.call, with those parameters
  for (gear_count in gears)
    #    for (sp_count in as.character(sp$species))
    for (sp_count in as.character(unique(model$param$sel_params[model$param$sel_params$gear==gear_count,"species"])))
    {
      selfunc <- as.character(param$sel_params[param$sel_params$gear==gear_count & param$sel_params$species==sp_count,"selfunc"][1])
      # Pretty horrible
      # Get parameters for selectivity function
      params <- param$sel_params[param$sel_params$gear==gear_count & param$sel_params$species==sp_count,c("param_name","param_value")]
      # Just checking that we only include the params we need for the function
      params <- params[params$param_name %in% names(formals(as.character(selfunc))),]
      param_list <- as.list(params$param_value)
      names(param_list) <- params$param_name
      param_list[["w"]] <- w
      #param_list <- list(unlist(param_list),w = seq(1,5000,length=20))
      model$selectivity[sp_count,,gear_count] <- do.call(selfunc,param_list)
    }
  
  #browser()
  
  # Effort is 2D array - time by gear - set to 1 for all time to start with
  itimemax  <- param$tmax / dt  #max index of time array
  model$effort <- array(1,dim=c(itimemax,length(gears)),dimnames=list(time=1:itimemax, gear=gears))
  
  return(model)
  # return everything...
}
#*******************************************************************************

# ngridfull is in model, not in param, is this right?


# Given a model with some inital values in N and nPP, project forward

Project <- function(model,sd_rrPP=0, sd_srr=0, backgroundAutocorrelationFactor = 0, fixedFeedingLevel = NA, fixedPredationMortality = NA)
{
  #browser()
  # Pull out some useful parameters - just a shortcut
  param <- model$param
  sp <- param$species
  ngrid <- param$ngrid  # ngridfull?
  nspp <- param$nspp
  dt <- param$dt
  w <- model$w
  wFull <- model$wFull
  dw <- model$dw
  dwFull <- model$dwFull
  alpha <- sp$alpha
  eRepro <- sp$eRepro
  theta <- param$theta
  
  # Handy stuff
  idx <- 2:ngrid
  itimemax  <- param$tmax / dt  #max index of time array
  
  # Matrices for solver
  # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
  A <- matrix(0,nrow=nspp,ncol=ngrid)
  B <- matrix(0,nrow=nspp,ncol=ngrid)
  S <- matrix(0,nrow=nspp,ncol=ngrid)
  
  # Temporary Matrices that get updated each time step
  # some of these saved for output
  N <- model$N[1,,]
  nPP <- model$nPP[1,]
  eSpawning <- matrix(0,nrow=nspp,ncol=ngrid)
  eSpawningPopulation <-  matrix(0,nrow=1,ncol=nspp)
  RDI <- matrix(0,nrow=1,ncol=nspp) # density independent recruitment  
  RDD <- matrix(0,nrow=1,ncol=nspp) # density dedependent recruitment    
  
  # Set up helpful matrices to avoid using sweep statements
  wmat <- matrix(w,byrow=T,nrow=nspp,ncol=ngrid)
  dwmat <- matrix(dw,byrow=T,nrow=nspp,ncol=ngrid)
  wFullmat <- matrix(wFull,byrow=T,nrow=nspp,ncol=length(wFull))
  dwFullmat <- matrix(dwFull,byrow=T,nrow=nspp,ncol=length(wFull))
  
  # Generate stochastic rrPP bits using Ranta and Kaitala's routine
  #backgroundAutocorrelationFactor <- 0.5 #  = 0 = white noise, -ve blue
  Znorm <- rnorm(itimemax,0,sd_rrPP)
  rrPP_noise <- numeric(itimemax)
  rrPP_noise[1] <- Znorm[1]
  for (i in 2:itimemax) rrPP_noise[i] <- backgroundAutocorrelationFactor*rrPP_noise[i-1] + Znorm[i]
  rrPP_noise <- exp(rrPP_noise)
  # par(mfrow=c(2,1))
  # plot(rrPP_noise, type="l")
  # compared to just standard rlnorm
  # plot(rlnorm(itimemax,0,sd_rrPP), type="l")
  
  #browser()
  
  # Which size element should eggs be dumped into
  w0idx <- as.vector(tapply(sp$Wmin,1:length(sp$Wmin),function(Wmin,wx) max(which(wx<=Wmin)),wx=w))
  # hacky shortcut so I can access the correct element of a 2D array using 1D notation
  w0idx_arrayref <- (w0idx-1) * nspp + (1:nspp)
  
  # BIG TIME LOOP
  for (itime in 1:itimemax)
  {
    # browser()
    # Calculate feeding level of all predators
    # This is the numbers of prey at a particular size that are exposed to each predator
    if (length(dim(theta))==2)
      Neff <- theta %*% N # use interaction matrix
    # Neff <- theta[,,1,1] %*% N
    # If theta is size disaggregated
    if (length(dim(theta))==4)
      #Neff <- apply(sweep(theta,c(2,4),N,"*"),c(1,3,4),sum)
      Neff <- colSums(aperm(sweep(theta,c(2,4),N,"*"),c(2,1,3,4))) # Apply is too slow
    
    # Neff is pred species x pred size x prey size
    
    # Neff: predator species x prey size
    # total number of prey at that size that are exposed to that predator
    
    # Update phiprey - not sure we can speed this up anymore
    # Fish eat from: background plankton spectrum + fish spectrum
    # predkernel: species x predator size x prey size
    # original, with 2D Neff
    if (length(dim(theta))==2){
      phiprey <- rowSums(sweep(model$predkernel,3,dwFull*wFull*nPP,"*"),dims=2) +
        rowSums(sweep(model$predkernel[,,model$idxGrid],c(1,3),Neff*dwmat*wmat,"*"),dims=2)
    }
    
    if (length(dim(theta))==4)
    {
      phiprey <- rowSums(sweep(model$predkernel,3,dwFull*wFull*nPP,"*"),dims=2) +
        rowSums(sweep(model$predkernel[,,model$idxGrid]*Neff,c(1,3),dwmat*wmat,"*"),dims=2)
    }
    
    # Encountered food
    encount <- model$SearchVol*phiprey
    ## Calculate feeding level:
    f <- encount/(encount + model$IntakeMax)
    # Or forget all the above and just fix feeding level!
    if (!is.na(fixedFeedingLevel))
      f[] <- fixedFeedingLevel
    
    # Predation mortality
    # pred kernel is how much each species (1) by mass (2), predates on other masses (3)
    # 12 x 100 x 130
    predrate <- sweep(model$predkernel,c(1,2),(1-f)*model$SearchVol*N*dwmat,"*")
    
    #browser()
    
    # Natural mortality
    if (length(dim(theta))==2)
      M2 <- theta %*% colSums(aperm(predrate, c(2,1,3)),dims=1)[,model$idxGrid]
    # If theta is size disaggretated
    if (length(dim(theta))==4)
      #M2 <- apply(sweep(theta, c(1,3,4), predrate[,,model$idxGrid], "*"),c(2,4),sum)
      M2 <- (rowSums(aperm(sweep(theta, c(1,3,4), predrate[,,model$idxGrid], "*"), c(2,4,1,3)),dims=2)) # Aieee
    # Or use the C version much faster
    #M2 <- M24D(theta,predrate[,,model$idxGrid])
    # Or forget the above and just fix the predation mortality for the fish (not the background)
    if(!is.na(fixedPredationMortality))
      M2[] <- fixedPredationMortality
    
    
    M2background <- 1 *colSums(predrate,dims=2)
    
    # Calc. assimilated intake:
    e <- sweep(f * model$IntakeMax,1,alpha,"*")
    
    # Subtract basal metabolism and activity (Activity has been set to 0 though)
    e <- e - model$StdMetab - model$Activity
    
    e[e<0] <- 0 ## Do not allow negative growth
    
    # Calculate the amount of energy allocated to gonads - same as gir in Appendix. e = alpha * f * Cmax - ks * m^p
    # energy left for reproduction * allocation to reproduction
    eSpawning <- model$psi * e     # as in matlab code Spectrum.m
    # subtract from assimilated food and use the rest for growth
    gg <- e - eSpawning
    
    # Update Fishing mortality based on effort at time itime
    # F by gear = sel * Q * effort
    
    Fmort <- sweep(model$selectivity,c(1,3),sweep(model$param$Q,2,model$effort[itime,],"*"),"*")
    # Total mortality - could be faster to make Z0 a matrix to make addition simple
    # Need to sum Fmort over gears
    Z = sweep(M2 + rowSums(Fmort,dims=2),1,model$Z0,"+")
    
    # Iterate species one time step forward:
    # Set up matrix:
    A[,idx] <- -gg[,idx-1]*dt/dwmat[,idx]
    B[,idx] <- 1 + gg[,idx]*dt/dwmat[,idx] + Z[,idx]*dt
    S[,idx] <- N[,idx]
    
    # Boundary condition upstream end (recruitment)
    eSpawningPopulation[] <- (eSpawning*N) %*% dw
    #**** Change w[1] here to w[eggsize] **********
    #RDI[] <- 0.5*(eSpawningPopulation*eRepro)/w[1] # Density independent recruitment, 0.5 from sex ratio
    RDI[] <- 0.5*(eSpawningPopulation*eRepro)/w[w0idx] # Density independent recruitment, 0.5 from sex ratio
    #**** Change gg[1] here to gg[eggsize], same for B and Z  **********
    #B[,1] <- 1+gg[,1]*dt/dw[1]+Z[,1]*dt
    B[w0idx_arrayref] <- 1+gg[w0idx_arrayref]*dt/dw[w0idx]+Z[w0idx_arrayref]*dt
    # Call the stock recruitment function in the parameter list
    # to calculate the density dependent recruitment
    # Where to add noise to SRR?
    # Add noise to R0
    #tempParam <- model$param
    #tempParam$species$R0 <- tempParam$species$R0 * rlnorm(nspp,0,sd_srr)
    #RDD[] <- model$param$SRR(eSpawningPopulation, RDI, tempParam)
    # Or add noise directly to realised recruitment 
    RDD[] <- model$param$SRR(eSpawningPopulation, RDI, model$param)
    # RDD <- RDD*rlnorm(nspp,0,sd_srr)
    # sd is sd on the log scale. So the following line gives noise with median of RDD and sd of log RDD = sd_srr
    # RDD[] <- rlnorm(nspp,log(RDD),sd_srr)
    # Then drop in the recruitment    
    #**** Change ref from 1 to eggsizref ********
    #N[,1] <- (N[,1] + RDD*dt/dw[1]) / B[,1]
    N[w0idx_arrayref] <- (N[w0idx_arrayref] + RDD*dt/dw[w0idx]) / B[w0idx_arrayref]
    # Add noise directly to realised density of eggs (as Andersen & Pedersen 2010, with sd_srr standardised from data) 
    N[w0idx_arrayref] <-  N[w0idx_arrayref]*rlnorm(nspp,0,sd_srr)  
    
    # Invert matrix
    # idx used to start at 2 because it was 1 after when the recruits arrived
    # It has to start 1 after Wmin, else numbers at Wmin get overwritten with 0
    # But, now that can differ between species
    # So we have the dreaded nested for loop
    
    # Using inlined C function for speed - see ../C/invert_matrix.r
    #N <- invert_matrixC(nspp,ngrid,w0idx,N,S,A,B);
    # Or use R code
    for (i in 1:nspp)
      for (j in (w0idx[i]+1):ngrid)
        N[i,j] <- (S[i,j] - A[i,j]*N[i,j-1]) / B[i,j]
    
    #browser()
    # Calc. background spectrum:
    # rrPP is the weight specific plankton growth rate
    # If stoch on growth rate
    #rrPPstoch <- model$rrPP * rrPP_noise[itime]
    #tmp <- (rrPPstoch*model$NinfPP / (rrPPstoch + M2background))
    #nPP[] <- tmp - (tmp - nPP) * exp( -(rrPPstoch+M2background)*dt)
    
    # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
    tmp <- (model$rrPP*model$NinfPP / (model$rrPP + M2background))
    nPP[] <- tmp - (tmp - nPP) * exp( -(model$rrPP+M2background)*dt)
    nPP <- nPP * rrPP_noise[itime]
    #
    
    # Save results:
    if ((itime %% param$isave)==0)
    {
      isav<-itime/param$isave
      model$f[isav,,]<- f
      model$gg[isav,,]<- gg
      model$N[isav,,]<-N
      model$RDD[isav,] <- RDD  # the recruitment after density dependence
      model$RDI[isav,] <- RDI # the density independent Recruitment
      model$F[isav,,]<-rowSums(Fmort,dims=2) # sum Fmort over gears
      model$M2[isav,,] <- M2
      model$M2background[isav,] <- M2background
      model$eSpawningPopulation[isav,] <- eSpawningPopulation
      model$eSpawning[isav,,] <- eSpawning
      model$nPP[isav,] <- nPP
    }
  }  #end of for time loop
  
  # Calculate these for easy reference  
  model$Yield[] <- rowSums(sweep((model$F * model$N),3,w*dw,"*"),dims=2)
  model$Biomass[] <- rowSums(sweep(model$N,3,w*dw,"*"),dims=2)
  
  #browser()
  
  #Fmort <- sweep(model$selectivity,c(1,3),sweep(model$param$Q,2,model$effort[itime,],"*"),"*")
  
  # Yield by time x species x size x gear
  #model$YieldByGear <- 
  
  return(model)
}

# Ask Julia about fleet structure

# Economic indicators
# Price per species by size - Jim
# Single gear - project forward under different constant effort to see yield / revenue against effort
# By species and total

# Net present value - discount rates
# 

# 
# Different gears will require costs per 
# Price per effort

# HCR - just calc total revenue on +-20% of effort
# (calc. revenue on last time step effort, +-20% of that effort - pick the highest projected revenue
# based on some cost assumptions)

# Sensitivity analysis - change in revenue / profit to changes in prices / costs
