# Mlla's copy from Dr Julia Blanchard 02/10/2017
# Plots
#library(odesolve)
library(deSolve)
library(reshape)
library(ggplot2)
library(grid)

# Feeding level against w
# M2 against w
# Biomass spectrum (with background community too) against w
# Biomass against time
# Relative weight (?) Not sure what that one is

# What was the growth curve one from Ken?
# Length against Age
# Check his matlab code...

plotFeedinglevel <- function(model, meantsteps = NA, plotleg = T)
{
  # get y data, either final point
  if(is.na(meantsteps))
    y <- model$f[dim(model$f)[1],,]
  else # or mean of last meantsteps
    y <- apply(model$f[(dim(model$f)[1]-meantsteps+1):dim(model$f)[1],,],c(2,3),mean)
  
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  #  lty <- rep(c(1,2), each = dim(model$f)[2]/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  # Get y lims
  ylim <- c(0,1)
  # Plot empty
  plot(x=model$w, y=y[1,], log="x",ylim=ylim, type="n", xlab = "mass (g)", ylab = "feeding level")
  # Plot the lines
  # Only plot w <= winf
  for (i in 1:model$param$nspp)
    points(x=model$w[model$w <= model$param$species$Winf[i]], y=y[i,model$w <= model$param$species$Winf[i]],
           col=col[i], type="l", lty=lty[i])
  
  if (plotleg == T)
    legend(x = "bottomright" , legend = as.character(model$param$species$species),
           col= col, lty=lty, cex=0.7, ncol=2)
  
}

plotM2 <- function(model, meantsteps = NA, plotleg=T)
{
  # Need to fix w range of M2
  #(length(NSmodel$wFull)-length(NSmodel$w)+1):length(NSmodel$wFull)
  # get y data, either final point
  if(is.na(meantsteps))
    #    y <- model$M2[dim(model$M2)[1],,(length(model$wFull)-length(model$w)+1):length(model$wFull)]
    y <- model$M2[dim(model$M2)[1],,]
  else # or mean of last meantsteps
    #    y <- apply(model$M2[(dim(model$M2)[1]-meantsteps+1):dim(model$M2)[1],,(length(model$wFull)-length(model$w)+1):length(model$wFull)],c(2,3),mean)
    y <- apply(model$M2[(dim(model$M2)[1]-meantsteps+1):dim(model$M2)[1],,],c(2,3),mean)
  
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  #  lty <- rep(c(1,2), each = dim(model$f)[2]/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  # Get y lims
  ylim <- c(0,max(y))
  
  # Plot empty
  plot(x=model$w, y=y[1,], log="x",ylim=ylim, type="n", xlab = "mass (g)", ylab = "natural mortality")
  # Plot the lines
  for (i in 1:model$param$nspp)
    points(x=model$w[model$w <= model$param$species$Winf[i]], y=y[i,model$w <= model$param$species$Winf[i]],
           col=col[i], lty=lty[i], type="l")
  
  if (plotleg == T)
    legend(x = "topright" , legend = as.character(model$param$species$species), col= col, lty=lty, cex=0.7, ncol=2)
}

# This needs to be corrected?
# Biomass is sum(N*w*dw), not just N * dw
# Should just be N * w because abundance is point, not abundance in that mass bucket
plotBioSpec <- function(model,meantsteps = NA, plotleg = T, main=NULL)
{
  #browser()
  # Calculate Biomass by w
  
  # w and dw are constant
  # What N are we using
  if(is.na(meantsteps))
    N <- model$N[dim(model$N)[1],,]
  else # or mean of last meantsteps
    N <- apply(model$N[(dim(model$N)[1]-meantsteps+1):dim(model$N)[1],,],c(2,3),mean)
  
  # w is in grams
  #spBiomass <- sweep(N,2,model$dw,"*") / 1e6
  # Want to calculate biomass at each point (not in each size class)
  # N is abundance at point, not abundance in that size class
  spBiomass <- sweep(N,2,model$w,"*") / 1e6
  PPBiomass <- model$nPP[dim(model$nPP)[1],] * model$wFull / 1e6
  refSpec <- model$param$kap*model$wFull^(1-model$param$lambda) / 1e6
  
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  #lty <- rep(c(1,2), each = dim(model$f)[2]/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  # Set xlim and ylim values
  xlim <- c(1e-2, max(model$wFull)) # don't plot the full resource spectrum
  ylim <- c(1e-6,max(PPBiomass,spBiomass,refSpec))
  
  # Set it up so x is full spectrum
  plot(x=model$wFull, y=PPBiomass, log="xy",type="n", ylim=ylim, xlim=xlim, xlab = "mass (g)", ylab = "biomass (t)", main=main)
  
  # Reference spectrum:
  points(x=model$wFull, y= refSpec, type="l", lty=3, col = 1)
  # Resource spectrum:
  points(x=model$wFull, y=PPBiomass, type="l", lty=3, col = 3)
  
  # Species spectrum
  for (i in 1:dim(spBiomass)[1])
    points(x=model$w, y = spBiomass[i,], col=col[i], type="l", lty=lty[i])
  
  if (plotleg==T) 
    legend(x="topright", legend = c("Reference", "Background", as.character(model$param$species$species)),
           lty = c(3,3,lty), col=c(1,3,col), cex=0.7, ncol=2) 
}

plotNSpec <- function(model,meantsteps = NA, plotleg = T, main=NULL, addCommSlope=FALSE, minw=NA)
{
  
  # w and dw are constant
  # What N are we using
  if(is.na(meantsteps))
    N <- model$N[dim(model$N)[1],,]
  else # or mean of last meantsteps
    N <- apply(model$N[(dim(model$N)[1]-meantsteps+1):dim(model$N)[1],,],c(2,3),mean)
  
  # w is in grams
  #spBiomass <- sweep(N,2,model$dw,"*") / 1e6
  # Want to calculate biomass at each point (not in each size class)
  # N is abundance at point, not abundance in that size class
  spN <- N / 1e6
  PPN <- model$nPP[dim(model$nPP)[1],] / 1e6
  refSpec <- model$param$kap*model$wFull^(-model$param$lambda) / 1e6
  
  # Fit and add the community slope if you want to
  if (addCommSlope){
    if (is.na(minw)) minw <- min(model$w)
    winc <- model$w >= minw
    # Not taking mean of last tsteps, just last timestep
    nt <- dim(model$N)[1]
    if (is.na(meantsteps))
      Nt <- model$N[nt,,]
    else
      Nt <- apply(model$N[(nt-meantsteps+1):nt,,],c(2,3),mean)
    
    Ntsum <- apply(Nt,2,sum)
    fit <- lm(log(Ntsum[winc]) ~ log(model$w[winc]))
    print(summary(fit))
  }
  
  
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  #lty <- rep(c(1,2), each = dim(model$f)[2]/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  # Set xlim and ylim values
  xlim <- c(1e-2, max(model$wFull)) # don't plot the full resource spectrum
  ylim <- c(1e-6,max(PPN,spN,refSpec))
  
  # Set it up so x is full spectrum
  plot(x=model$wFull, y=PPN, log="xy",type="n", ylim=ylim, xlim=xlim, xlab = "mass (g)", ylab = "Abundance", main=main)
  
  # Reference spectrum:
  points(x=model$wFull, y= refSpec, type="l", lty=3, col = 1)
  # Resource spectrum:
  points(x=model$wFull, y=PPN, type="l", lty=3, col = 3)
  
  # Species spectrum
  for (i in 1:dim(spN)[1])
    points(x=model$w, y = spN[i,], col=col[i], type="l", lty=lty[i])
  
  #browser()
  if (addCommSlope){
    inter <- log(mean(Nt[, min(which(winc))])) - fit$coefficients[2] * log(model$w[min(which(winc))]) # sort this out
    commSpec <-  exp(inter + fit$coefficients[2] * log(model$w[winc]))
    #commSpec <-  exp(fit$coefficients[2] * log(model$w[winc]))
    lines(x=model$w[winc], y =commSpec / 1e6, lwd = 2)
  }
  
  if (plotleg==T)
    legend(x="topright", legend = c("Reference", "Background", as.character(model$param$species$species)),
           lty = c(3,3,lty), col=c(1,3,col), cex=0.7, ncol=2)
}


# Plot Total Biomass through time
plotBioTime <- function(model, trange=NA, plotleg=T, main=NULL)
{
  #browser()
  # Someway of specifying trange
  if (is.na(trange))
    trange <- 1:dim(model$N)[1]
  # Check trange is OK
  if (!all(trange %in% (1:dim(model$N)[1])))
    stop("User specified trange is outside time range of model")
  
  # N * w then summed across w to get total biomass through time
  #spBiomass <- apply(sweep(model$N,c(3,2),model$w / 1e6,"*"),c(1,2),sum)
  spBiomass <- rowSums(sweep(model$N,3,model$w*model$dw,"*"),dims=2)
  
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  
  ylim <- c(1,max(spBiomass))
  plot(x=trange,y=trange,type="n",ylim=ylim,log="y",xlab = "timestep", ylab = "total biomass (t)", main=main)  
  
  for (i in 1:dim(spBiomass)[2])
    points(x=trange, y = spBiomass[trange,i], col=col[i], lty=lty[i], type="l")
  
  if (plotleg == T)
    legend(x = "bottomright" , legend = as.character(model$param$species$species), col= col, lty=lty, cex=0.7, ncol=2)
  return(invisible(spBiomass))
}


# Some kind of SRR plot?
# R0 is the max R
# R is the DI R
# Rtemp is the actual R
#Rtemp <-  sp$R0 * model$R / (sp$R0+model$R)


#R0 <- 3.3
#R <- seq(from=0, to = 100, by = 0.1)
#Rt <- (R0 * R) / (R0 + R)
#plot(R,Rt)

plotSRRratios <- function(model, meantsteps=NA, plotleg=T)
{
  # get y data, either final point
  if(is.na(meantsteps))
  {
    RDI <- model$RDI[dim(model$RDI)[1],]
    RDD <- model$RDD[dim(model$RDD)[1],]    
  }
  else # or mean of last meantsteps
  {
    RDI <- apply(model$RDI[(dim(model$RDI)[1]-meantsteps+1):dim(model$RDI)[1],],2,mean)
    RDD <- apply(model$RDD[(dim(model$RDD)[1]-meantsteps+1):dim(model$RDD)[1],],2,mean)
  }
  
  Rmax <- model$param$species$R0
  DDRmaxrat <- RDD / Rmax # What proportion of the maximum recruitment was achieved
  DDDIrat <- RDD / RDI # By how much was the potential recruitment achieved
  
  ylim <- c(0,1.2)
  
  plot(x=1:length(RDI), y=1:length(RDI), type="n", ylim=ylim, xlab="species", ylab="ratio", xaxt="n")
  axis(1,at=1:length(RDI),labels=as.character(model$param$species$species),cex.axis=0.7)
  # Line at 1
  points(x=c(1,length(RDI)),y=c(1,1), type="l", lty=2)
  # ratio DD recruitment with max recruitment  
  points(x=1:length(RDI), y = DDRmaxrat, pch=0)
  # ratio DD rec with DI rec
  points(x=1:length(RDI), y = DDDIrat, pch=16)
  
  if(plotleg==T)
    legend(x="right",legend=c("DD rec. / max rec.", "DD rec. / DI rec."), pch=c(0,16), cex=0.7)
  
}



plotResults <- function(model, trange=NA, meantsteps=NA)
{
  #browser()
  par(mfrow=c(3,2))
  plotFeedinglevel(model,meantsteps, plotleg=F)
  plotM2(model,meantsteps, plotleg=F)
  #plotBioSpec(model,meantsteps, plotleg=F)
  plotNSpec(model,meantsteps, plotleg=F)
  plotBioTime(model,trange, plotleg=F)
  #  plotSRRratios(model, meantsteps, plotleg=T)
  plotF(model,plotleg=F)
  
  # Add legend as final plot
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  #lty <- rep(c(1,2), each = dim(model$f)[2]/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  plot(x=1:40,y=1:40,type="n",axes=F,xlab="",ylab="")
  legend(x="left", legend = c("Reference", "Background", as.character(model$param$species$species)),
         lty = c(3,3,lty), col=c(1,3,col), cex=1, ncol = 2, bty="n") 
  
  
}


# Selectivity Plots by gear and species
plotSelectivity <- function(model, plotleg=T)
{
  col <- rep(rainbow(ceiling(model$param$nspp/2),start=0,end=5/6),2)
  #lty <- rep(c(1,2), each = model$param$nspp/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  # Need to adapt for multiple gears
  gear <- 1
  # Set up axes
  plot(model$w,model$selectivity[1,,gear],log="x",type="n",ylim=c(0,1),xlab="log mass (g)",ylab="Proportion")
  #plot(model$w,model$selectivity[1,,gear],type="n",ylim=c(0,1),xlab="mass (g)",ylab="Proportion")
  #plot(x = log10(model$w),model$selectivity[1,,gear],type="n",ylim=c(0,1),xlab="log10 mass (g)",ylab="Proportion")
  for (i in 1:model$param$nspp)
    points(x=model$w,y=model$selectivity[i,,gear],col=col[i], lty=lty[i], type="l")
  #points(x=log10(model$w[model$w<=model$param$species$Winf[i]]),y=model$selectivity[i,model$w<=model$param$species$Winf[i],gear],col=col[i], lty=lty[i], type="l")
  
  if (plotleg == T)
    legend(x = "bottomright" , legend = as.character(model$param$species$species), col= col, lty=lty, cex=0.7, ncol=2)
}

# Total Fishing mortality (summed over gears)
plotF <- function(model,plotleg=T)
{
  col <- rep(rainbow(ceiling(model$param$nspp/2),start=0,end=5/6),2)
  #lty <- rep(c(1,2), each = model$param$nspp/2)
  lty <- rep(c(1,2), each = ceiling(dim(model$f)[2]/2))
  ymax <- max(model$F)
  
  # Set up axes
  plot(x=model$w,y=model$F[dim(model$F)[1],1,],log="x",type="n",ylim=c(0,ymax),xlab="mass (g)",ylab="Total F")
  for (i in 1:model$param$nspp)
    points(x=model$w[model$w<=model$param$species$Winf[i]],y=model$F[dim(model$F)[1],i,model$w<=model$param$species$Winf[i]],col=col[i], lty=lty[i], type="l")
  
  if (plotleg == T)
    legend(x = "bottomright" , legend = as.character(model$param$species$species), col= col, lty=lty, cex=0.7, ncol=2)
}

# Implements the somatic growth equation through time
# Assumes equilibrium to get the feeding level
# So growth is determined by the simulated feeding level at equib
# Question is: does this match the growth according to VB?
# dw/dt = g(m) = (food intake - metabolism) (allocation to growth)
# dw/dt = g(m) = (1 - psi) ( a f(m) Cmax - ks m^p)
# where:
#   Cmax = h m^n
#   psi = allocation to reproduction
somaticGrowth <- function(t,ww,parms)
{
  #browser()
  model = parms[["model"]]
  sp = parms[["sp"]]
  # must call correct t - iPlot? f at t? No
  ix <- which(model$w >= ww)[1]
  if (is.na(ix))
    return(dw=list(0))
  # final feeding level, assumed to be at equib, pick last value
  #model$f[dim(model$f)[3],sp,ix]
  return(list(dw=(model$param$species$alpha[sp] * model$f[dim(model$f)[1],sp,ix] * model$param$species$h[sp] *
                    ww^model$param$n - model$param$species$ks[sp] * ww^model$param$p) *
                (1 - model$psi[sp,ix])))
}

# Not sure where this comes from
# Some sort of approximation using maintenance = growth (and asymptotic size)
Tmax <- function(model)
{
  return(60 * model$param$species[,"Wmat"] ^ (1-model$param$n) /
           (model$param$species[,"alpha"] * model$param$f0est * (model$param$species[,"h"] -model$param$species[,"ks"])))
}

calcGrowth <- function(model,sp)
{
  #  growth <- list()
  #  for (sp in 1:model$param$nspp)
  #  {
  t <- seq(0,Tmax(model)[sp],length=20)
  ww <- model$w[1]
  tw <- lsoda(ww,t,somaticGrowth,parms=list(model=model,sp=sp))
  #    growth[[sp]] <- tw
  #  }
  #  return(growth)
  return(tw)
}

# lsoda
#sp <- 3
#t <- seq(0,Tmax(NSmodel)[sp],length=100)
##t <- seq(0,40,length=100)
#ww <- NSmodel$w[1]
#test <- lsoda(ww,t,somaticGrowth,parms=list(model=NSmodel,sp=sp))
#plot(x=test[,1],y=test[,2])
#

# VB:
# Lt = Linf(1 - exp(-K(t-t0)))
#


plotGrowthAll <- function(model,plotleg=T)
{
  #browser()
  # Make this more general - for any number of species
  #par(mfrow=c(4,5))
  par(mfrow=c(3,4))
  species <- model$param$species
  # Get t0 from VB equation
  if(is.null(species$Linf))
    species$Linf <- (species$Winf / species$a)^(1/species$b)
  # why is this -ve at front (orig Matlab code)?
  # Was to do with the way the VB growth was phrased. t0 is time that length = 0
  #  t0 <- -log(1-(model$param$w0/species$a)^(1/species$b) / species$Linf) / species$k_vb
  
  t0 <- log(1-(model$param$w0/species$a)^(1/species$b) / species$Linf) / species$k_vb
  for (sp in 1:model$param$nspp)
  {
    tw <- calcGrowth(model,sp)
    weightVB = species$a[sp] * (species$Linf[sp] * (1 - exp(-species$k_vb[sp] * (tw[,1]-t0[sp])))) ^ species$b[sp]
    #plot(x=tw[,1],y=tw[,2],type="l",xlab="time", ylab="mass", ylim=c(0,species$Winf[sp]))
    plot(x=tw[,1],y=tw[,2],type="n",xlab="time", ylab="mass", ylim=c(0,max(c(tw[,2], weightVB))))
    lines(x=tw[,1],tw[,2],lty=1,col=2)
    lines(x=tw[,1],y=weightVB,lty=2,col=2)
    # calculate mean error relative to Winf
    err <- sqrt(mean(((tw[,2]-weightVB)/species$Winf[sp])^2))
    title(paste(species$species[sp],signif(err,3),sep=" "))
  }
}

#**************************************************************************
plotWinf <- function(model,plotleg=T)
{
  col <- rep(rainbow(ceiling(dim(model$f)[2]/2),start=0,end=5/6),2)
  pch <- rep(c(16,17), each = dim(model$f)[2]/2)
  par(mfrow=c(2,1))
  #  par(mfrow=c(1,3))
  # R0 against Winf
  plot(x=model$param$species$Winf,y=model$param$species$R0*1e-6,
       log="xy", type="n", xlab="Winf (g)", ylab="Rmax (t)")
  for (i in 1:model$param$nspp)
    points(x=model$param$species$Winf[i],y=model$param$species$R0[i]*1e-6,col=i,pch=pch[i])
  if (plotleg == T)
    legend(x = "bottomleft" , legend = as.character(model$param$species$species), col= col, pch=pch, cex=0.7, ncol=2)
  
  # eRepro against Winf
  plot(x=model$param$species$Winf,y=model$param$species$eRepro,
       log="xy", type="n", xlab="Winf (g)", ylab="eRepro")
  for (i in 1:model$param$nspp)
    points(x=model$param$species$Winf[i],y=model$param$species$eRepro[i],col=i,pch=pch[i])
  
  #    plot(1, 1, xlim=c(0,20),ylim=c(0,1.2),type="n", lwd=3, bty="n",axes="F",xlab="",ylab="")
  #  if (plotleg == T)
  #    legend(x = "bottomleft" , legend = as.character(model$param$species$species), col= col, pch=pch, cex=0.7, ncol=2)
  
  
}

#*******************************************************************************
# Plots the predation rates by predator on prey mass
plotPredationRates <- function(model, include_background=TRUE)
{
  # Set up some short cuts
  f <- model$f[dim(model$f)[1],,]
  sv <- model$SearchVol
  N <- model$N[dim(model$N)[1],,]
  dwmat <- matrix(rep(model$dw,dim(N)[1]),nrow=dim(N)[1],byrow=TRUE)
  pk <- model$predkernel
  m2 <- model$M2[dim(model$M2)[1],,]
  
  # predation by each species, by predator mass on prey mass
  predratefull <- sweep(pk,c(1,2),(1-f)*sv*N*dwmat,"*")
  
  if (!include_background)
  {
    # Collapse over the predator mass to get total predation by each predator on each prey mass
    predrate <- apply(predratefull,c(1,3),sum)[,model$idxGrid] # slower than the colsums / aperm method
    dimnames(predrate) <- list(Species = model$param$species$species,
                               Mass = model$w)
  }
  if (include_background)
  {
    # All mass sizes including background
    predrate <- apply(predratefull,c(1,3),sum) # slower than the colsums / aperm method
    dimnames(predrate) <- list(Species = model$param$species$species,
                               Mass = model$wFull)
  }
  
  # Turn into a dataframe and remove 0s
  preddf <- melt(predrate)
  preddf[preddf$value==0,"value"] <- NA
  # Add a column of log10 Mass to get the widths
  preddf <- cbind(preddf,log10Mass = log10(preddf$Mass))
  # Add column of dws to act as widths for tile plot - needed for including background
  preddf <- ddply(preddf, .(Species), transform, dlog10Mass = c(diff(log10Mass)[1],diff(log10Mass)))
  
  # And plot
  d <- ggplot(preddf) + geom_tile(aes(x=log10Mass,y=Species,fill=(value),width=dlog10Mass))
  # log the values if you want
  #d <- d + scale_fill_gradient2(name="Predation rate", trans="log")
  d <- d + scale_fill_gradient2(name="Predation rate")
  d <- d + scale_y_discrete(name="Predator") +
    scale_x_continuous(name=expression(paste(Log[10]," Prey Mass (g)", sep="")))
  
  return(d)
}

# Plots the predation rates on each prey by predator
plotPredationRatesbyPrey <- function(model,sp=NULL)
{
  #browser()
  
  # Set up some short cuts
  f <- model$f[dim(model$f)[1],,]
  sv <- model$SearchVol
  N <- model$N[dim(model$N)[1],,]
  dwmat <- matrix(rep(model$dw,dim(N)[1]),nrow=dim(N)[1],byrow=TRUE)
  pk <- model$predkernel
  m2 <- model$M2[dim(model$M2)[1],,]
  theta <- model$param$theta
  
  # Names in theta may not be same as in species file so rename
  dimnames(theta)[[1]] <- model$param$species$species
  dimnames(theta)[[2]] <- model$param$species$species
  
  # predation by each species, by predator mass on prey mass
  predratefull <- sweep(pk,c(1,2),(1-f)*sv*N*dwmat,"*")
  # Collapse over the predator mass to get total predation by each predator on each prey mass
  predrate <- apply(predratefull,c(1,3),sum)[,model$idxGrid] # slower than the colsums / aperm method
  dimnames(predrate) <- list(Species = model$param$species$species,
                             Mass = model$w)
  
  # Turn into a dataframe and remove 0s
  preddf <- melt(predrate)
  preddf[preddf$value==0,"value"] <- NA
  # Add a column of log10 Mass to get the widths
  preddf <- cbind(preddf,log10Mass = log10(preddf$Mass))
  # And change the name of Species to Predator
  names(preddf)[names(preddf)=="Species"] <- "Predator"
  
  # name theta
  dimnames(theta) = list(Predator=dimnames(theta)[[1]], Prey=dimnames(theta)[[2]])
  thetadf <- melt(theta)
  names(thetadf)[names(thetadf)=="value"] <- "Interaction"
  # Combine the dataframes so we get interaction and predation rates
  preypred <- join(preddf,thetadf,by="Predator")
  preypred <- cbind(preypred, scaledpred = preypred$Interaction * preypred$value)
  
  #browser()
  
  # Awesome!
  if (is.null(sp))
    d <- ggplot(preypred) + geom_tile(aes(x=log10Mass,y=Predator,fill=(scaledpred))) + facet_wrap(~Prey)
  if (!is.null(sp))
    d <- ggplot(preypred[preypred$Prey %in% sp,]) + geom_tile(aes(x=log10Mass,y=Predator,fill=(scaledpred))) + facet_wrap(~Prey)
  # log the values if you want
  #d <- d + scale_fill_gradient2(name="Predation rate", trans="log")
  d <- d + scale_fill_gradient2(name="Predation rate")
  d <- d + scale_y_discrete(name="Predator") +
    scale_x_continuous(name=expression(paste(Log[10]," Prey Mass (g)", sep="")))
  
  return(d)
}

# Plots the predation rates on each prey by predator
# Fixed for theta disaggregated by size
plotPredationRatesbyPrey4Dtheta <- function(model,sp=NULL)
{
  # Set up some short cuts
  f <- model$f[dim(model$f)[1],,]
  sv <- model$SearchVol
  N <- model$N[dim(model$N)[1],,]
  dwmat <- matrix(rep(model$dw,dim(N)[1]),nrow=dim(N)[1],byrow=TRUE)
  pk <- model$predkernel
  m2 <- model$M2[dim(model$M2)[1],,]
  theta <- model$param$theta
  
  # predation by each species, by predator mass on prey mass
  predratefull <- sweep(pk,c(1,2),(1-f)*sv*N*dwmat,"*")
  # Collapse over the predator mass to get total predation by each predator on each prey mass
  #predrate <- apply(predratefull,c(1,3),sum)[,model$idxGrid] # slower than the colsums / aperm method
  predrate <- predratefull[,,model$idxGrid]
  browser()
  dimnames(predrate) <- list(Predator = model$param$species$species,
                             PredMass = model$w,
                             PreyMass = model$w)
  
  # Turn into a dataframe and remove 0s
  preddf <- melt(predrate)
  preddf[preddf$value==0,"value"] <- NA
  # Add a column of log10 Mass to get the widths
  preddf <- cbind(preddf,log10PredMass = log10(preddf$PredMass),log10PreyMass = log10(preddf$PreyMass))
  # And change the name of Species to Predator
  #names(preddf)[names(preddf)=="Species"] <- "Predator"
  
  # name theta
  #dimnames(theta) = list(Predator=dimnames(theta)[[1]], Prey=dimnames(theta)[[2]])
  dimnames(theta) = list(Predator=model$param$species$species, Prey = model$param$species$species, PredMass = model$w, PreyMass = model$w)
  thetadf <- melt(theta)
  names(thetadf)[names(thetadf)=="value"] <- "Interaction"
  # Combine the dataframes so we get interaction and predation rates
  preypred <- join(preddf,thetadf,by=c("Predator", "PredMass", "PreyMass"))
  preypred <- cbind(preypred, scaledpred = preypred$Interaction * preypred$value)
  # free up some memory
  rm("theta")
  rm("thetadf")
  rm("preddf")
  
  #browser()
  if (is.null(sp)) sp <- as.character(model$param$species$species)
  # Plot what?
  # Lots of pred - prey plots?
  # Looks awesome! But takes a long time to plot
  #	 plotall <- ggplot(preypred[preypred$Prey %in% sp,]) + geom_tile(aes(x = log10PredMass, y = log10PreyMass, fill = scaledpred)) + facet_wrap(Predator ~ Prey)
  #plotall <- plotall + scale_fill_gradient2(name="Predation rate")
  #plotall <- plotall + scale_x_continuous(name=expression(paste(Log[10]," Predator Mass (g)", sep=""))) +
  #				 scale_y_continuous(name=expression(paste(Log[10]," Prey Mass (g)", sep="")))
  
  # Sum total predation from each predator
  preypredsum <- ddply(preypred, .(Prey, log10PreyMass, Predator), summarise, sumpred = sum(scaledpred, na.rm=T))
  
  
  #d <- ggplot(preypredsum) + geom_tile(aes(x=log10PreyMass,y=Predator,fill=(scaledpred))) + facet_wrap(~Prey)
  #if (!is.null(sp))
  d <- ggplot(preypredsum[preypredsum$Prey %in% sp,]) + geom_tile(aes(x=log10PreyMass,y=Predator,fill=(sumpred))) + facet_wrap(~Prey)
  # log the values if you want
  #d <- d + scale_fill_gradient2(name="Predation rate", trans="log")
  d <- d + scale_fill_gradient2(name="Predation rate")
  d <- d + scale_y_discrete(name="Predator") +
    scale_x_continuous(name=expression(paste(Log[10]," Prey Mass (g)", sep="")))
  
  return(d)
}




#*******************************************************************************


#*******************
# Yield + SSB plots
#*******************

plotYieldSSBCompare <- function(model, meantsteps = 10, yieldname = "Catch_8595", ssbname = "SSB_8595", log10 = TRUE, addRankings=TRUE){
  # Sort out Yield
  Yieldhat <- apply(model$Yield[(dim(model$Yield)[1]-meantsteps+1):dim(model$Yield)[1],],2,mean)/1e6
  Yieldobs <- model$param$species[,yieldname]
  ydf <- data.frame(species = model$param$species$species, hat = Yieldhat, obs = Yieldobs)
  # Sort out SSB
  Nhat <- apply(model$N[(dim(model$N)[1]-meantsteps+1):dim(model$N)[1],,],c(2,3),mean)
  SSBhat <- apply(sweep(model$psi * Nhat,2,model$w * model$dw,"*"),1,sum) /1e6
  SSBobs <- model$param$species[,ssbname]
  ssbdf <- data.frame(species = model$param$species$species, hat = SSBhat, obs = SSBobs)
  # Put SSB and Yield into single df
  y_ssb_df <- rbind(cbind(ssbdf,measure="SSB"),cbind(ydf,measure="Yield"))
  # Cut out the rows with NA
  y_ssb_df <- y_ssb_df[!is.na(y_ssb_df$obs),]
  
  # Plot these
  p <- ggplot(y_ssb_df) + geom_point(aes(x = (obs), y = (hat), colour = species)) + facet_wrap(~ measure, scales = "free")
  # Add a ref line
  p <- p + geom_abline(intercept=0,slope=1)
  if (log10)
    p <- p + scale_x_continuous(name = "Log 10 observed", trans="log10") + scale_y_continuous(name = "Log 10 predicted", trans="log10")
  else
    p <- p + scale_x_continuous(name = "Observed") + scale_y_continuous(name = "Predicted")
  
  if (addRankings){
    rankdfmelt <- ddply(melt(y_ssb_df, id.vars = c("species","measure")), .(measure, variable), function(x) data.frame(species = x$species, rank = rank(x$value)))
    rankdf <- cast(rankdfmelt, species + measure ~ variable, value = "rank")
    pr <- ggplot(rankdf) + geom_point(aes(x = (obs), y = (hat), colour = species)) + facet_wrap(~ measure, scales = "free")
    # Add a ref line
    pr <- pr + geom_abline(intercept=0,slope=1)
    pr <- pr + scale_x_continuous(name = "Observed rank") + scale_y_continuous(name = "Predicted rank")
    vplayout <- function(x,y)
      viewport(layout.pos.row=x,layout.pos.col=y)
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2,1)))
    print(p,vp=vplayout(1,1))
    print(pr,vp=vplayout(2,1))
  }
  #return(p)
}




#*******************************************************************************
