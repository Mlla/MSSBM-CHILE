# Mlla's copy from Dr Julia Blanchard 02/10/2017
# Selectivity functions and notes

# Selectivity reflects gear-fish interaction and therefore varies with age or size
# e.g. sigmoid shape
# goes between 0 and 1
# Catchability (Q) represents fishing 'conditions' e.g. season (spawning season?),
# location (on spawning grounds?), environmental conditions, skill of skipper etc.
# Effort is effort
# Fishing mortality varies with age-size
# The relationship is assumed to be:
# F = E * Q * Sel
# (F and Sel are vectors with size)

# Here E=1 in the base year. Q can just be absorbed by scaling of Sel?
# E and Q are dynamic and may change with time
# Sel is between 0 and 1, so selectivity is modified by Q and effort to give F

# Selectivity functions
# Knife edge (original)

knife_edge <- function(w,wFstart,xi)
  return(1/(1+exp(((wFstart-w)/xi)*wFstart)))

# Test
#sel_knife <- knife_edge(w=NSmodel$w,wFstart=NSmodel$param$species$wFstart[12],xi=NSmodel$param$xi)
#plot(x=NSmodel$w,y=sel_knife,log="x",type="l")

#sel_func <- "knife_edge"
#sel_knife <- do.call(knife_edge,list(w=NSmodel$w,wFstart=NSmodel$param$species$wFstart[12],xi=NSmodel$param$xi))
#sel_knife <- do.call("knife_edge",list(w=NSmodel$w,wFstart=NSmodel$param$species$wFstart[12],xi=NSmodel$param$xi))
#sel_knife <- do.call(sel_func,list(w=NSmodel$w,wFstart=NSmodel$param$species$wFstart[12],xi=NSmodel$param$xi))

#******************************************************************************

# From FAO and Lembo
trawl_length <- function(L,L25,L50)
{
  SR <- L50 - L25
  S1 <- L50*log(3)/SR
  S2 <- S1 / L50
  return(1 / (1 + exp(S1 - S2*L)))
}
# test - see FAO p187
#FAO use L50 and L75. As L25 and L75 are symmetrical around L50, we can use L25
#L <- 0:25
#L50 <- 13.2
#L75 <- 14.7
#L25 <- L50-(L75-L50)
#sel <- trawl_length(L=L, L25=L25, L50=L50)
#plot(L,sel,type="l")
#lines(x=c(0,100),y=c(0.5,0.5),lty=2)
#lines(x=c(L50,L50),y=c(0,1),lty=2)
#lines(x=c(L75,L75),y=c(0,1),lty=2)
#lines(x=c(L25,L25),y=c(0,1),lty=2)
#

# What does this mean for selectivity at weight?
# w <- a * L ^b
trawl <- function(w,L25,L50,a,b)
{
  L <- (w/a)^(1/b)
  return(trawl_length(L,L25,L50))
}

# CHECK a in species list
#sp <- 12
#w <- NSmodel$w
#a <- NSmodel$param$species$a[sp]
#b <- NSmodel$param$species$b[sp]
#L <- (w/a)^(1/b)
#L25 <- 25
#L50 <- 30
#w25 <- a*L25^b
#w50 <- a*L50^b
#sel <- trawl(w,L25,L50,a,b)
#par(mfrow=c(2,1))
#plot(L,sel,type="l")
#lines(x=c(L50,L50),y=c(0,1),lty=2)
#lines(x=c(L25,L25),y=c(0,1),lty=2)
#plot(w,sel,log="x",type="l")
#lines(x=c(w50,w50),y=c(0,1),lty=2)
#lines(x=c(w25,w25),y=c(0,1),lty=2)

#******************************************************************************
# From Pope et al 2006
# Popes values are:
# gamma (slope) = 0.2
# delta (50% selected length at delta proportion of Linf) = 0.33
# Linfcentre = 70cm - used to estimate Q

pope_logistic_length <- function(L,Linf,gamma,delta)
{
  #lambda <- 1 + kappa * (Linf - Linfcentre)
  sel <- 1 / (1 + exp(gamma * (delta * Linf - L)))
  return(sel)
}

pope_logistic <- function(w,a,b,Linf,gamma,delta)
{
  L <- (w/a)^(1/b)
  return(pope_logistic_length(L,Linf,gamma,delta))
}

#delta <- 0.33
#gamma <- 0.2
##kappa <- -0.0035
#L <- 0:150
#Linf <- 13
#selL <- pope_logistic_length(L,Linf,gamma,delta)
#a <- 0.0439
#b <- 2.91
#w <- a * L ^b
#plot(l,selL,type="l",log="x")
#
#w <- 0:10000 #g
#l <- (w/a)^(1/b)
#sel <- pope_logistic_length(l,Linf,gamma,delta)
#plot(w,sel,type="l",log="x")
#
## 1 g
#l <- (1/a)^(1/b)
## a 1g BOF is 3cm long - seems dodgy - hence selectivity is too high
#pope_logistic_length(l,Linf,gamma,delta)
## Plot looks out - check selectivity calculation
