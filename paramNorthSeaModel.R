# Mlla's copy from Dr Julia Blanchard 02/10/2017
# ==============================================

# Function to set up parameter list for the NS model

paramNSModel <- function(fileSpecies,
                         theta,
                         kap=1e11)
{
#browser()
  # Use ... to set values?
  # Set up list
  # Set everything to NA to start with (or populate with defaults)
  # then pull out ... and overwrite?

  # theta now passed in as the matrix - not a file name

  # Read in species specific parameter from file:
  species = read.delim(fileSpecies)
  #species = read.csv(fileSpecies)

  # round the Winfs - false precision otherwise
  species$Winf <- round(species$Winf,1)

  # Set up parameter list
  param=list( species = species,                     # species parameters
              theta   = theta,  # interaction matrix

              # Model Parameters
              nspp    = dim(species)[1],             # Number of species
              ngrid   = 100,                         # no. of size classes
              ngridPP = 30,
              tmax    = 100,                         # no. years
              dt      = 1,                           # time step
              isave   = 1,                           # how often to save results
              w0      = 0.001,                           # mass of eggs i.e. minimum size in g
              wMax    = max(species$Winf)*1.1,           # max weight in grid
              slope0  = NA,                          # initial slope, calculated below

              # Primary production
              wPPcut  = 10,                          # cutoff size of background spectrum
              rPP     = 10,                          # Growth rate Primary productivity
              kap     = kap,                        # Carrying capacity of resource spectrum  (height of plankton)
              lambda  = NA,                          # Exponent of resource spectrum - calculated below
             
              # Growth
              n       = 2/3,                         # scaling of maximum intake (the prefactor h is calculated below, 
              									     # Hartvig et al. 2011 Theor Pop Biol use 3/4, Andersen & Pedersen 2010 Proc Roy Soc B use 2/3)
              p       = 0.7,                         # scaling of std metabolism
              SRR     = NA,                          # the stock recruitment function

              # Feeding
              q       = 0.8,                         #search volume exponent

              # Gamma calculation
              f0est   = 0.6,                         # equilibrium feeding level, for which h-bar was estimated

              # Intrinsic natural mortality
              Z0pre   = 0.6                         # intrinsic natural mort. prefactor
               )

  # After list calculations
    
  # Exponent for intrinsic natural mortality 
  # (Andersen & Pedersen 2010, Hartvig et al. 2011 both set this as  n-1, when n is 3/4 this is consistent with metabolic theory of Brown 2004, at the species-level)
  param$Z0exp <- param$n - 1              	     

  # From SI information IMAGE, the prefactor of maximum intake is caluclated as: 
  # h = ((3 * K) / (alpha * f0))  * W^(1/3)
  # where:
  # f0 = 0.6, initial feeding level
  # alpha = 0.6, assimilation efficiency
  # K = k_vb, species specific empirical von bertalanffy K parameters
  # W = Winf, species-specific asymptotic weight
  
  #alpha, the assimilation efficiency, is held fixed across species
  param$species$alpha <- 0.6
  
  # species-specific prefactor of maximum intake
  param$species$h <- ((3 * param$species$k_vb) / (param$species$alpha * param$f0est)) * param$species$Winf^(1/3)
  
  
  # Add standard metabolism and activity to species list
  param$species$ks <- 0.2 * param$species$h # Standard metabolism
  param$species$k <- rep(0,dim(species)[1]) # Activity

  # Exponent of resource spectrum
  param$lambda <- 2+param$q-param$n

  # Add search volume to species list 
  param$species$gamma <- (param$f0est*param$species$h*param$species$beta^(2-param$lambda)) /
                        ((1-param$f0est)*sqrt(2*pi)*param$kap*param$species$sigma)
  #species specific gamma estimated from kappa from equilibrium assumptions Andersen & Beyer 2006 Am Nat.168: 54-61

  # Recruitment
  param$species[,"eRepro"] <- 1
  #param$species$eRepro= rep(0.02, param$nspp)           # efficiency of gonad production
  #param$species$R0 = species$R0 #* species$alpha*species$h*param$f0*param$w0^param$n
                                # Boundary condition= Growth rate at egg size*Number density at egg size

  # Just set R0 to anything to start with - from JB
  # But we do need an equilibrium to start the calibration
  # This comes from Ken's theory (abundance and Winf scaling with -1.5
  # Need multiplier. Keep increasing until stable. Looks Ok
  # kap = 1e11 and 1e12, interaction = twostage2DnoSaithe
  # kap = 1e11 and 1e12, interaction = Schoener_twostage2D_noSaithe
  param$species$R0 <- 1000*param$kap * param$species$Winf^(-1.5)
  
  # only use the below equilibrium values of R0 if calibrating to kappa ( because the are dependent on param$kap)
  # give  R0 values the same as theoretical ones

  # alphap.alphah<-(param$species$h*param$species$beta^(2*param$n-param$q-1)*exp((2*param$n*(param$q-1)-param$q^2+1)*param$species$sigma^2/2))/(param$species$alpha*param$species$h)
  # kapsp<-param$kap*param$species[,"Winf"]^(2*param$n-param$q-3+(alphap.alphah))
  # RmaxAB<-(kapsp*param$species[,"Wmin"]^(-param$n-(alphap.alphah)))*(1-(param$species[,"Wmin"]/param$species[,"Winf"])^(param$p-param$n))^(((alphap.alphah)*(param$p-param$n)) -1)
  # param$species$R0<-RmaxAB

  # Used to set initial values
  param$species[,"N0"] <- param$species[,"Winf"]^(0.7-3+0.5+1)*1e10
  #param$species$N0 = 1e-13 * species$R0 / (param$alpha*param$f0est*param$species$h) * param$w0^(-param$n-1)
  # Initial spectrum  - just used to set initial abundances in SizeBasedModel.r
  param$slope0 <- -param$n - 0.5 # estimated initial slope

  # Set the SRR to be a Beverton Holt esque relationship
  # This should be species specific
  param$SRR <- function(SSB, RDI, param){return(param$species$R0 * RDI / (param$species$R0+RDI))}

  # Fishing parameters ands settings
  # Selectivity is set in the Setup function not here.
  # This is because we need w which is not specified until setup
  # Here we set the selectivity parameters: gear, function name and the parameter names and values
  # 

    # Knife edge selectivity (original)
    #param$sel_params <- data.frame(gear = 1,species=param$species$species,
    #                      selfunc = "knife_edge", param_name = rep(c("xi","wFstart"),each=param$nspp),
    #                      param_value = c(rep(0.01,param$nspp),c(4,5,10,100,165,115,165,115,260,175,500,988)))

    # Trawl (sigmoid shape - FAO book and Lembo paper
#    param$sel_params <- data.frame(gear = 1,species=param$species$species,
#                      selfunc = "trawl", param_name = rep(c("L25","L50","a","b"),each=param$nspp),
#		    param_value = c(param$species$L25,param$species$L50,param$species$a,param$species$b))

  # Different gear for each species - so we muck about with historical Fs later on
  #browser()
    param$sel_params <- data.frame(gear = param$species$species,
				   species=param$species$species,
                      selfunc = "trawl", param_name = rep(c("L25","L50","a","b"),each=param$nspp),
		    param_value = c(param$species$L25,param$species$L50,param$species$a,param$species$b))

    gears <- unique(param$sel_params$gear)

    # Set the catchability. This is a 2D array. Catchability is a scalar by stock and gear
    param$Q <- array(0,dim=c(param$nspp,length(gears)),dimnames=list(species=param$species$species,gear=gears))
    # set equal to F0 (so an effort of 1 gives max Fmort of F0)
    # Use mean F0 1985 - 1995
    #param$Q[] <- param$species$F0_8595
    # If using different gear for each species
    diag(param$Q) <- param$species$F0_8595


    # Effort is set by gear only
    # Set this in setup routine so it can change with time (i.e. add time dimemsion)
  return(param)
}


#************************************************************************
# Same as above but includes 4 'fleets'
paramNSModel_multifleet <- function(fileSpecies,
                         fileInteraction,
                         kap=1e11)
{


  # Read in species specific parameter from file:
  species = read.delim(fileSpecies)
  # tweak Winf
  species$Winf <- species$Winf * 1.1

  # Read interaction matrix file
  imat <- read.delim(fileInteraction)
  theta   <- as.matrix(imat[,-1])  # interaction matrix
  rownames(theta) <- names(imat)[-1]

  # Set up parameter list
  param=list( species = species,                     # species parameters
              theta   = theta,  # interaction matrix

              # Model Parameters
              nspp    = dim(species)[1],             # Number of species
              ngrid   = 100,                         # no. of size classes
              ngridPP = 30,
              tmax    = 100,                         # no. years
              dt      = 1,                           # time step
              isave   = 1,                           # how often to save results
              w0      = 0.001,                           # mass of eggs i.e. minimum size in g
              wMax    = max(species$Winf),           # max weight in grid
              slope0  = NA,                          # initial slope, calculated below

              # Primary production
              wPPcut  = 10,                          # cutoff size of background spectrum
              rPP     = 10,                          # Growth rate Primary productivity
              kap     = kap,                        # Carrying capacity of resource spectrum  (height of plankton)
              lambda  = NA,                          # Exponent of resource spectrum - calculated below

              # Growth
              n       = 2/3,                         # scaling of intake
              p       = 0.7,                         # scaling of std metabolism
              SRR     = NA,                          # the stock recruitment function

              # Feeding
              q       = 0.8,                          #search volume exponent

              # Gamma calculation
              f0est   = 0.6,                         # equilibrium feeding level, for which h-bar was estimated

              # Mortality
              Z0pre   = 0.6,                        # intrinsic natural mort. prefactor
              Z0exp   = (2/3) - 1                  # exponent for intrinsinc nat mort (here = n-1 as in Andersen & Perdersen 2010, but alternatively could be -1/4)

  )
  # After list calculations
  # From SI information IMAGE
  # h = ((3 K) / (alpha f0))  * W^(1/3)
  # f0 = 0.6
  # K = vb K
  # W = Winf


  # Add standard metabolism and activity to species list
  param$species$ks <- 0.2 * param$species$h # Standard metabolism
  param$species$k <- rep(0,dim(species)[1]) # Activity

  # adjust param$sp$h
  #param$species$h <- param$species$h * 1.5

  # Exponent of resource spectrum
  param$lambda <- 2+param$q-param$n

  # Add search volume to species list
  param$species$gamma <- (param$f0est*param$species$h*param$species$beta^(2-param$lambda)) /
                        ((1-param$f0est)*sqrt(2*pi)*param$kap*param$species$sigma)
  #species specific gamma estimated from kappa from equilibrium assumptions Andersen & Beyer 2006

  # Recruitment
  param$species[,"eRepro"] <- 1
  #param$species$eRepro= rep(0.02, param$nspp)           # efficiency of gonad production
  #param$species$R0 = species$R0 #* species$alpha*species$h*param$f0*param$w0^param$n
                                # Boundary condition= Growth rate at egg size*Number density at egg size
  param$species[,"N0"] <- param$species[,"Winf"]^(0.7-3+0.5+1)*1e10
  #param$species$N0 = 1e-13 * species$R0 / (param$alpha*param$f0est*param$species$h) * param$w0^(-param$n-1)

  # Initial spectrum
  param$slope0 <- -param$n - 0.5 # estimated initial slope

  # Set the SRR to be a Beverton Holt esque relationship
  # This should be species specific
  param$SRR <- function(SSB, RDI, param){return(param$species$R0 * RDI / (param$species$R0+RDI))}

  # Fishing parameters ands settings
  # Selectivity is set in the Setup function not here.
  # This is because we need w which is not specified until setup
  # Here we set the selectivity parameters: gear, function name and the parameter names and values
  #

  #********* Fleet bits ****************
  # Give everything a trawl shape - but set up 4 gears
  #browser()

  
    # Trawl (sigmoid shape - FAO book and Lembo paper
    param$sel_params <- data.frame(gear = 1,species=param$species$species,
                      selfunc = "trawl", param_name = rep(c("L25","L50","a","b"),each=param$nspp),
                      param_value = c(param$species$L25,param$species$L50,param$species$a,param$species$b))

    # Industrial gear: sprat, sandeel, norway pout
    gear1sp <- param$species$species[1:3]
    # Pelagic: herring
    gear2sp <- param$species$species[4]
    # Beam trawl: dab, sole, plaice, spurdog
    gear3sp <- param$species$species[c(5,7,9, 13)]
    # Otter demersal: gurnard, whiting, haddock, cod, saithe, wolffish
    gear4sp <- param$species$species[c(8, 6, 10, 11, 12, 14)]

    param$sel_params[param$sel_params$species %in% gear1sp,"gear"] <- 1
    param$sel_params[param$sel_params$species %in% gear2sp,"gear"] <- 2
    param$sel_params[param$sel_params$species %in% gear3sp,"gear"] <- 3
    param$sel_params[param$sel_params$species %in% gear4sp,"gear"] <- 4

    gears <- unique(param$sel_params$gear)

    # Set the catchability. This is a 2D array. Catchability is a scalar by stock and gear
    param$Q <- array(NA,dim=c(param$nspp,length(gears)),dimnames=list(species=param$species$species,gear=gears))
    # set equal to F0 (so an effort of 1 gives max Fmort of F0)
    # Use mean F0 1985 - 1995
    param$Q[] <- param$species$F0_8595

    # Effort is set by gear only
    # Might need to set this in setup routine so it can change with time (i.e. add time dimemsion)
    #param$effort <- array(1,dim=c(1,length(gears)),dimnames=list(all="all", gear=gears))

  # Redundant stuff?
  #param$alphae = sqrt(2*pi)*param$gamma*species$sigma*species$beta^(param$lambda-2)* exp((param$lambda-2)^2*species$sigma^2/2)
  #param$f0 = 1 / (1 + species$h/(param$kap*param$alphae[1]))

  # 'Diagnostics'
  #param$fc = param$ks/(species$alpha*species$h)     #critical feeding level - only enough food eaten to meet standard metabolism
  param$Volumecubicmetres=5.5e13    #unit of volume. Here total volume of North sea is used (Andersen & Ursin 1977)

  return(param)
}
