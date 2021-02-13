
baserun <- function(nSpecies = 18 ,F0 = matrix(0.3,3,2), S= NA, Parameterset = 'Generic', wsize = c(100,3000)){
  #
  # Set the parameter set
  #
  if (Parameterset == 'Generic'){
    wInf = c(4,50000)
    W <- 10^seq(log10(wInf[1]),log10(wInf[2]),length.out = nSpecies) # 10 species in logspace
    param <- baseparameters(W,kappa = 1.3,h = 20)
  }
  
  if (Parameterset == 'North Sea')
    data(NorthSea)

  if (Parameterset == 'Benguela Current')
    data(Benguela)
  
  if (Parameterset == 'Baltic Sea')
    data(Baltic)

  if (Parameterset == 'Northeast US Cont. Shelf')
    data(NEUSCS)

  if (Parameterset == 'Barents Sea')
    data(Barents)

  W <- param$wInf
  nSpecies <- param$nSpecies
  #
  # Set fishing:
  #
  param$F0 <- rep(NA,nSpecies)
  param$F0[W <= wsize[1]] <- F0[1,1]
  param$F0[W <= wsize[2] & W > wsize[1]] <- F0[2,1]
  param$F0[W > wsize[2]] <- F0[3,1]
  
  param$fishing <- "Trawl" # See "fishing.R" for definitions
  #
  # Run the "before" simulation
  #
  param$tEnd <- 40
  param$dt <- 0.5 # run the model faster
  SF <- IterateSpectrum(param,S = S) # Add S here to start at initial conditions from before
  #
  # run the "after" simulation
  #
  param$F0[W <= wsize[1]] <- F0[1,2]
  param$F0[W <= wsize[2] & W > wsize[1]] <- F0[2,2]
  param$F0[W > wsize[2]] <- F0[3,2]
  
  SF2 <- IterateSpectrum(param,SF)

  return(list(SFpre = SF,SFpost = SF2,param = param))
}



