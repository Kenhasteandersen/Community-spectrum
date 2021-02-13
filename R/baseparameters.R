#' Get a standard parameter set
#'
#' @param wInf Array of asymptotic sizes to use
#' @param kappa Carrying capacity of the resource
#' @param h Factor for the maximum consumption rate
#'
#' @details See the source code for description of fields
#' 
#' @return A structure with parameters
#' 
#' @export
baseparameters  = function(wInf = 10^seq(log10(4),log10(30000),length.out = 18), kappa = 1.3,h = 20){

beta <- 100.0; # Pred prey mass ratio

param <- list()

nSpecies <- length(wInf);
param$nSpecies <- nSpecies;   # no$ of species
param$fGridexp <- 0.2;     # Expansion of grid
param$nGridPP <- 50;    # Grid cells in the resource spectrum
param$dt <- 0.2;        # Time step 0.1
param$tEnd <- 40;       # End time # originally 40
param$iSave <- 1;        # how often to save results (in iterations)

# Species-specific:
param$wInf <- t(wInf);#logspace(log10(10), log10(1000000), nSpecies)';
param$wMax <- 2*max(param$wInf);

# Growth parameters:
param$alpha <- 0.6;                    # Assimilation efficiency
param$n <- 3/4; 
param$h <- rep(h, nSpecies)
param$f0 <- 0.6; # Feeding level used for calculation of equ$ spectra$
param$k <- 0
param$ks <- 0.12*param$h;         # Activity
param$p <- 3/4;                  # Scaling of std$ metabolism
param$fc <- param$ks[1]/(param$alpha*param$h);

# Encounter:
param$q <- 0.8;                        # Scaling of search volume
param$beta <- beta;                     # Predation/prey mass ratio
param$sigma <- 1.3;                            # Width of size preference (1$3)
lambda <- 2+param$q-param$n;
alphae <- sqrt(2*pi)*param$sigma*param$beta^(lambda-2)*exp((lambda-2)^2*param$sigma^2/2);
param$gamma <- (param$f0*param$h / (alphae*kappa*(1-param$f0)));
param$v <- 1*matrix(1,nSpecies)  # Vulnerabilities

# Mortality
param$mu0prefactor <-2; #2      # Natural mortality
param$mu0exponent <- -1/4; # -1/4
param$muS0 <- 0;                        # Prefactor for starvation mortality

# Recruitment:
param$w0 <- 0.001                   # Weigth to recruit at
param$R0 <- 0/param$wInf                   # "Background" recruitment flux
param$rho <- param$wInf^(0)                # Egg survival
param$rho <- param$rho/param$rho[1]
param$alphaMature <- 0.25             # Fraction of wInf to mature at 0$25
param$eRepro <- 0.1                    # Eff$ of gonad production - 0$1, CC 0$05

# Resource spectrum
param$typeResource <- 1 # semi-chemostat
param$rR <- 4
param$kappaR <- kappa;
param$PPmin <- 0.001
param$kR <- -2-param$q+param$n;
param$lR <- 0.25;
param$wRcut <- 1;                  # Cut off of background spectrum
#
# Calculate Rmax for recruitment function:
#
n <- param$n;
q <- param$q;
# Physiological mortality
alphap <- param$f0 * param$h*param$beta^(2*n-q-1) * exp((2*n*(q-1)-q^2+1)*param$sigma^2/2)
param$a <- alphap / (param$alpha*param$h*param$f0-param$ks) 
# Rmax:
tmpA = param$wInf[1];
tmpB = (log10(param$wInf[param$nSpecies])-log10(param$wInf[1]))/(param$nSpecies-1)
dW = tmpB*tmpA*10^(tmpB*((1:param$nSpecies)-1))
param$Rmax = 1.7*param$kappaR*(param$alpha*param$f0*param$h*param$w0^param$n-param$ks*
                                 param$w0^param$p)*param$wInf^(2*param$n-param$q-3+param$a)*param$w0^(-param$n-param$a[1])*dW
#
# Make some reasonable initial conditions:
#
w <- makegrid(param)[[1]];
param$Ninit <- matrix(0,param$nSpecies, length(w))
for (i in 1:nSpecies) {
  param$Ninit[i,] <- w^(-param$n-param$a[i])
  param$Ninit[i, w>param$alphaMature*param$wInf[i]] = 0
  param$Ninit[i,] = 1 * param$Ninit[i,]/param$Ninit[i,1]*param$Rmax[i]/
    (param$alpha*param$f0*param$h[i]*param$w0^param$n-param$ks[i]*param$w0^param$p)
}
# ----------------------------------------------------
# Params for trawling function
# ----------------------------------------------------- 
param$F0 <- 0 * matrix(1,nSpecies); # Overall fishing pressure on

param$c <- 0.01
param$myF <- 10
#param$aF <- 0$2;
param$nF <- 0.05
param$gSigma <- 1.5 #

return(param)
}

#
# Parameters for Lake Kariba. Does not quite work yet.
#
paramKariba <- function() {
  param <- baseparameters(wInf = 10^seq(log10(20), log10(0.01*110^3), length.out = 20), h=28, kappa = 3.3e9)
  kappa = 3.3e9
  nSpecies <- param$nSpecies
  param$q <- 0.75
  lambda <- 2+param$q-param$n;
  alphae <- sqrt(2*pi)*param$sigma*param$beta^(lambda-2)*exp((lambda-2)^2*param$sigma^2/2);
  param$gamma <- (param$f0*param$h / (alphae*kappa*(1-param$f0)));
  param$kR <- -2-param$q+param$n;
  
  #
  # Calculate Rmax for recruitment function:
  #
  n <- param$n;
  q <- param$q;
  # Physiological mortality
  alphap <- param$f0 * param$h*param$beta^(2*n-q-1) * exp((2*n*(q-1)-q^2+1)*param$sigma^2/2)
  param$a <- alphap / (param$alpha*param$h*param$f0-param$ks) 
  # Rmax:
  tmpA = param$wInf[1];
  tmpB = (log10(param$wInf[param$nSpecies])-log10(param$wInf[1]))/(param$nSpecies-1)
  dW = tmpB*tmpA*10^(tmpB*((1:param$nSpecies)-1))
  param$Rmax = 1.5e13 #1.7*param$kappaR*(param$alpha*param$f0*param$h*param$w0^param$n-param$ks*
                       #            param$w0^param$p)*param$wInf^(2*param$n-param$q-3+param$a)*param$w0^(-param$n-param$a[1])*dW
  #
  # Make some reasonable initial conditions:
  #
  w <- makegrid(param)[[1]];
  param$Ninit <- matrix(0,param$nSpecies, length(w))
  for (i in 1:nSpecies) {
    param$Ninit[i,] <- w^(-param$n-param$a[i])
    param$Ninit[i, w>param$alphaMature*param$wInf[i]] = 0
    param$Ninit[i,] = 1 * param$Ninit[i,]/param$Ninit[i,1]*param$Rmax[i]/
      (param$alpha*param$f0*param$h[i]*param$w0^param$n-param$ks[i]*param$w0^param$p)
  }
  return(param)
}
