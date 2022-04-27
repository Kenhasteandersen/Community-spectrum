#' Run the community spectrum simulator manually
#'
#' @param nSpecies Number of asymptotic size groups ("species").  Only relevant for trait-based model.
#' @param F0 Array with fishing mortality on each species
#' @param S Output from a previous simulation to use as initial conditions
#' @param param Set of parameters to use
#'
#' @return A structure with the output of the simulation with the main fields being:
#' @return  - t: The time steps where the solution is saved
#' @return  - nSpecies: no. of asymptotic size groups
#' @return  - w: the weight classes
#' @return  - g: growth rate as a function of: time, species, weight
#' @return  - f: feeding level (consumption divided by max consumption) as a function of: time, species, weight
#' @return  - M2R: predation mortality on resource
#' @return  - M2, Ms, Fin: mortality from predation, starvation and fishing
#' @return  - Ntot: community spectrum as a function of time and weight
#' @return  - N: species spectra as a function of time, species and weight
#' @return  - SSB: spawning stock biomass as a function of time and species
#' @return  - Yield: fisheries yield as a function of time and species
#'
#' @examples runTraitBasedModel()
#' 
#' 
#' @export
runTraitBasedModel <- function(nSpecies = 18 ,
                               F0 = 0.3+0*(1:nSpecies), 
                               S= NA, 
                               param = baseparameters(10^seq(log10(4),log10(30000),length.out = nSpecies),kappa = 1.3,h = 20))
{
  param$F0 <- F0
  param$fishing <- "Trawl" # See "fishing.R" for definitions
  
  S <- IterateSpectrum(param,S = S)
  S$param <- param
  
  return(S)
}

#' Run the food web size spectrum model manually
#'
#' @param Parameterset Text string with the name of the parameter set to use
#' @param F0 Array with fishing mortality on each fleet, or an array of fishig mortality on each species.
#' @param S Output from a previous simulation to use as initial conditions
#'
#' @return A structure with the output of the simulation
#'
#' @examples runFoodWebModel('North Sea', F0 = c(0.1, 0.3, 0.7))
#' @examples runFoodWebModel('Northeast US Cont. Shelf')
#' @examples runFoodWebModel('Benguela Current')
#' @examples runFoodWebModel('Baltic Sea')
#' @examples runFoodWebModel('Barents Sea')
#' 
#' @export
runFoodWebModel <- function(Parameterset='TraitBased', F0 = c(0.3, 0.3, 0.3), S = NA) 
{
  #
  # Choose the parameter set
  #
  if (Parameterset == 'Traitbased'){
    W <- 10^seq(log10(wInf[1]),log10(wInf[2]),length.out = nSpecies) # 10 species in logspace
    param <- baseparameters(W,kappa = 0.005,h = 20)
  }
  
  if (Parameterset == 'North Sea')
    data(NorthSea)
  
  if (Parameterset == 'Benguela Current')
    data('Benguela')
  
  if (Parameterset == 'Baltic Sea')
    data(Baltic)
  
  if (Parameterset == 'Northeast US Cont. Shelf')
    data(NEUSCS)
  
  if (Parameterset == 'Barents Sea')
    data(Barents)
  
  if (Parameterset == 'Kariba')
    param <- paramKariba()
  
  W <- param$wInf
  nSpecies <- param$nSpecies
  
  if (length(F0)==3) {
    param$F0 <- rep(NA,nSpecies)
    param$F0[W <= 100] <- F0[1]
    param$F0[W <= 3000 & W > 100] <- F0[2]
    param$F0[W > 3000] <- F0[3]
  } else {
    param$F0 <- F0
  }
  
  param$fishing <- "Trawl" # See "fishing.R" for definitions
  
  param$tEnd <- 40
  S <- IterateSpectrum(param,S = S) # Add S here to start at initial conditions from before
  S$param <- param
  
  return(S)
}
## side:  1 for x-axis, 2 for y-axis (3, 4 secondary x and y axes)
## labels: logical, if TRUE, adds them
## las: integer in [0,1,2,3] see ?par
logaxes <- function(side = 1, pow = -1:1, labels=TRUE, las = 1, col = 1) {
  at <- as.vector(sapply(seq(range(pow)[1], range(pow)[2]), function(p) (2:9)*10^p))
  axis(side, at = 10^pow, lwd=0, labels = FALSE, tcl=0.5, lwd.ticks = 1, col = col)
  if(labels) axis(side, at = 10^pow, lwd=0, labels = FALSE, tcl=-0., lwd.ticks = 1)
  if(labels) {
    for(i in pow) {
      mtext(side = side, at = 10^i, text=bquote(10^.(i)), line = 0.5, las = las, cex=0.6)  
    }
  }
  axis(side, at = at, labels = FALSE, tcl=0.25, lwd = 0, lwd.ticks=1, col = col)
} 

#' Plot the results from the food web model or the trait-based model.
#'
#' \strong{Left column}.
#' The left column shows output as a function of body weight of individuals.
#' Top panel: biomass spectra. Thick solid line is the community spectrum, thin lines
#' are the individual species or asymptotic size groups. The dashed line is the theoretical
#' expectation. The thick dashed line is the resource spectru.
#' Middle panel: feeding level, i.e. consumption rate divided by maximum consumption rate.
#' The dashed line is the theoretical expectation and the dotted line is the level
#' where the fish would starve.
#' Lower panel: mortality. Lines show predation mortality (black) and fishing mortality (red). 
#' Dashed line is the theoretical expectation. Blue circles is the constant background 
#' mortality on each species or asymptotic size group.
#' \strong{Right column}.
#' The right column shows out as a functio of the asymptotic size of species
#' (from the food web model) or asymptotic size groups (from the trait-based model).
#' Top panel: Spawning stock biomass
#' Middle panel: Life time reproroductive output per egg in the absence of stock-recruitment
#' relation. This is a measure of the recruitment of the species or asymptotic size group: 
#' a value below 1 indicates that the stock is collapsed (dotted line).
#' Lower panel: fisheries yield.
#'
#' @param S Output from a previous simulation to plot
#'
#' @return -
#'
#' @examples plotResults( runTraitBasedModel() )
#' 
#' @export
plotResults <- function( S ) 
{
  param <- S$param
  idxEnd <- param$tEnd/param$dt
  w <- S$w
  W <- param$wInf
  xlimit <- c(0.02, max(w))
  
  par(mfcol = c(3, 2), 
      cex = 0.6,
      mar = c(0, 4, 0, 0), 
      oma = c(.2, 0, 0.5, 0.5),
      tcl = -0.25,
      mgp = c(2, 0.6, 0))
  #layout(mat=matrix(c(1,1,2,3), 2, 2,byrow=TRUE),
  #       widths=c(3,1), heights=c(1,2))
  # -------------------------------
  # Size spectrum
  # -------------------------------
  B <- S$Ntot[idxEnd,]*w
  B[B == 0] <- NA # For plotting on log scale 
  
  yl <- param$kappaR * c(0.1*max(w)^(1+param$kR) , 0.02^(1+param$kR))
  #
  # Community spectrum:
  #
  plot(w,B, log = 'xy', type = 'l', col = 'black', lwd=3,
       ylab = 'Biomass density (-)', 
       xlab = 'Weight (g)', 
       ylim=yl, xlim=xlimit,
       axes=FALSE)
  logaxes(1, -2:5, labels=FALSE)
  logaxes(2, seq(floor(min(log10(yl))), ceiling(max(log10(yl))),by=2))
  box()
  #
  # Theoretical solution:
  #
  lines( w, param$kappaR*w^(1+param$kR), lty=2)
  #
  # Species spectra
  #
  N <- S$N[idxEnd,,]
  for (i in 2:param$nSpecies)
    lines(w,N[i,]*w, col = alpha('black',alpha = 0.3))
  #
  # Resource
  #
  wPP <- S$wPP
  lines(wPP, S$nPP[idxEnd,1,]*wPP, col='black', lwd=3, lty=2)
  # -------------------------------
  # Feeding levels
  # -------------------------------
  f <- S$f[idxEnd,,]
  plot(w, rep(param$fc[1],length(w)), type = 'l', log="x", lty=3,  # critical feeding level
       ylab = 'Feeding level',
       ylim=c(0,1),
       xlim= xlimit,
       axes=FALSE)
  lines(w, rep(param$f0, length(w)), lty=2) # Theoretical feeding level
  for (i in 2:param$nSpecies) {
    ix <- w <= W[i]
    lines(w[ix], f[i,ix], lty=1, lwd=3)
  }
  logaxes(1, -2:5, labels=FALSE)
  axis(2)
  box()
  # -------------------------------
  # Mortality
  # -------------------------------
  M2 <- S$M2[idxEnd,,]
  #mar = c(4, 4, 0, 0)
  par(mar = c(3, 4, 0, 0))
  plot(w, mean(param$a)*param$alpha*mean(param$h)*(param$f0-mean(param$fc)) * w^(param$n-1), type='l',lty=2,log="x",
       ylab = 'Mortality (1/yr)',
       xlab = "Weight (g)",
       ylim = c(0,4),
       xlim = xlimit, axes=FALSE)
  for (i in 2:param$nSpecies) {
    ix <- w<=param$wInf[i]
    lines(w[ix], M2[i,ix], type='l', lty=1, col=alpha('black', alpha =  0.5))
    lines(w[ix], S$Fi[i,ix], lty=1,col=alpha('red', alpha =  0.5))
    points(param$wInf[i], S$Z0[i], col='blue')
  }
  logaxes(1, -2:5)
  axis(2)
  box()
  
  
  # -------------------------------
  # SSB
  # -------------------------------
  SSB <- S$SSB[idxEnd,]
  par(mar = c(0, 4, 0, 0))
  plot(W, SSB, type='l', log='xy', lty=1,
       ylab = 'SSB (ton)', axes=FALSE)
  points(W, SSB)
  logaxes(1, min(floor(log10(param$wInf))):max(ceiling(log10(param$wInf))), labels=FALSE)
  logaxes(2, min(floor(log10(SSB))):max(ceiling(log10(SSB))))
  box()
  # -------------------------------
  # R0
  # -------------------------------
  R0 <- S$Rp[idxEnd,]/S$R[idxEnd,]
  plot(W, R0, type='l', log='xy', lty=1,
       ylab = 'R0',
       ylim = c(0.5,max(R0)), axes=FALSE)
  points(W, R0)
  lines(W, rep(1,length(W)), lty=3)
  logaxes(1, min(floor(log10(param$wInf))):max(ceiling(log10(param$wInf))), labels=FALSE)
  logaxes(2, -1:max(ceiling(log10(R0))))
  box()
  
  # -------------------------------
  # Yield
  # -------------------------------
  par(mar = c(3, 4, 0, 0))
  plot(W, S$Yield, type='l', log='xy', lty=1,
       ylab = 'Yield (ton/yr)',
       xlab = "Asymptotic weight (g)", axes=FALSE)
  points(W, S$Yield)
  logaxes(1, min(floor(log10(param$wInf))):max(ceiling(log10(param$wInf))))
  logaxes(2, min(floor(log10(S$Yield))):max(ceiling(log10(S$Yield))))
  box()
}


alpha <- function(col, alpha) {
  adjustcolor(col, alpha.f = alpha)
}

#' Shiny ui
#'
#' @return
#'
#' @examples 
#' 
#' @export
shinyui <- fluidPage(
  # Google analytics tracking
  HTML(paste("<script>
             (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
             (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
             m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
             })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
             
             ga('create', 'UA-96462600-1', 'auto');
             ga('send', 'pageview');
             
             </script>")),
  
  #
  # Header
  #
  h1('Community size spectrum simulator'),
  p('Calculate the expected ecosystem effect of a management plan involving changing the fishing mortality 
    on one aspect of the fish community from an "initial" fishing pattern to "new" fishing pattern. '),
  fluidRow(
    column(1),
    column(1,actionButton(inputId = 'click', label = 'Start simulation')))
  ,
  p(' - wait a few seconds and scroll down to see results of the simulation.')
  ,
  #
  # Inputs
  #
  wellPanel(tabsetPanel(
    tabPanel('Define ecosystem',
             p('Select an ecosystem. Use "generic" for running the trait-based model which models a generic fish community.'),
             fluidRow(
               column(2, br(), p('Ecosystem')),
               column(2, selectInput(inputId = 'Parameterset', label = '',
                                     choices =  c('Generic', 'North Sea', 'Baltic Sea', 
                                                  'Benguela Current', 'Northeast US Cont. Shelf',
                                                  'Barents Sea'))),
               column(8, htmlOutput("EcosystemDescription") )
             ))
    ,
    tabPanel('Fishing mortalities',
             fluidRow(
               column(2, br(),br(), p('Forage fleet')),
               column(3,sliderInput(inputId = 'Fsmall', label = '   Initial fishing level', value = 0.5, min = 0,
                                    max = 3, step = 0.1)),
               column(3,sliderInput(inputId = 'Fsmall.after', label = '   New fishing level', value = 0.5, min = 0,
                                    max = 3, step = 0.1))
             )
             ,
             fluidRow(
               column(2, br(),br(), p('Pelagic fleet')),
               column(3,sliderInput(inputId = 'Fmedium', label = '', value = 0.5, min = 0,
                                    max = 3, step = 0.1)),
               column(3,sliderInput(inputId = 'Fmedium.after', label = '', value = 0.5, min = 0,
                                    max = 3, step = 0.1))
             )
             ,
             fluidRow(
               column(2, br(),br(), p('Demersal fleet')),
               column(3,sliderInput(inputId = 'Flarge', label = '', value = 0.5, min = 0,
                                    max = 3, step = 0.1)),
               column(3,sliderInput(inputId = 'Flarge.after', label = '', value = 0.5, min = 0,
                                    max = 3, step = 0.1))))
               ,
    tabPanel('Define fishing fleets',
             p('There are three fleets: 1) forage fleet targeting small species; 
      2) pelagic fleet tageting medium sized species, and 3) a demersal fleet targeting large species. 
      The fleets are defined by largest fish targeted by the forage fish fleet, 
      and the smallest size targeted by the demersal fleet:'),
             fluidRow(
               column(2),
               column(3,sliderInput(inputId = 'wMiddle', label = 'Max. weight of forage fish', value = 100, min = 0,
                                    max = 1000, step = 50)),
               column(3,sliderInput(inputId = 'wLarge', label = 'Min. weight of large fish', value = 3000, min = 1000,
                                    max = 10000, step = 500))
             ))))
  ,
  #
  # Plots:
  #
  h3('Simulation results')
  ,
  fluidRow(
    column(6,
           plotOutput(outputId = 'plotBiomass')),
    column(6,
           plotOutput(outputId = 'plotF')))
  ,
  fluidRow(
    column(6,
           plotOutput(outputId = 'plotYield')),
    column(6,
           plotOutput(outputId = 'plotSSB')))
  ,
  fluidRow(plotOutput(outputId = 'plotSpectrum', height = 600))
  ,
  wellPanel(
    p('The calculations are based on the model in ',
      a('Andersen et al (2016).', href='https://www.researchgate.net/publication/284514316_The_theoretical_foundations_for_size_spectrum_models_of_fish_communities'), 
      'The calibrations for the specific ecosystem are described in ',
      a('Jacobsen et al (2016).', href='https://www.researchgate.net/publication/305110273_Efficiency_of_fisheries_is_increasing_at_the_ecosystem_level?_iepl%5BviewId%5D=dbdozbH76RI7ERE6VF6iTz8O&_iepl%5BprofilePublicationItemVariant%5D=default&_iepl%5Bcontexts%5D%5B0%5D=prfpi&_iepl%5BinteractionType%5D=publicationTitle'),
      'Code by ',
      a("Nis Sand Jacobsen", href='mailto:nisjac@uw.edu'),
      ' and ',
      a('Ken Haste Andersen.', href = 'mailto:kha@aqua.dtu.dk'),
      'Output of the model should not be used in a practical management setting without first consulting the authors. '),
    p('The model can be downloaded and run from git:'),
    p('install.packages(\'devtools\')'),
    p('install.packages(\'shiny\')'),
    p('devtools::install_git(\'https://github.com/Kenhasteandersen/Community-spectrum.git\') '),
    p('library(cspectrum)'),
    p('runCommunitySpectrum()')
    #,
    #img(src="Logo.png")
  )
)

#' Shiny ui
#'
#' @param input
#' @param output
#'
#' @return
#'
#' @examples 
#' 
#' @export
shinyserver <- function(input,output){
  #
  # Define fleet sizes:
  #
  wMiddle <- 100
  wLarge <- 3000
  #
  # Make the description of the selected ecosystem:
  #
  output$EcosystemDescription <- renderText({
    
    if (input$Parameterset=='Generic')
      return("The generic trait-based model containing species with asymptotic sizes between 4 g and 50 kg.")
    
    if (input$Parameterset=='North Sea')
      return("North Sea ecosystem. Contains the following species:<br> 
             sandeel, sprat, norway pout, herring, sole, plaice, haddock, pollock, cod, whiting.")
    
    if (input$Parameterset == 'Benguela Current')
      return("Benguela current ecosystem. Contains the following species:<br>
             anchovy, sardine, kingklip, shallow-water hake, deep-water hake.")
    
    if (input$Parameterset == 'Baltic Sea')
      return("Baltic Sea ecosystem. Contains the following species:<br>
             sprat, herring, cod.")
    
    if (input$Parameterset == 'Northeast US Cont. Shelf')
      return("Northeast us continential shelf ecosystem. Contains the following species:<br>
             Atlantic butterfish, herring, menhaden, yellowtail flounder, acadian redfish, witch flounder, 
             Atlantic croaker, winter flounder, black sea bass, american plaice, weakfish
             spiny dogfish, bluefish, summer flounder, haddock, monkfish, cusk, tilefish, pollock,
             white hake, striped bass, atlantic cod, atlantic halibut.")
    
    if (input$Parameterset == 'Barents Sea')
      return("Barents sea ecosystem. Contains the following species:<br>
             capelin, pollock, golden redfish, greenland halibut, haddock, atlantic cod.")
  })
  #
  #   Run simulation when the button is clicked
  #
  simResults <- eventReactive(input$click,{
    #
    # Set sizes of fleets:
    #
    wMiddle <- input$wMiddle
    wLarge <- input$wLarge
    if (wMiddle > wLarge){
      wMiddle <- wLarge}
    wsize <- c(wMiddle, wLarge)
    #
    # Set fishing mortalities:
    #
    F0 <- matrix(NA,3,2)
    F0[1,1] <- input$Fsmall
    F0[2,1]<- input$Fmedium
    F0[3,1]<- input$Flarge
    F0[1,2]<- input$Fsmall.after
    F0[2,2]<- input$Fmedium.after
    F0[3,2]<- input$Flarge.after
    #
    # Run simulation:
    #
    SF <- baserun(
      nSpecies = 27,
      F0 = F0,
      S= NA, Parameterset = input$Parameterset, wsize)
    
    return(SF)
  })
  
  
  output$plotBiomass <- renderPlot(
    {
      
      SF <- simResults()[[1]] # Run before
      SF2 <- simResults()[[2]]  # Run after
      param <- simResults()[[3]] # params
      wMiddle <- input$wMiddle
      wLarge <- input$wLarge
      
      idx.biomass <- which(names(SF) == 'SSB')
      Biomass <- SF[[idx.biomass]]
      
      idx.time <- which(names(SF) == 't')
      time <- SF[[idx.time]]
      
      # Sum biomass in small medium and large
      if (length(which(param$wInf < wMiddle)) > 1){
        Biomass.small <- rowSums(Biomass[,which(param$wInf < wMiddle)])
      }else{Biomass.small <- Biomass[,which(param$wInf < wMiddle)]}
      
      if (length(which(param$wInf <= wLarge & param$wInf > wMiddle)) > 1){
        Biomass.medium <- rowSums(Biomass[,which(param$wInf <= wLarge & param$wInf > wMiddle)])
      }else{Biomass.medium <- Biomass[,which(param$wInf <= wLarge & param$wInf > wMiddle)]}
      
      if (length(which(param$wInf > wLarge)) > 1){
        Biomass.large <- rowSums(Biomass[,which(param$wInf > wLarge)])
      }else{Biomass.large <- Biomass[,which(param$wInf > wLarge)]}
      
      if(length(Biomass.small) == 0){
        Biomass.small <- rep(NA, length(time))
      }
      
      idx.biomass <- which(names(SF2) == 'SSB') # Overwrite and plot the after new fishing
      Biomass <- SF2[[idx.biomass]]
      # Sum biomass in small medium and large
      if (length(which(param$wInf < wMiddle)) > 1){
        Biomass.small.a <- rowSums(Biomass[,which(param$wInf < wMiddle)])
      }else{Biomass.small.a <- Biomass[,which(param$wInf < wMiddle)]}
      
      if (length(which(param$wInf <= wLarge & param$wInf > wMiddle)) > 1){
        Biomass.medium.a <- rowSums(Biomass[,which(param$wInf <= wLarge & param$wInf > wMiddle)])
      }else{Biomass.medium.a <- Biomass[,which(param$wInf <= wLarge & param$wInf > wMiddle)]}
      
      if (length(which(param$wInf > wLarge)) > 1){
        Biomass.large.a <- rowSums(Biomass[,which(param$wInf > wLarge)])
      }else{Biomass.large.a <- Biomass[,which(param$wInf > wLarge)]}
      
      if(length(Biomass.small.a) == 0){
        Biomass.small.a <- rep(NA, length(time))
      }
      # Plot the time varying biomass
      
      # axis limits
      minYl <- min(c(Biomass.small,Biomass.medium,Biomass.large,Biomass.small.a,Biomass.medium.a,Biomass.large.a), na.rm = T)
      maxYl <- max(c(Biomass.small,Biomass.medium,Biomass.large,Biomass.small.a,Biomass.medium.a,Biomass.large.a), na.rm = T)
      
      title <- 'Development in spawning stock biomass'
      plot(time-time[length(time)], Biomass.small, log = 'y', xlab = 'Time (years)', ylab = 'Biomass (ton)', type = 'l', ylim =c(minYl,maxYl),
           xlim = c(-20, time[length(time)]), col = alpha('black', alpha =  0.5), main = title)
      if(length(Biomass.medium > 0)){
        lines(time-time[length(time)],Biomass.medium, col = alpha('black', alpha =  0.5), lwd = 2)
      }
      lines(time-time[length(time)],Biomass.large, col = alpha('black', alpha =  0.5), lwd = 3)
      
      #time2 <- seq(time[length(time)]+param$dt,2*time[length(time)], length.out = length(time))
      lines(time,Biomass.large.a, col = alpha('red', alpha =  0.3), lwd = 3)
      lines(time,Biomass.small.a, col = alpha('red', alpha =  0.3))
      if(length(Biomass.medium > 0)){
        lines(time,Biomass.medium.a, col = alpha('red', alpha =  0.3), lwd = 2)
      }
      legend('bottomright', legend = c('Small', 'Medium', 'Large'), lty = c(1,1,1), lwd = c(1,2,3), bty = 'n')
      # Add diagonal line
      lines(rep(0, 100), seq(min(minYl), max(maxYl), length.out = 100), lty = 2, col = 'black', lwd = 3)
    })
  
  
  output$plotF <- renderPlot({
    SF <- simResults()[[1]]
    param <- simResults()[[3]]
    SF2 <- simResults()[[2]]
    
    idx.fishing <- which(names(SF) == 'Fin')
    fishing <- SF[[idx.fishing]]
    
    idx.fishing <- which(names(SF2) == 'Fin')
    fishing.after <- SF2[[idx.fishing]]
    
    yl <- c(0,max(c(fishing,fishing.after)))
    title <- 'Fisheries selectivity'
    plot(SF$w,fishing[1,], log = 'x', type = 'l', lwd = 2, xlim = c(0.1,max(param$wInf)),
         col = alpha('black',alpha = 0.5), ylab = 'Fishing mortality (per year)',
         xlab = 'Weight (g)', ylim = yl, main = title)
    lines(SF2$w, fishing.after[1,], lwd = 2, col = alpha('red', alpha =  0.3))
    for (i in 2:param$nSpecies){
      ix = SF$w < param$wInf[i]
      lines(SF$w[ix], fishing[i,ix], lwd = 2, col = alpha('black',alpha = 0.5))
      lines(SF2$w[ix], fishing.after[i,ix], col = alpha('red', alpha =  0.5))
    }
    legend('topleft', legend = c('Before', 'After'), lty = c(1,1), col = c('black','red'), bty = 'n')
    
  })
  
  output$plotYield <- renderPlot({
    
    SF <- simResults()[[1]]
    SF2 <- simResults()[[2]]      #param <- MM[[which(MM == 'param')]]
    param <- simResults()[[3]]
    
    idx.Yield <- which(names(SF) == 'Yield')
    Yield <- SF[[idx.Yield]]
    idx.Yield <- which(names(SF2) == 'Yield')
    Yield2 <- SF2[[idx.Yield]]
    
    # Sum biomass in small medium and large
    yl <- c(min(c(Yield[Yield>0],Yield2[Yield2>0])),
            max(c(Yield[Yield>0],Yield2[Yield2>0])))
    
    if (input$Parameterset == 'Generic'){
      plot(param$wInf,Yield, log = 'xy', col = alpha('black',alpha = 0.5), type = 'l', 
           xlab = 'Asymptotic weight (g)', 
           ylab = 'Yield (ton/km2/yr)',
           ylim = yl,
           xlim = c(min(param$wInf),max(param$wInf)),
           main = 'Yield (ton/km2/yr)', lwd = 3)
      lines(param$wInf,Yield2, col = alpha('red', alpha = 0.3), lwd = 3)
      lines(rep(input$wMiddle,100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      lines(rep(input$wLarge, 100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      legend('bottomright', legend = c('Before', 'After'), lty = c(1,1), 
             col = c(alpha('black', alpha = 0.5),'red'), 
             bty = 'n')
      fmt <- '%3.2f '
    }else{
      plot(param$wInf,Yield, log = 'xy', col = alpha('black',alpha = 0.5),
           xlab = 'Asymptotic weight (g)', 
           ylab = 'Yield (ton/yr)',
           ylim = yl,
           main = 'Yield (ton/yr)', pch = 16, lwd = 3, cex = 2)
      points(param$wInf,Yield2, col = alpha('red', alpha = 0.5), cex = 2, lwd = 3)
      lines(rep(input$wMiddle,100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      lines(rep(input$wLarge,100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      legend('bottomleft', legend = c('Before', 'After'), pch = c(16,1), col = c(alpha('black', alpha = 0.5),'red'), bty = 'n')
      # Print summed yields:
      fmt <- '%0.2e '
    }
    
    # Print summed yields:
    text(x=min(param$wInf)*sqrt(input$wMiddle/min(param$wInf)), 
         y=yl[1]*(yl[2]/yl[1])^0.95, adj=0.5,
         labels=sprintf(fmt, sum(Yield2[param$wInf<input$wMiddle])), col='red')
    text(x=input$wMiddle*sqrt(input$wLarge/input$wMiddle), 
         y=yl[1]*(yl[2]/yl[1])^0.95, adj=0.5,
         labels=sprintf(fmt, sum(Yield2[param$wInf>input$wMiddle & param$wInf<input$wLarge])), col='red')
    text(x=input$wLarge*sqrt(max(param$wInf)/input$wLarge), 
         y=yl[1]*(yl[2]/yl[1])^0.95, adj=0.5,
         labels=sprintf(fmt, sum(Yield2[param$wInf>=input$wLarge])), col='red')
    
    text(x=min(param$wInf)*sqrt(input$wMiddle/min(param$wInf)), 
         y=yl[1]*(yl[2]/yl[1])^0.85, adj=0.5,
         labels=sprintf(fmt, sum(Yield[param$wInf<input$wMiddle])), col='black')
    text(x=input$wMiddle*sqrt(input$wLarge/input$wMiddle), 
         y=yl[1]*(yl[2]/yl[1])^0.85, adj=0.5,
         labels=sprintf(fmt, sum(Yield[param$wInf>input$wMiddle & param$wInf<input$wLarge])), col='black')
    text(x=input$wLarge*sqrt(max(param$wInf)/input$wLarge), 
         y=yl[1]*(yl[2]/yl[1])^0.85, adj=0.5,
         labels=sprintf(fmt, sum(Yield[param$wInf>=input$wLarge])), col='black')
  })
  
  output$plotSSB <- renderPlot({
    SF <- simResults()[[1]]
    SF2 <- simResults()[[2]]      #param <- MM[[which(MM == 'param')]]
    param <- simResults()[[3]]
    
    SSB <- SF$SSB[length(SF$t),]
    SSB2 <- SF2$SSB[length(SF2$t),]
    # Sum biomass in small medium and large
    yl <- c(min(c(SSB,SSB2)),
            max(c(SSB,SSB2)))
    
    title <- 'Spawning stock biomass'
    
    if (input$Parameterset == 'Generic'){
      plot(param$wInf,SSB, log = 'xy', col = alpha('black',alpha = 0.5), type = 'l',
           xlab = 'Asymptotic weight (g)', 
           ylab = 'Spawning stock biomass',
           ylim = yl,
           xlim = c(min(param$wInf),max(param$wInf)),
           main = title, lwd = 3)
      lines(param$wInf,SSB2, col = alpha('red', alpha = 0.3), lwd = 3)
      lines(rep(input$wMiddle,100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      lines(rep(input$wLarge,100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      legend('bottomright', legend = c('Before', 'After'), lty = c(1,1), col = c(alpha('black', alpha = 0.5),'red'), bty = 'n')
    }else{
      plot(param$wInf,SSB, log = 'xy', col = alpha('black',alpha = 0.5),
           xlab = 'Asymptotic weight (g)', 
           ylab = 'Spawning stock biomass (ton)',
           ylim = yl,
           main = title, pch = 16, lwd = 3, cex = 2)
      points(param$wInf,SSB2, col = alpha('red', alpha = 0.5), cex = 2, lwd = 3)
      lines(rep(input$wMiddle,100), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      lines(rep(input$wLarge,wMiddle), seq(1e-15,yl[2]+1000, length.out = 100), lty = 2)
      legend('bottomleft', legend = c('Before', 'After'), pch = c(16,1), col = c(alpha('black', alpha = 0.5),'red'), bty = 'n')
    }
  })
  
  output$plotSpectrum <- renderPlot(expr={
    
    SF <- simResults()[[1]]
    SF2 <- simResults()[[2]]      #param <- MM[[which(MM == 'param')]]
    param <- simResults()[[3]]
    
    idx.N <- which(names(SF) == 'N')
    N <- SF[[idx.N]]
    
    idxEnd <- param$tEnd/param$dt
    Spectrum1 <- N[idxEnd,,]
    Spectrum1[Spectrum1 == 0] <- NA # For plotting on log scale 
    
    idx.N <- which(names(SF2) == 'N')
    N <- SF2[[idx.N]]
    Spectrum2 <- SF2$N[idxEnd,,]
    Spectrum2[Spectrum2 == 0] <- NA # For plotting
    
    w <- SF$w
    title <- 'Biomass density at equilibrium'
    
    #yl <- c(2*min(c(t(replicate(param$nSpecies, w))*Spectrum1,t(replicate(param$nSpecies, w))*Spectrum2), na.rm = T), 
    #        100*max(c(t(replicate(param$nSpecies, w))*Spectrum1,t(replicate(param$nSpecies, w))*Spectrum2), na.rm = T))
    yl <- param$kappaR * c(0.001*max(w)^(1+param$kR) , 0.1*0.02^(1+param$kR))
    
    plot(w,Spectrum1[1,]*w, log = 'xy', type = 'l', col = alpha('black',alpha = 0.3),
         ylab = 'Biomass density (-)', 
         xlab = 'Weight (g)', 
         main = title, 
         ylim=yl, xlim=c(0.02, max(w)))
    lines(w,Spectrum2[1,]*w, col = alpha('red',alpha = 0.3))
    
    for (i in 2:param$nSpecies){
      lines(w,Spectrum1[i,]*w, col = alpha('black',alpha = 0.3))
      lines(w,Spectrum2[i,]*w, col = alpha('red',alpha = 0.3))
    }
    #
    # Community spectrum:
    #
    lines(w, SF$Ntot[idxEnd,]*w, col=alpha('black',alpha = 0.3), lwd=3)
    lines(w, SF2$Ntot[idxEnd,]*w, col=alpha('red',alpha = 0.3), lwd=3)
    #
    # Resource
    #
    wPP <- SF$wPP
    lines(wPP, SF$nPP[idxEnd,1,]*wPP, col=alpha('black',alpha = 0.3), lwd=3, lty=2)
    lines(wPP, SF2$nPP[idxEnd,1,]*wPP, col=alpha('red',alpha = 0.3), lwd=3, lty=2)
  })
}


#' Run the community spectrum simulator
#'
#' @param ... Not used at the moment
#'
#' @return Nothing, it is used to run the Shiny app
#' @export
#'
#' @examples
#' runCommunitySpectrum()
runCommunitySpectrum <- function(...) {
  runApp(shinyApp(ui = shinyui, server = shinyserver))
}