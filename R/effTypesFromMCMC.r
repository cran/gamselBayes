########## R script: effTypesFromMCMC ##########

# For classification of a generalised additive
# model fit using Markov chain Monte Carlo (MCMC)
# samples from the Bayesian gamsel-type model fit.

# Last changed: 03 AUG 2023

effTypesFromMCMC <- function(gammaBetaMCMC,gammaUMCMC,lowerMakesSparser)
{
   # Determine the number of predictors and the MCMC sample size:

   numPred <- ncol(gammaBetaMCMC)
   if (is.null(gammaUMCMC)) dGeneral <- 0
   if (!is.null(gammaUMCMC)) dGeneral <- ncol(gammaUMCMC)
   dLinear <- numPred - dGeneral
   nMCMC <- nrow(gammaBetaMCMC)   

   # Set up the effect type character string with "zero" 
   # starting values:

   effectTypesHat <- rep("zero",numPred)

   # Loop through the "beta" MCMC sample values
   # and update according to the spike at zero 
   # not being below a threshold value:

   for (iPred in 1:numPred)
   {
      if (mean(gammaBetaMCMC[,iPred])>(1-lowerMakesSparser))          
         effectTypesHat[iPred] <- "linear" 
   }

   if (dGeneral>0)
   {
      for (jNon in 1:dGeneral)
      {
         if (mean(gammaUMCMC[,jNon])>(1-lowerMakesSparser))         
            effectTypesHat[dLinear+jNon] <- "nonlinear"    
      }
   }

   # Return the vector of estimated effect types:

   return(effectTypesHat)
}

############ End of effTypesFromMCMC ############
