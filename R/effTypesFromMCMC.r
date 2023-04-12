########## R script: effTypesFromMCMC ##########

# For classification of a generalised additive
# model fit using Markov chain Monte Carlo (MCMC)
# samples from the Bayesian gamsel-type model fit.

# Last changed: 06 OCT 2021

effTypesFromMCMC <- function(betaMCMC,uMCMC,dLinear,lowerMakesSparser,effectiveZero)
{
   # Determine the number of predictors and the MCMC sample size:

   numPred <- ncol(betaMCMC)
   dGeneral <- numPred - dLinear
   nMCMC <- nrow(betaMCMC)   

   # Set up the effect type character string with "zero" 
   # starting values:

   effectTypesHat <- rep("zero",numPred)

   # Loop through the "beta" MCMC sample values
   # and update according to the spike at zero 
   # not being below a threshold value:

   for (iPred in 1:numPred)
   {
      betaCurr <- betaMCMC[,iPred]
      if (!coefZeroMCMC(betaCurr,lowerMakesSparser,effectiveZero))  
         effectTypesHat[iPred] <- "linear"
   }

   if (dGeneral>0)
   {
      iColSttPos <- 1

      # Loop through the remaining "u" MCMC sample values
      # and update according to *any* of the spikes at
      # zero not being below a threshold value:

      for (jNon in 1:dGeneral)
      {
         # Determine column numbers for the current predictor:

         for (iCol in  iColSttPos:ncol(uMCMC[[jNon]]))
         {
            ukCurr <- uMCMC[[jNon]][,iCol]
            if (!coefZeroMCMC(ukCurr,lowerMakesSparser,effectiveZero))  
               effectTypesHat[dLinear+jNon] <- "nonlinear"
         }
      }
   }

   # Return the vector of estimated effect types:

   return(effectTypesHat)
}

############ End of effTypesFromMCMC ############
