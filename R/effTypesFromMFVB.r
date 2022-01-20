########## R script: effTypesFromMFVB ##########

# For classification of a generalised additive
# model fit using mean field variational Bayes (MFVB)
# samples from the Bayesian gamsel-type model fit.

# Last changed: 06 OCT 2021

effTypesFromMFVB <- function(mu.q.betaTilde,sigsq.q.betaTilde,mu.q.gamma.beta,
                            mu.q.uTilde,sigsq.q.uTilde,mu.q.gamma.u,
                            lowerMakesSparser,effectiveZero)
{
   # Determine dimension varibles:

   numPred <- length(mu.q.gamma.beta)
   dGeneral <- length(mu.q.gamma.u)
   dLinear <- numPred - dGeneral

   # Set up the effect type character string 
   # with "zero" starting values:

   effectTypesHat <- rep("zero",numPred)

   # Loop through the q*(beta_j) values and update
   # according to its coefficient being non-zero:

   for (iPred in 1:numPred)
   {
      isBetaZero <- coefZeroMFVB((1-mu.q.gamma.beta[iPred]),lowerMakesSparser,effectiveZero)

      if (!isBetaZero)
         effectTypesHat[iPred] <- "linear"
   }

   if (dGeneral>0)
   {
      iColSttPos <- 1

      # Loop through the remaining "u" mu.q.gamma.u
      # and update according to *any* of the spikes at
      # zero not being below a threshold value:

      for (jNon in 1:dGeneral)
      {
         for (iCol in iColSttPos:length(mu.q.gamma.u[[jNon]]))
         {
            isUkZero <- coefZeroMFVB((1-mu.q.gamma.u[[jNon]][iCol]),lowerMakesSparser,effectiveZero)
            if (!isUkZero)
               effectTypesHat[dLinear+jNon] <- "nonlinear"
         }
      }
   }

   # Return the vector of estimated effect types:

   return(effectTypesHat)
}

############ End of effTypesFromMFVB ############
