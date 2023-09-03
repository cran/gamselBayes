########## R script: effTypesFromMFVB ##########

# For classification of a generalised additive
# model fit using mean field variational Bayes (MFVB)
# samples from the Bayesian gamsel-type model fit.

# Last changed: 03 AUG 2023

effTypesFromMFVB <- function(mu.q.gamma.beta,mu.q.gamma.u,lowerMakesSparser)
{
   # Determine dimension varibles:

   numPred <- length(mu.q.gamma.beta)
   if (is.null(mu.q.gamma.u))  dGeneral <- 0
   if (!is.null(mu.q.gamma.u)) dGeneral <- length(mu.q.gamma.u)
   dLinear <- numPred - dGeneral

   # Set up the effect type character string 
   # with "zero" starting values:

   effectTypesHat <- rep("zero",numPred)

   # Loop through the q*(beta_j) values and update
   # according to its coefficient being non-zero:

   for (iPred in 1:numPred)
   {
      if (mu.q.gamma.beta[iPred]>(1-lowerMakesSparser))
         effectTypesHat[iPred] <- "linear"
   }

   if (dGeneral>0)
   {
      iColSttPos <- 1

      for (jNon in 1:dGeneral)
      {
         if (mu.q.gamma.u[jNon]>(1-lowerMakesSparser))
            effectTypesHat[dLinear+jNon] <- "nonlinear"
      }
   }

   # Return the vector of estimated effect types:

   return(effectTypesHat)
}

############ End of effTypesFromMFVB ############
