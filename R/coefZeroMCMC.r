########## R-function: coefZeroMCMC ##########

# Decides whether a coefficient is zero or not
# based on its MCMC samples.

# Last changed: 25 OCT 2021

coefZeroMCMC <- function(coefMCMC,lowerMakesSparser,effectiveZero)
{
   numZeroes <- sum(as.numeric(abs(coefMCMC)<sqrt(.Machine$double.eps)))
   propZeroes <- numZeroes/length(coefMCMC)
   maxVal <- max(c(lowerMakesSparser,effectiveZero))
   return(propZeroes>maxVal)
}

########## End of coefZeroMCMC ############

