########## R-function: coefZeroMFVB ##########

# Decides whether a coefficient is zero or not
# based on its Normal-Zero q-density.

# Last changed: 25 OCT 2021

coefZeroMFVB <- function(probZero,lowerMakesSparser,effectiveZero)
{
   maxVal <- max(c(lowerMakesSparser,effectiveZero))
   return(probZero>maxVal)
}

########## End of coefZeroMFVB ############

