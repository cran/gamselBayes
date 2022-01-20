########## R-function: rNormalZero ##########

# For drawing a random sample from a Normal-Zero 
# distribution.

# Last changed: 15 JAN 2021

rNormalZero <- function(n,mu=0,sigma=1,probNonZero=0.5)
{
   # Draw the auxiliary mixing variables:

   auxMix <- rbinom(n,1,probNonZero)
   nNonZero <- sum(auxMix)

   # Obtain and return Normal-Zero draws:

   x <- rep(0,n)
   if (nNonZero>0)
      x[auxMix==1] <- rnorm(nNonZero,mu,sigma)

   return(x)   
}

########## End of rNormalZero ############

