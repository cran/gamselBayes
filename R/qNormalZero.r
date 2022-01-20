########## R-function: qNormalZero ##########

# Draws a sample from a Normal-Zero distribution.

# Last changed: 15 JAN 2021

qNormalZero <- function(p,mu=0,sigma=1,probNonZero=0.5)
{
   # Determine change points:

   CPlow <- probNonZero*pnorm(-mu/sigma)
   CPupp <- CPlow + 1 - probNonZero

   # Determine lower and upper functions:

   argLow <- p/probNonZero
   qnormLow <- rep(0,length(argLow))
   legalInds <- (1:length(argLow))[(argLow>0)&(argLow<1)]
   qnormLow[legalInds] <- qnorm(argLow[legalInds])
   lowFun <- (mu + sigma*qnormLow)*as.numeric(p<CPlow)

   argUpp <- (p+probNonZero-1)/probNonZero
   qnormUpp <- rep(0,length(argUpp))
   legalInds <- (1:length(argLow))[(argUpp>0)&(argUpp<1)]
   qnormUpp[legalInds] <- qnorm(argUpp[legalInds])
   uppFun <- (mu + sigma*qnormUpp)*as.numeric(p>=CPupp)

   # Return the sum of the lower and upper functions:
    
   return(lowFun + uppFun)
}

########## End of qNormalZero ############

