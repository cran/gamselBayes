########## R function: summGBargsProc ##########

# For conducting argument processing for the 
# summary.gamselBayes() function:

# Last changed: 22 JUL 2021

summGBargsProc <- function(credLev,sigFigs,nMC)
{
   # Check the legality of the inputted "credLev" value:

   if (credLev<=0) 
   {
      warnStr1 <- "The inputted credible level is zero or negative."
      warnStr2 <- "The default value of 0.95 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      credLev <- 0.95
   }
   if (credLev>=1) 
   {
      warnStr1 <- "The inputted credible level is 1 or higher."
      warnStr2 <- "The default value of 0.95 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      credLev <- 0.95
   }

   # Check the legality of the inputted "sigFigs" value:

   sigFigs <- round(sigFigs)
   if (sigFigs<0) 
   {
      warnStr1 <- "The inputted number of significant figures is zero or negative."
      warnStr2 <- "The default value of 5 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      sigFigs <- 5
   }

   # Deal with the case where the inputted "sigFigs" exceeds 7:

   if (sigFigs>7) 
   {
      warnStr1 <- "The inputted number of significant figures exceeds 7. The"
      warnStr2 <- "following command may be required to achieve desired result:"
      warnStr3 <- paste("options(digits=",sigFigs,")",sep="")
      warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
              immediate.=TRUE)
   }

   # Check the legality of the inputted "nMC" value:

   nMC <- round(nMC)
   if (nMC<51)
   {
      warnStr1 <- "The inputted number of Monte Carlo samples for variational."
      warnStr2 <- "inference is too low. The default value of 10000 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      nMC <- 10000
   }

   # Return the processed inputs:

   return(list(credLev=credLev,sigFigs=sigFigs,nMC=nMC))
}

############ summGBargsProc ############
