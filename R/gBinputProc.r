########## R function: gBinputProc ##########

# For conducting input processing for the gamselBayes() function:

# Last changed: 07 AUG 2023

gBinputProc <- function(y,family,method,lowerMakesSparser)
{
   # Check legality of "family":

   if (!any(family==c("binomial","gaussian")))
   {
      warnStr1 <- "The inputted family is not one of the available options."
      warnStr2 <- "The \"gaussian\" default method was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      family <- "gaussian"
   }

   # Check legality of "method":

   if (!any(method==c("MCMC","MFVB")))
   {
      warnStr1 <- "The inputted method is not one of the available options."
      warnStr2 <- "The \"MCMC\" default method was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      method <- "MCMC"
   }

   # Check the legality of "y" for the family="binomial":

   if (family=="binomial")
   {
      if (!vectorIsBinary(y))
         stop("All entries of y must be 0 or 1 when family is \"binomial\".\n")
   }  

   # Set the default values of the "lowerMakesSparser" parameter:

   if (is.null(lowerMakesSparser))
   {
      if (method=="MCMC") lowerMakesSparser <- 0.5
      if (method=="MFVB") lowerMakesSparser <- 0.1
   }

   # Make sure that an inputted "lowerMakesSparser" is legal:

   if (!is.na(lowerMakesSparser))
   {
      if ((lowerMakesSparser<0)|(lowerMakesSparser>1))
      {
         warnStr1 <- "The inputted parameter for encouraging sparsity"
         warnStr2 <- "(lowerMakesSparser) is negative or exceeds 1."
         warnStr3 <- "The default value for the specified method was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n  ",warnStr3,"\n",sep=""),
                 immediate.=TRUE)
         if (method=="MCMC") lowerMakesSparser <- 0.5
         if (method=="MFVB") lowerMakesSparser <- 0.1
      }
   }

   # Return the processed inputs:

   return(list(family=family,method=method,lowerMakesSparser=lowerMakesSparser))
}

############ End of gBinputProc ############
