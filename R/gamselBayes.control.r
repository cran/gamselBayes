########## R function: gamselBayes.control ##########

# Control function for gamselBayes().

# Last changed: 03 AUG 2023

gamselBayes.control <- function(numIntKnots=25,truncateBasis=TRUE,numBasis=12,
                                sigmabeta0=100000,sbeta=1000,sepsilon=1000,su=1000,rhoBeta=0.5,rhoU=0.5,
                                nWarm=1000,nKept=1000,nThin=1,maxIter=1000,toler=1e-8,msgCode=1)
{
   # Make sure that numIntKnots is legal:

   numIntKnots <- round(numIntKnots)
   if ((numIntKnots<7)|(numIntKnots>50))
   {
      warnStr1 <- "The inputted number of interior knots is not between 8."
      warnStr2 <- "and 50. The default value of 25 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      numIntKnots <- 25
   }

   # Make sure that truncateBasis is legal:

   if (!(any(truncateBasis==c(FALSE,TRUE))))
   {
      warnStr1 <- "The inputted flag for truncation of the spline basis is neither"
      warnStr2 <- "FALSE nor TRUE. The default value of TRUE was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      truncateBasis <- TRUE
   }

   # Make sure that numBasis is not too low:

   numBasis <- round(numBasis)
   if (numBasis<5)
   {
      warnStr1 <- "The inputted number of basis functions is negative or too low."
      warnStr2 <- "The default value of 12 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      numBasis <- 12
   }

   # Make sure that numBasis functions does not exceed (numIntKnots+2):
  
   if (numBasis>(numIntKnots+2))
   {
      warnStr1 <- "The inputted number of basis functions is too high compared with"
      warnStr2 <- "the value of numIntKnots. The default value of 12 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      numBasis <- 12
   }

   # Make sure that sigmabeta0 is legal:
  
   if (sigmabeta0<=0)
   {
      warnStr1 <- "The inputted sigma_beta0 standard deviation hyperparameter is non-positive."
      warnStr2 <- "The default value of 100000 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      sigmabeta0 <- 100000
   }

   # Make sure that sbeta is legal:
  
   if (sbeta<0)
   {
      warnStr1 <- "The inputted sigma_beta scale hyperparameter is non-positive."
      warnStr2 <- "The default value of 1000 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      sbeta <- 1000
   }

   # Make sure that sepsilon is legal:
  
   if (sepsilon<0)
   {
      warnStr1 <- "The inputted sigma_epsilon scale hyperparameter is non-positive."
      warnStr2 <- "The default value of 1000 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      sepsilon <- 1000
   }

   # Make sure that su is legal:
  
   if (su<0)
   {
      warnStr1 <- "The inputted sigma_u scale hyperparameter is non-positive."
      warnStr2 <- "The default value of 1000 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      su <- 1000
   }

   # Make sure that rhoBeta is legal:
  
   if ((rhoBeta<=0)|(rhoBeta>=1))
   {
      warnStr1 <- "The inputted rho_beta probability hyperparameter is outside (0,1)."
      warnStr2 <- "The default value of 0.5 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      rhoBeta <- 0.5
   }

   # Make sure that rhoU is legal:
  
   if ((rhoU<=0)|(rhoU>=1))
   {
      warnStr1 <- "The inputted rho_U probability hyperparameter is outside (0,1)."
      warnStr2 <- "The default value of 0.5 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      rhoU <- 0.5
   }

   if (!is.null(nWarm))
   {
      # Make sure that nWarm is legal:

      nWarm <- round(nWarm)
      if (nWarm<0) 
      {
         warnStr1 <- "The inputted number of warmup Markov chain Monte Carlo iterations is negative."
         warnStr2 <- "The default value of 1000 was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
         nWarm <- 1000
      }
   }

   # Make sure that nKept is legal:

   if (!is.null(nKept))
   {
      nKept <- round(nKept)
      if (nKept<0) 
      {
         warnStr1 <- "The inputted number of kept Markov chain Monte Carlo iterations is negative."
         warnStr2 <- "The default value of 1000 was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
         nKept <- 1000
      }
   }

   # Make sure that nThin is legal:

   nThin <- round(nThin)
   if (nThin<0) 
   {
      warnStr1 <- "The inputted Markov chain Monte Carlo thinning factor is negative."
      warnStr2 <- "The default value of 1 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      nThin <- 1
   }
   
   if (!is.null(nKept))
   {
      if (nThin>nKept) 
      {
         warnStr1 <- "The inputted Markov chain Monte Carlo thinning factor exceeds the"
         warnStr2 <- "number of kept values. The default value of 1 was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
         nThin <- 1
      }
   }

   if (!is.null(maxIter))
   {
      # Make sure that maxIter is legal:

      maxIter <- round(maxIter)
      if (maxIter<0) 
      {
         warnStr1 <- "The inputted maximum number of mean field variational Bayes iterations is negative."
         warnStr2 <- "The default value of 1000 was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
         maxIter <- 1000
      }
   }


   if (!is.null(toler))
   {
      # Make sure that toler is legal:

      if ((toler<0)|(toler>0.5)) 
      {
         warnStr1 <- "The inputted mean field variational Bayes convergence tolerance is negative."
         warnStr2 <- "The default value of 0.00001 was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
         toler <- 0.00001
      }
   }

   # Make sure that msgCode is legal:

   msgCode <- round(msgCode)
   if (!any(msgCode==(0:2)))
   {
      warnStr1 <- "The inputted message code number is not 0, 1 or 2."
      warnStr2 <- "The default value of 1 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      msgCode <- 1
   }

   return(list(numIntKnots=numIntKnots,truncateBasis=truncateBasis,numBasis=numBasis,
               sigmabeta0=sigmabeta0,sbeta=sbeta,sepsilon=sepsilon,
               su=su,rhoBeta=rhoBeta,rhoU=rhoU,nWarm=nWarm,nKept=nKept,nThin=nThin,
               maxIter=maxIter,toler=toler,msgCode=msgCode))
}

############ End of gamselBayes.control ############
