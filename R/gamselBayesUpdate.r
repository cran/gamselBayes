########## R function: gamselBayesUpdate ##########

# For updating a gamselBayes() fit object after possible adjustment
# to parameters controlling effect type estimation.

# Last changed: 03 AUG 2023

gamselBayesUpdate <- function(fitObject,lowerMakesSparser=NULL)
{
   # Check legality of non-NULL "lowerMakesSparser" values:

   if (!is.null(lowerMakesSparser))
   {
      if ((lowerMakesSparser<0)|(lowerMakesSparser>1))
      {
         warnStr1 <- "The inputted parameter for encouraging sparsity"
         warnStr2 <- "(lowerMakesSparser) is negative or exceeds 1."
         warnStr3 <- "The default value for the specified method was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n  ",warnStr3,"\n",sep=""),
                 immediate.=TRUE)
         if (fitObject$method=="MCMC") lowerMakesSparser <- 0.1
         if (fitObject$method=="MFVB") lowerMakesSparser <- 0
      }
   }

   # Obtain relevant sub-objects of the fit object:

   if (is.null(fitObject$Xlinear)) dLinear <- 0
   if (!is.null(fitObject$Xlinear))  dLinear <- ncol(fitObject$Xlinear)
   if (is.null(fitObject$Xgeneral)) dGeneral <- 0
   if (!is.null(fitObject$Xgeneral)) dGeneral <- ncol(fitObject$Xgeneral)
   method <- fitObject$method
   MCMCobj <- fitObject$MCMC
   MFVBobj <- fitObject$MFVB

   # Estimate the effect type:
 
   if (method=="MCMC")
   {
      # Extract the MCMC samples for the coefficients:
      
      betaTildeMCMC <- MCMCobj$betaTilde
      gammaBetaMCMC <- MCMCobj$gammaBeta
      uTildeMCMC <- MCMCobj$uTilde
      gammaUMCMC <- MCMCobj$gammaU
      betaMCMC <- gammaBetaMCMC*betaTildeMCMC

       
      if (dGeneral>0)
      {
         uMCMC <- vector("list",dGeneral)
         for (j in 1:dGeneral)
            uMCMC[[j]] <- gammaUMCMC[,j]*uTildeMCMC[[j]]
      }
      effectTypesHat <- effTypesFromMCMC(gammaBetaMCMC,gammaUMCMC,lowerMakesSparser)
   }

   if (method=="MFVB")
   {
      # Extract the MFVB parameters for the coefficients:

      mu.q.betaTilde <- MFVBobj$betaTilde$mu.q.betaTilde
      sigsq.q.betaTilde <- diag(MFVBobj$betaTilde$Sigma.q.betaTilde)
      mu.q.gamma.beta <- MFVBobj$gammaBeta
      mu.q.uTilde <- MFVBobj$uTilde$mu.q.uTilde
      sigsq.q.uTilde <- MFVBobj$uTilde$sigsq.q.uTilde
      mu.q.gamma.u <- MFVBobj$gammaU

      effectTypesHat <- effTypesFromMFVB(mu.q.gamma.beta,mu.q.gamma.u,lowerMakesSparser)
   }

   # Update the "effectTypesHat" component of the fit object:

   outObj <- fitObject
   outObj$effectTypesHat <- effectTypesHat
   
   # Return the updated gamseBayes object:

   class(outObj) <- "gamselBayes"

   return(outObj)
}   

############ End gamselBayesUpdate ###########

