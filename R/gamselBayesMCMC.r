########## R function: gamselBayesMCMC ##########

# For performing the Markov chain Monte Carlo
# iterations for fitting and inference for the 
# Bayesian gamsel-type model.

# Last changed: 31 JUL 2023

gamselBayesMCMC <- function(y,X,Z,ncZvec,family,XTy,XTX,ZTy,ZTX,ZTZ,
                            hyperPars,nWarm,nKept,nThin,msgCode)
{
   # Unpack hyperpameters:

   sigmaBeta0HYP <- hyperPars[1]
   sepsHYP <- hyperPars[2]
   sbetaHYP <- hyperPars[3]
   suHYP <- hyperPars[4]
   rhoBetaHYP <- hyperPars[5]
   rhoUHYP <- hyperPars[6]
 
   # Set MCMC dimension variable: 
  
   numMCMC <-  nWarm + nKept 

   # Set dimension values:

   X <- as.matrix(X)
   numPred <- ncol(X)
   dGeneral <- length(ncZvec)
   if (dGeneral>0) ncZmax <- max(ncZvec)
 
   # Obtain the "ZsttInds" and "ZendInds" vectors based on the
   # R array index convention:
   
   ZendIndsR <- cumsum(ncZvec)
   ZsttIndsR <- c(1,1+head(ZendIndsR,dGeneral-1))

   # Obtain the "ZsttInds" and "ZendInds" vectors base on the
   # C++ array index convention:

   ZendInds <- ZendIndsR - 1
   ZsttInds <- ZsttIndsR - 1

   # Allocate family number:

   if (family=="gaussian") familyNum <- 1
   if (family=="binomial") familyNum <- 2

   if (dGeneral==0)  # Set dummy values of Z-related quantities to avoid C++ problems:
   {
      Z <- matrix(0,1,1) ; ncZvec <- 1 ; ncZmax <- 1
      ZsttInds <- matrix(0,1,1) ; ZendInds <- matrix(0,1,1)
      ZTy <- matrix(0,1,1) ; ZTX <- matrix(0,1,1) 
      ZTZ <- matrix(0,1,1)
   }

   # Obtain MCMC samples:

   innerObj <- gamselBayesMCMCinner(y,X,Z,familyNum,ncZvec,ncZmax,dGeneral,
                                    ZsttInds,ZendInds,XTy,XTX,ZTy,ZTX, ZTZ,
                                    sigmaBeta0HYP,sepsHYP,sbetaHYP,suHYP,
                                    rhoBetaHYP,rhoUHYP,numMCMC,msgCode)

   # Extract samples and discard the warm-up values:

   beta0MCMC <- innerObj$beta0[-(1:nWarm)]
   betaTildeMCMC <- t(as.matrix(innerObj$betaTilde[,-(1:nWarm),drop=FALSE]))
   gammaBetaMCMC <- t(as.matrix(innerObj$gammaBeta[,-(1:nWarm),drop=FALSE]))
   sigmaBetaMCMC <- 1.0/sqrt(innerObj$recipSigsqBeta[-(1:nWarm)])
   sigmaEpsMCMC <- 1.0/sqrt(innerObj$recipSigsqEps[-(1:nWarm)])

   if (dGeneral>0)
   {
      gammaUMCMC <- t(as.matrix(innerObj$gammaU[,-(1:nWarm)]))
      uTildeArrayMCMC <- innerObj$uTilde[,,-(1:nWarm),drop=FALSE]
   }

   # Do thinning if 'nThin' exceeds unity:

   if (nThin>1)
   {
      thinnedInds <- seq(1,nKept,by=nThin)
      beta0MCMC <- beta0MCMC[thinnedInds]
      betaTildeMCMC <- as.matrix(betaTildeMCMC[thinnedInds,])
      gammaBetaMCMC <- as.matrix(gammaBetaMCMC[thinnedInds,])
      sigmaBetaMCMC <- sigmaBetaMCMC[thinnedInds]
      sigmaEpsMCMC <- sigmaEpsMCMC[thinnedInds]
      if (dGeneral>0)
      {
         gammaUMCMC <- as.matrix(gammaUMCMC[thinnedInds,])
         uTildeArrayMCMC <- uTildeArrayMCMC[,,thinnedInds,drop=FALSE]
      }
   }

   # Convert "sigmaEpsMCMC" to NULL for the Binomial family:

   if (family=="binomial")
      sigmaEpsMCMC <- NULL

   if (dGeneral>0)
   {
      # Organise the "uTilde" MCMC samples into a list
      # with each list entry corresponding to a particular
      # predictor and a matrix with rows corresponding to the
      # MCMC samples:

      uTildeMCMC <- vector("list",dGeneral)
      for (j in 1:dGeneral)
         uTildeMCMC[[j]] <- t(uTildeArrayMCMC[,j,])

   }

   if (dGeneral==0)
      uTildeMCMC <- NULL
  
   # Return kept samples:

    return(list(beta0=beta0MCMC,betaTilde=betaTildeMCMC,
                gammaBeta=gammaBetaMCMC,sigmaBeta=sigmaBetaMCMC,
                uTilde=uTildeMCMC,gammaU=gammaUMCMC,
                sigmaEps=sigmaEpsMCMC))
}

############ End of gamselBayesMCMC ############
