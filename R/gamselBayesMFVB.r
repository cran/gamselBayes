########## R function: gamselBayesMFVB ##########

# For performing the mean field variational Bayes
# iterations for fitting and inference for the 
# Bayesian gamsel-type model.

# Last changed: 22 AUG 2023

gamselBayesMFVB <- function(y,X,Z,ncZvec,family,XTy,XTX,ZTy,ZTX,ZTZ,hyperPars,
                            maxIter,toler,msgCode)
{
   # Unpack hyperpameters:

   sigmaBeta0HYP <- hyperPars[1]
   sepsHYP <- hyperPars[2]
   sbetaHYP <- hyperPars[3]
   suHYP <- hyperPars[4]
   rhoBetaHYP <- hyperPars[5]
   rhoUHYP <- hyperPars[6]

   # Set dimension  variables:

   n <- length(y)
   ncX <- ncol(X)
   dGeneral <- length(ncZvec)
   if (dGeneral>0) ncZmax <- max(ncZvec)

   # Obtain the "ZsttInds" and "ZendInds" vectors:

   ZendInds <- cumsum(ncZvec)
   ZsttInds <- c(1,1+head(ZendInds,dGeneral-1))

   # Make an adjustment for the C++ array index convention:

   ZendInds <- ZendInds - 1
   ZsttInds <- ZsttInds - 1

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
 
   # Carry out mean field variational Bayes iterations:

   innerObj <- gamselBayesMFVBinner(y,X,Z,familyNum,ncZvec,ncZmax,dGeneral,
                                    ZsttInds,ZendInds,XTy,XTX,ZTy,ZTX,ZTZ,
                                    sigmaBeta0HYP,sepsHYP,sbetaHYP,suHYP,
                                    rhoBetaHYP,rhoUHYP,maxIter,toler,msgCode) 

   # Extract q-density parameters:

   mu.q.beta.0 <- innerObj$muqBetaZero
   sigsq.q.beta.0  <- innerObj$sigsqqBetaZero
   mu.q.betaTilde <- as.vector(innerObj$muqBetaTilde)
   Sigma.q.betaTilde <- innerObj$SigmaqBetaTilde
   mu.q.gamma.beta <- as.vector(innerObj$muqgammaBeta)
   kappa.q.sigsq.beta <- innerObj$kappaqSigsqBeta 
   lambda.q.sigsq.beta <- innerObj$lambdaqSigsqBeta 
   kappa.q.sigsq.eps <- innerObj$kappaqSigsqEps 
   lambda.q.sigsq.eps <- innerObj$lambdaqSigsqEps 
   logMargLik <- as.vector(innerObj$logMargLik)
   numIters <- innerObj$numIters
   relErr <- innerObj$relErr
   stopType <- innerObj$stopType

   # If the MFVB iterations were stopped, due to convergence being attained, before
   # the maximum number of iterations was reached then print a message about this:

   if ((msgCode>0)&(stopType==0))
   {
      msgStr1 <- "\n   The mean field variational Bayes algorithm converged"
      msgStr2 <- paste(" after ",numIters," iterations.\n\n",sep="")
      cat(paste(msgStr1,msgStr2,sep=""))
   }

   # If the MFVB iterations were stopped, due to "maxIter" being reached, before
   # convergence was achieved then print a warning message:

   if ((msgCode!=0)&(stopType==1))
   {
      msgStr1 <- "\n   NOTE: The mean field variational Bayes iterations did not attain"
      msgStr2 <- "   the specified convergence tolerance before the maximum number of"
      msgStr3 <- paste("   iterations was reached. The relative error is ",signif(relErr,4)," which",sep="")
      msgStr4 <- paste("   exceeds the tolerance of ",signif(toler,4),".",sep="")
      cat(paste(msgStr1,"\n",msgStr2,"\n",msgStr3,"\n",msgStr4,"\n",sep=""))
   }

   # Truncate the logMargLik vector to match the number of iterations
   # before convergence:

   if ((stopType==0)&(length(logMargLik)>=(numIters+1)))  
      logMargLik <- logMargLik[-((numIters+1):length(logMargLik))]   
   
   if (dGeneral>0)
   {
      mu.q.gamma.u <- as.vector(innerObj$muqgammaU)
      muqUtilde <- innerObj$muqUtilde
      sigsqqUtilde <- innerObj$sigsqqUtilde
      mu.q.uTilde <- vector("list",dGeneral)
      sigsq.q.uTilde <- vector("list",dGeneral)
      for (j in 1:dGeneral)
      {
         mu.q.uTilde[[j]] <- muqUtilde[,j]
         sigsq.q.uTilde[[j]] <- sigsqqUtilde[,j]
      }
   }

   if (dGeneral==0)
   {
      sigsq.q.uTilde <- NULL
      mu.q.gamma.u <- NULL
      mu.q.uTilde <- NULL
   }

   return(list(beta0=c(mu.q.beta.0=mu.q.beta.0,sigsq.q.beta.0=sigsq.q.beta.0),
               betaTilde=list(mu.q.betaTilde=mu.q.betaTilde,Sigma.q.betaTilde=Sigma.q.betaTilde),
               gammaBeta=mu.q.gamma.beta,sigmaBeta=c(kappa.q.sigsq.beta,lambda.q.sigsq.beta),
               uTilde=list(mu.q.uTilde=mu.q.uTilde,sigsq.q.uTilde=sigsq.q.uTilde),
               gammaU=mu.q.gamma.u,sigmaEps=c(kappa.q.sigsq.eps,lambda.q.sigsq.eps),
               logMargLik=logMargLik))
}

############ End of gamselBayesMFVB ############
