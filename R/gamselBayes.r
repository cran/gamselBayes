########## R function: gamselBayes ##########

# For performing generalized additive model selection via a 
# Bayesian inference engine approach.

# Last changed: 03 AUG 2023

gamselBayes <- function(y,Xlinear=NULL,Xgeneral=NULL,method="MCMC",lowerMakesSparser=NULL,
                        family="gaussian",verbose=TRUE,control=gamselBayes.control())
{
   # Process some of the inputs, mainly checking their legalities:

   processedInput <- gBinputProc(y,family,method,lowerMakesSparser)

   family <- processedInput$family
   method <- processedInput$method
   lowerMakesSparser <- processedInput$lowerMakesSparser

   # Process some of the inputs, mainly checking their legalities:

   processedData <- gBdataPreProc(y,Xlinear,Xgeneral,family)
  
   y <- processedData$y
   Xlinear <- processedData$Xlinear
   Xgeneral <- processedData$Xgeneral
   X <- processedData$X
   dLinear <- processedData$dLinear
   dGeneral <- processedData$dGeneral
   meany <- processedData$meany
   sdy <- processedData$sdy
   meanXlinear <-  processedData$meanXlinear
   sdXlinear <-  processedData$sdXlinear
   meanXgeneral <-  processedData$meanXgeneral
   sdXgeneral <-  processedData$sdXgeneral

   # Extract the sample size:

   sampSize <- length(y)

   # Pass the verbose specification to the msgCode variable:

   msgCodeVerbose <- as.numeric(verbose)

   # Unpack the control values:

   numIntKnots  <- control$numIntKnots
   truncateBasis <- control$truncateBasis
   numBasis <- control$numBasis
   sigmaBeta0HYP  <- control$sigmabeta0
   sbetaHYP  <- control$sbeta
   sepsHYP  <- control$sepsilon
   suHYP  <- control$su
   rhoBetaHYP <- control$rhoBeta
   rhoUHYP <- control$rhoU
   nWarm <- control$nWarm
   nKept <- control$nKept
   nThin <- control$nThin
   maxIter <- control$maxIter
   toler <- control$toler
   msgCodeControl <- control$msgCode

   # Sort out conflicts between the "verbose" specification and the 
   # "msgCode" specification.

   if ((msgCodeVerbose==1)&(msgCodeControl==0))
   {              
      warnStr1 <- "The verbose and msgCode specifications conflict with"
      warnStr2 <- "each other. The default value of msgCode=1 was used."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      msgCodeControl <- 1 
   }

   # Determine the value of "msgCode" based on current values of 
   # "msgCodeVerbose" and "msgCodeControl".

   if (msgCodeVerbose==0) msgCode <- 0
   if (msgCodeVerbose>0) msgCode <- msgCodeControl
   
   # Obtain Z matrix containing spline basis functions 
   # and auxiliary quantities: 
 
   ZmatObj <- gamselBayesZproc(X,sampSize,dLinear,dGeneral,numIntKnots,
                               truncateBasis,numBasis)
   Z <- ZmatObj$Z
   ncZvec <- ZmatObj$ncZvec
   rangexList <- ZmatObj$rangexList
   intKnotsList <- ZmatObj$intKnotsList
   OStoDRmatList <- ZmatObj$OStoDRmatList

   # Obtain sufficient statistic quantities:

   suffStatObj <- gBsuffStats(y,X,Z)
   XTy <- suffStatObj$XTy
   XTX <- suffStatObj$XTX
   ZTy <- suffStatObj$ZTy
   ZTX <- suffStatObj$ZTX
   ZTZ <- suffStatObj$ZTZ
 
   # Organise the hyperparameters:
 
   hyperPars <- c(sigmaBeta0HYP,sepsHYP,sbetaHYP,suHYP,rhoBetaHYP,
                  rhoUHYP)

   # Obtain fits:

   if (msgCode>0) cat("\n")
   if (method=="MCMC")
   {
      MCMCobj <- gamselBayesMCMC(y,X,Z,ncZvec,family,XTy,XTX,ZTy,ZTX,ZTZ,
                                 hyperPars,nWarm,nKept,nThin,msgCode)
      MFVBobj <- NULL 
   }

   if (method=="MFVB")
   {
      MFVBobj <- gamselBayesMFVB(y,X,Z,ncZvec,family,XTy,XTX,ZTy,ZTX,ZTZ,
                                 hyperPars,maxIter,toler,msgCode)
      MCMCobj <- NULL
   }
   if (msgCode>0) cat("\n")

   # Estimate the effect type:
 
   if (method=="MCMC")
   {
      # Extract the MCMC samples for the "gamma" variables:
      
      gammaBetaMCMC <- MCMCobj$gammaBeta
      gammaUMCMC <- MCMCobj$gammaU

      effectTypesHat <- effTypesFromMCMC(gammaBetaMCMC,gammaUMCMC,lowerMakesSparser)
   }    

   if (method=="MFVB")
   {
      # Extract the MFVB parameters for the "gamma" variables"

      mu.q.gamma.beta <- MFVBobj$gammaBeta
      mu.q.gamma.u <- MFVBobj$gammaU

      effectTypesHat <- effTypesFromMFVB(mu.q.gamma.beta,mu.q.gamma.u,lowerMakesSparser)
   }

   # Prepare and return the output object:

   outObj <- list(method=method,family=family,Xlinear=Xlinear,Xgeneral=Xgeneral,rangex=rangexList,
                  intKnots=intKnotsList,OStoDRmat=OStoDRmatList,truncateBasis=truncateBasis,
                  numBasis=numBasis,MCMC=MCMCobj,MFVB=MFVBobj,
                  effectTypesHat=effectTypesHat,meany=meany,sdy=sdy,meanXlinear=meanXlinear,
                  sdXlinear=sdXlinear,meanXgeneral=meanXgeneral,sdXgeneral=sdXgeneral)

   class(outObj) <- "gamselBayes"

   return(outObj)
}   

############ End gamselBayes ###########

