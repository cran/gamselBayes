########## R function: plot.gamselBayes ##########

# For plotting a gamselBayes() fit object.

# Last changed: 31 JUL 2023

plot.gamselBayes <- function(x,credLev=0.95,gridSize=251,nMC=5000,varBand=TRUE,shade=TRUE,
                             yscale="response",rug=TRUE,rugSampSize=NULL,estCol="darkgreen",
                             varBandCol=NULL,rugCol="dodgerblue",mfrow=NULL,xlim=NULL,ylim=NULL,
                             xlab=NULL,ylab=NULL,mai=NULL,pages=NULL,cex.axis=1.5,cex.lab=1.5,...)
{
   fitObject <- x

   # Extract the sample size value:

   if (!is.null(fitObject$Xlinear))  sampSizeVal <- nrow(fitObject$Xlinear)
   if (!is.null(fitObject$Xgeneral)) sampSizeVal <- nrow(fitObject$Xgeneral)

   # Process the other inputs, with the exception of "xlim" and "ylim":

   processedInputs <- plotGBinputProc(credLev,gridSize,nMC,varBand,yscale,shade,rug,
                                      rugSampSize,estCol,varBandCol,rugCol,mai,pages,
                                      sampSizeVal)
   credLev <- processedInputs$credLev
   gridSize <- processedInputs$gridSize
   nMC <- processedInputs$nMC
   varBand <- processedInputs$varBand
   shade <- processedInputs$shade
   yscale <- processedInputs$yscale
   rug <- processedInputs$rug
   rugSampSize <- processedInputs$rugSampSize
   estCol <- processedInputs$estCol
   varBandCol <- processedInputs$varBandCol
   rugCol <- processedInputs$rugCol
   mai <- processedInputs$mai
   pages <- processedInputs$pages

   # Obtain effect type and parameter and estimates depending
   # on the method type:

   effectTypesHat <- fitObject$effectTypesHat

   # Determine the indices of the predictors which have a 
   # non-linear effect estimate:

   indsEstNonlinEff <- (1:length(effectTypesHat))[effectTypesHat=="nonlinear"]

   # Set the number of panels, which corresponds to the number of 
   # predictors which are estimated to have a non-linear effect:

   numPanels <- length(indsEstNonlinEff)

   # Check the legality of the inputted "xlim" and "ylim" values, starting with
   # the number of panels to be plotted:

   processedLimits <- plotGBlimsProc(xlim,ylim,numPanels)
   xlimMat <- processedLimits$xlim
   ylimMat <- processedLimits$ylim
  
   if (numPanels==0)
   {
      message(paste(" None of the predictors have estimated non-linear effects\n",
                   "and so there is nothing to plot.\n"))
   }
   if (numPanels>0)
   {
      # Extract the method:

      method <- fitObject$method

      # Extract the predictor values:

      Xlinear <- fitObject$Xlinear
      Xgeneral <- fitObject$Xgeneral

      # Obtain the combined predictor matrix:

      if (is.null(Xlinear)) X <- as.matrix(Xgeneral)
      if (is.null(Xgeneral)) X <- as.matrix(Xlinear)
      if ((!is.null(Xlinear))&(!is.null(Xgeneral))) 
         X <- as.matrix(cbind(Xlinear,Xgeneral))

      # Extract the "dLinear", "dGeneral" values:

      if (is.null(Xlinear))     dLinear <- 0
      if (!is.null(Xlinear))    dLinear <- ncol(Xlinear)
      if (is.null(Xgeneral))  dGeneral <- 0
      if (!is.null(Xgeneral)) dGeneral <- ncol(Xgeneral)
   
      # Obtain the vectors of means and standard deviations required to 
      # convert axes to the orginal units:

      meany <- fitObject$meany
      sdy <- fitObject$sdy
      if (dLinear>0)
      {
         meanXall <- c(fitObject$meanXlinear,fitObject$meanXgeneral)
         sdXall <- c(fitObject$sdXlinear,fitObject$sdXgeneral)         
      } 
      if (dLinear==0)
      {
         meanXall <- fitObject$meanXgeneral
         sdXall <- fitObject$sdXgeneral         
      } 

      meanX <- meanXall[indsEstNonlinEff]
      sdX <- sdXall[indsEstNonlinEff]

      # Determine the vector of x-axis labels for each panel:

      if (is.null(xlab))
      {
         if (dLinear==0) xLabsVecFull <- names(Xgeneral)
         if (dGeneral==0) xLabsVecFull <- names(Xlinear)
         if ((dLinear>0)&(dGeneral>0)) xLabsVecFull <- c(names(Xlinear),names(Xgeneral))
         xLabsVec <- xLabsVecFull[indsEstNonlinEff]
      } 
      if (!is.null(xlab))
      {
         xlab <- as.character(xlab)
         if (length(xlab)!=numPanels) 
         {
            stopStr1 <- "xlab must be a vector with number of entries equal to the number of panels."
            stopStr2 <- paste("The number of panels for this fit object is ",numPanels,".",sep="")
            stop(paste(stopStr1,"\n",stopStr2,"\n",sep=""))
         }
         xLabsVec <- xlab
      }

      # Determine the vector of y-axis labels for each panel:

      if (is.null(ylab))
      {
         # Determine the default y-axis label: 
        
         if ((fitObject$family=="gaussian")&(yscale=="response"))
            ylabDefault <- "mean response"
     
         if ((fitObject$family=="binomial")&(yscale=="response"))
            ylabDefault <- "probability"

         if ((fitObject$family=="binomial")&(yscale=="link"))
            ylabDefault <- "log odds"

         yLabsVec <- rep(ylabDefault,numPanels)

      } 

      if (!is.null(ylab))
      {
         ylab <- as.character(ylab)
         if (length(ylab)!=numPanels) 
         {
            stopStr1 <- "ylab must be a vector with number of entries equal to the number of panels."
            stopStr2 <- paste("The number of panels for this fit object is ",numPanels,".",sep="")
            stop(paste(stopStr1,"\n",stopStr2,"\n",sep=""))
         }
         yLabsVec <- ylab
      }

      # Subset the columns of X to correspond to those predictors for 
      # which there is an estimated non-linear effect:

      XestNonlin <- as.matrix(X[,indsEstNonlinEff])

      # Extract spline basis function parameters:

      rangexList <- fitObject$rangex
      intKnotsList <- fitObject$intKnots
      OStoDRmatList <- fitObject$OStoDRmat
      truncateBasis <- fitObject$truncateBasis
      numBasis <- fitObject$numBasis
 
      if (method=="MCMC") 
      {
         # Extract the required MCMC samples:

         MCMCobj <- fitObject$MCMC

         beta0MCMC <-  MCMCobj$beta0
         betaTildeMCMC <- MCMCobj$betaTilde
         gammaBetaMCMC <- MCMCobj$gammaBeta
         uTildeMCMC <- MCMCobj$uTilde
         gammaUMCMC <- MCMCobj$gammaU      
      }

      if (method=="MFVB")
      {
         # Extract the MFVB q-density parameters that are required for 
         # graphical display:
 
         MFVBobj <- fitObject$MFVB

         mu.q.beta0 <- MFVBobj$beta0[1]
         sigsq.q.beta0 <- MFVBobj$beta0[2]

         mu.q.betaTilde <- MFVBobj$betaTilde[[1]]
         Sigma.q.betaTilde <- MFVBobj$betaTilde[[2]]
         sigsq.q.betaTilde <- diag(Sigma.q.betaTilde)
         mu.q.gamma.beta <- MFVBobj$gammaBeta
      
         mu.q.uTilde <- MFVBobj$uTilde[[1]]
         sigsq.q.uTilde <- MFVBobj$uTilde[[2]]
         mu.q.gamma.u <- MFVBobj$gammaU
      }
   
      # Obtain foundational plotting  objects:

      if (method=="MCMC")
      {
         # Obtain the gridwise beta0 MCMC sample:

         beta0MC <- matrix(rep(beta0MCMC,gridSize),gridSize,length(beta0MCMC),byrow=TRUE)
      } 

      if (method=="MFVB")
      {
         # Obtain the gridwise beta0 MFVB sample:

         beta0MC <- rnorm(nMC,mu.q.beta0,sqrt(sigsq.q.beta0))

      } 

      xgList <- vector("list",numPanels)
      xMediansBetaMC <- vector("list",numPanels)
      xGridsBetaMC <- vector("list",numPanels)
      ZmediansUMC <- vector("list",numPanels)
      ZgridsUMC <- vector("list",numPanels)

      for (j in 1:numPanels)
      {
         # Obtain current "j index" value:

         jIndBetaCurr <- indsEstNonlinEff[j] 
         jIndUcurr <- jIndBetaCurr - dLinear

         # Set up plotting grid for current panel:

         xCurr <- XestNonlin[,j]
         xgCurr <- seq(min(xCurr),max(xCurr),length=gridSize) 
         xgList[[j]] <- xgCurr  

         # Obtain the current median-wise grid:

         xMediangCurr <- rep(median(xCurr),gridSize)

         if (method=="MCMC")
         {
            # Obtain the beta MCMC sample for the "j"th panel:

            betaMCcurr <- gammaBetaMCMC[,jIndBetaCurr]*betaTildeMCMC[,jIndBetaCurr] 
         } 

         if (method=="MFVB")
         {
            # Obtain the beta MFVB sample for the "j"th panel:

            betaMCcurr <- rNormalZero(nMC,mu.q.betaTilde[jIndBetaCurr],
                                      sqrt(Sigma.q.betaTilde[jIndBetaCurr,jIndBetaCurr]),
                                      mu.q.gamma.beta[jIndBetaCurr])
         } 

         # Obtain xMediansBetaMC[[j]] :

         xMediansBetaMC[[j]] <- tcrossprod(xMediangCurr,betaMCcurr)
          
         # Obtain xGridsBetaMC[[j]] :

         xGridsBetaMC[[j]] <- tcrossprod(xgCurr,betaMCcurr)

         if (method=="MCMC")
         {
            # Obtain the u MCMC sample for the "j"th panel:

            uMCcurr <- gammaUMCMC[,jIndUcurr]*uTildeMCMC[[jIndUcurr]]
         }  
         if (method=="MFVB")
         {
            # Obtain the u MFVB sample for the "j"th panel:

            tmpMat <- NULL
            for (k in 1:length(mu.q.uTilde[[jIndUcurr]]))
               tmpMat <- cbind(tmpMat,rNormalZero(nMC,mu.q.uTilde[[jIndUcurr]][k],
                                  sqrt(sigsq.q.uTilde[[jIndUcurr]][k]),mu.q.gamma.u[jIndUcurr]))
            uMCcurr <- tmpMat
         } 

         # Extract the current Z matrix attributes:

         rangexCurr <- rangexList[[jIndUcurr]]
         intKnotsCurr <- intKnotsList[[jIndUcurr]]
         OStoDRmatCurr <- OStoDRmatList[[jIndUcurr]]

         # Obtain the "j"th median-wise "Z" matrix:

         ZOSgCurr <- ZOSull(xMediangCurr,range.x=rangexCurr,intKnots=intKnotsCurr)
         COSgCurr <- cbind(1,xMediangCurr,ZOSgCurr)
         CcDRgCurr <- tcrossprod(COSgCurr,t(OStoDRmatCurr))
         ZgCurr <- CcDRgCurr[,-c(1,2)]

         if ((truncateBasis)&(ncol(ZgCurr)>=numBasis))
            ZgCurr <- ZgCurr[,1:numBasis]

         # Obtain ZmediansUMCMC[[j]]: :

         ZmediansUMC[[j]] <- tcrossprod(ZgCurr,uMCcurr[,1:ncol(ZgCurr)])

         # Obtain the "j"th gridwise "Z" matrix:
  
         ZOSgCurr <- ZOSull(xgCurr,range.x=rangexCurr,intKnots=intKnotsCurr)
         COSgCurr <- cbind(1,xgCurr,ZOSgCurr)
         CcDRgCurr <- tcrossprod(COSgCurr,t(OStoDRmatCurr))
         ZgCurr <- CcDRgCurr[,-c(1,2)]
         if ((truncateBasis)&(ncol(ZgCurr)>=numBasis))
            ZgCurr <- ZgCurr[,1:numBasis]

         # Obtain ZgridsUMC[[j]]:

         ZgridsUMC[[j]] <- tcrossprod(ZgCurr,uMCcurr[,1:ncol(ZgCurr)])
         
         # Obtain the "base median-wise" Monte Carlo sample:

         baseMediangMC <- beta0MC
         for (j in 1:numPanels)
            baseMediangMC <- baseMediangMC + xMediansBetaMC[[j]] + ZmediansUMC[[j]]
      }

      estgList <- vector("list",numPanels)
      lowgList <- vector("list",numPanels)
      uppgList <- vector("list",numPanels)
  
      # Obtain plotting vectors:

      for (j in 1:numPanels)
      {
         # Obtain plotting vectors:
  
         fHatMCcurr <- baseMediangMC + (xGridsBetaMC[[j]] - xMediansBetaMC[[j]])
         fHatMCcurr <- fHatMCcurr + (ZgridsUMC[[j]] - ZmediansUMC[[j]])

         if ((fitObject$family=="binomial")&(yscale=="response"))
            fHatMCcurr <- pnorm(fHatMCcurr)

         estgList[[j]] <- apply(fHatMCcurr,1,mean)			  
         lowgList[[j]] <- apply(fHatMCcurr,1,quantile,(1-credLev)/2) 
         uppgList[[j]] <- apply(fHatMCcurr,1,quantile,(1+credLev)/2) 
  
      }

      # Transform the "estgList", "lowgList" and "uppgList" vectors to correspond to 
      # the original input data:

      for (j in 1:numPanels)
      {
         estgList[[j]] <- sdy*estgList[[j]] + meany
         lowgList[[j]] <- sdy*lowgList[[j]] + meany
         uppgList[[j]] <- sdy*uppgList[[j]] + meany
      }
   
      # Transform the "xgList" vectors to correspond to the original input data:

      for (j in 1:numPanels)
      {
         aVal <- meanX[j]
         bVal <- sdX[j]
         XestNonlin[,j] <- aVal + bVal*XestNonlin[,j]
         xgList[[j]] <- aVal + bVal*xgList[[j]]
      }
   
      # Display all of the generalized additive model nonlinear fits:

      colVec <- rep(estCol,numPanels)
      shadeColVec <- rep(varBandCol,numPanels)
      varBandColVec <- rep(varBandCol,numPanels)

      displayAllGAMfits(XestNonlin,xgList,estgList,lowgList,uppgList,xLabsVec,yLabsVec,
                        colVec,varBandColVec,shade,shadeColVec,rug,rugSampSize,rugCol,
                        mfrow,xlimMat,ylimMat,mai,pages,cex.axis,cex.lab)
   }
 
   return(invisible())
}

############ End of plot.gamselBayes ############
