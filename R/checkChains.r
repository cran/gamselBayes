########## R function: checkChains ##########

# For graphically checking the Markov chain Monte Carlo 
# samples ("chains") of a gamselBayes() fit object.

# Last changed: 14 JAN 2022

checkChains <- function(fitObject,colourVersion=TRUE,paletteNum=1)
{
   # Ensure that the "colourVersion" input is legal:

   if (!(any(colourVersion==c(FALSE,TRUE))))
   {
      warnStr1 <- "The inputted flag for producing a colour version is neither"
      warnStr2 <- "  FALSE nor TRUE. The default value of TRUE was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      colourVersion <- TRUE
   }

   # Ensure that the "paletteNum" input is legal:

   if (!(any(paletteNum==c(1,2))))
   {
      warnStr1 <- "The inputted palette number is neither 1 nor 2. The default"
      warnStr2 <- "  value of 1 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      paletteNum <- 1
   }

   # Extract the "method" value:

   method <- fitObject$method

   # If "method" is MFVB then print an informational messages,
   # otherwise obtain relevant MCMC samples and graphically
   # display them:

   if (method=="MFVB")
   {
      message("The checkChains() function only applies to the method equalling \"MCMC\".")
   }

   if (method=="MCMC")
   {
      # Obtain effect type and parameter and estimates depending
      # on the method type:

      effectTypesHat <- fitObject$effectTypesHat

      # Determine the indices of the predictors which have a 
      # linear effect estimate:

      indsEstLinEff <- (1:length(effectTypesHat))[effectTypesHat=="linear"]

      # Determine the indices of the predictors which have a 
      # non-linear effect estimate:

      indsEstNonlinEff <- (1:length(effectTypesHat))[effectTypesHat=="nonlinear"]

      # Determine the numbers of MCMC samples of each type:

      numMCMCsampsLin <- length(indsEstLinEff)
      numMCMCsampsNonlin <- length(indsEstNonlinEff)

      # Determine the total number of MCMC samples:

      numMCMCsamps <- numMCMCsampsLin + numMCMCsampsNonlin

      if (numMCMCsamps==0)
      {
         message("There are no Markov chain Monte Carlo samples to view.")
      }
      if (numMCMCsamps>0)
      {
         # Extract the predictor values:

         Xlinear <- fitObject$Xlinear
         Xgeneral <- fitObject$Xgeneral

         # Obtain the combined predictor matrix:

         if (is.null(Xlinear)) X <- as.matrix(Xgeneral)
         if (is.null(Xgeneral)) X <- as.matrix(Xlinear)
         if ((!is.null(Xlinear))&(!is.null(Xgeneral))) 
            X <- as.matrix(cbind(Xlinear,Xgeneral))

         # If the sample size exceeeds the nominal value
         # of 1000 observations then thin the data for
         # chain assessment purposes:

         sampSizeVal <- nrow(X)
         if (sampSizeVal>1000)
         {
            subInds <- sample(1:sampSizeVal,1000,replace=FALSE)
            X <- as.matrix(X[subInds,])  
         }     
      
         # Extract the "dLinear", "dGeneral" values:

         if (is.null(Xlinear))     dLinear <- 0
         if (!is.null(Xlinear))    dLinear <- ncol(Xlinear)
         if (is.null(Xgeneral))  dGeneral <- 0
         if (!is.null(Xgeneral)) dGeneral <- ncol(Xgeneral)
 
         # Determine the vector of labels for each MCMC sample:

         if (dLinear==0) predNamesFull <- names(Xgeneral)
         if (dGeneral==0) predNamesFull <- names(Xlinear)
         if ((dLinear>0)&(dGeneral>0)) predNamesFull <- c(names(Xlinear),names(Xgeneral))
  
         predNames <- NULL
         if (numMCMCsampsLin>0)
            predNames <-  c(predNames,predNamesFull[indsEstLinEff])
         if (numMCMCsampsNonlin>0)
            predNames <-  c(predNames,predNamesFull[indsEstNonlinEff])

         # Obtain the "beta" samples for the (non-zero) linear effects, assuming
         # that there are some such samples:

         if (numMCMCsampsLin==0)
            MCMCsampsToPlotLin <- NULL 

         if (numMCMCsampsLin>0)
         {
            betaTildeMCMC <- fitObject$MCMC$betaTilde[,indsEstLinEff]
            gammaBetaMCMC <- fitObject$MCMC$gammaBeta[,indsEstLinEff]
            MCMCsampsToPlotLin <- as.matrix(gammaBetaMCMC*betaTildeMCMC)
         }
   
         # Obtain the vertical slice at median samples for the non-linear
         # effects, assuming that there are some such samples:

         if (numMCMCsampsNonlin==0)
            MCMCsampsToPlotNonlin <- NULL 

         if (numMCMCsampsNonlin>0)
         {
        
            # Subset the columns of X to correspond to those predictors for 
            # which there is an estimated non-linear effect:

            X <- as.matrix(X[,indsEstNonlinEff])
   
            # Obtain version of predictor design matrix for plotting, noting
            # that only those corresponding to the effectTypesHat=="nonlinear"
            # indices are needed:

            # Extract spline basis function parameters:

            rangexList <- fitObject$rangex
            intKnotsList <- fitObject$intKnots
            truncateBasis <- fitObject$truncateBasis
            numBasis <- fitObject$numBasis

            # Determine the index corresponding to the median position
            # of the sorted predictor samples:

            indMedian <- round((nrow(X)+1)/2)

            # Form the spline basis functions:

            Zlist <- vector("list",numMCMCsampsNonlin)
            ncZvec <- NULL
            for (j in 1:numMCMCsampsNonlin)
            {
               xCurr <- sort(X[,j])
               Zlist[[j]] <- ZcDR(xCurr,rangexList[[indsEstNonlinEff[j]-dLinear]],
                                  intKnotsList[[indsEstNonlinEff[j]-dLinear]])
               if ((truncateBasis)&(ncol(Zlist[[j]])>=numBasis))
                  Zlist[[j]] <- Zlist[[j]][,1:numBasis]
               ncZvec <- c(ncZvec,ncol(Zlist[[j]]))
            }

            # Extract the required MCMC samples:

            MCMCobj <- fitObject$MCMC
            beta0MCMC <-  MCMCobj$beta0
            betaTildeMCMC <- MCMCobj$betaTilde[,indsEstNonlinEff]
            gammaBetaMCMC <- MCMCobj$gammaBeta[,indsEstNonlinEff]
            betaMCMC <- as.matrix(gammaBetaMCMC*betaTildeMCMC)
            uTildeMCMC <- MCMCobj$uTilde
            gammaUMCMC <- MCMCobj$gammaU
 
            # Form the required MCMC samples for the coefficients:

            uMCMC <- vector("list",numMCMCsampsNonlin)
            for (j in 1:numMCMCsampsNonlin)
               uMCMC[[j]] <- gammaUMCMC[[indsEstNonlinEff[j]-dLinear]]*uTildeMCMC[[indsEstNonlinEff[j]-dLinear]]
      
            # Set up the matrix for storing the MCMC samples to be graphically assessed:
       
            numMCMC <- length(beta0MCMC)
            MCMCsampsToPlotNonlin <- matrix(NA,numMCMC,numMCMCsampsNonlin)

            for (j in 1:numMCMCsampsNonlin)
            {
               # Obtain current abscissae vector:

               xCurr <- sort(X[,j])

               # Obtain plotting vectors:
  
               fHatMCMCcurr <- tcrossprod(xCurr,betaMCMC[,j]) 
               fHatMCMCcurr <- fHatMCMCcurr + tcrossprod(Zlist[[j]],uMCMC[[j]])
       
               MCMCsampsToPlotNonlin[,j] <- fHatMCMCcurr[indMedian,] 
            }
         }
      }

      # Combine the matrices of MCMC samples of both type:

      if (is.null(MCMCsampsToPlotLin)) MCMCsampsToPlot <- MCMCsampsToPlotNonlin
      if (is.null(MCMCsampsToPlotNonlin)) MCMCsampsToPlot <- MCMCsampsToPlotLin
      if ((!is.null(MCMCsampsToPlotLin))&(!is.null(MCMCsampsToPlotNonlin)))
         MCMCsampsToPlot <- as.matrix(cbind(MCMCsampsToPlotLin,MCMCsampsToPlotNonlin))
  
      # Re-order the samples alphabetically:

      alphaOrder <- order(predNames)
      predNames <- predNames[alphaOrder]
      MCMCsampsToPlot <- MCMCsampsToPlot[,alphaOrder]

      # Determine number of views for plotting:

      numMCMCsampsToPlot <- ncol(MCMCsampsToPlot)
      numSampsPerView <- 6
      numViews <- ceiling(numMCMCsampsToPlot/numSampsPerView)

      # Loop through the views:

      indsStt <- 1
      for (iView in 1:numViews)
      {
         indsEnd <- indsStt + numSampsPerView - 1
         if (indsEnd>numMCMCsampsToPlot)
            indsEnd <- numMCMCsampsToPlot
         indsCurr <- indsStt:indsEnd
         summChainsGamsel(as.matrix(MCMCsampsToPlot[,indsCurr]),predNames[indsCurr],
                          colourVersion=colourVersion,paletteNum=paletteNum)
         indsStt <- indsEnd + 1
         if (indsEnd<numMCMCsampsToPlot)
            readline("Hit Enter to view the next checkChains() plot.\n")
      }
   }

   return(invisible())
}

############ End of checkChains ############
