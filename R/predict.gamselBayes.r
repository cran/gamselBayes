########## R function: predict.gamselBayes ##########

# For prediction based on new data from a gamselBayes() 
# fit object.

# Last changed: 11 APR 2023

predict.gamselBayes <- function(object,newdata,type="response",...)
{
   fitObject <- object

   # Process the newdata input:
   
   processedInput <- predictGBnewdataProc(newdata,fitObject)
   XlinearNew <- processedInput$XlinearNew
   XgeneralNew <- processedInput$XgeneralNew
   XNew <- processedInput$XNew 
   dLinear <- processedInput$dLinear
   dGeneral <- processedInput$dGeneral
   lenPred <- processedInput$lenPred
   
   # Check the legality of the "type" argument:

   if (!any(type==c("link","response","terms")))
      stop("The type must be one of \"link\", \"response\" or \"terms\"\n.")
   
   # Extract the method and family:

   method <- fitObject$method
   family <- fitObject$family
  
   if (method=="MCMC") # Extract the number of MCMC samples:
      numMCMC <- length(fitObject$MCMC$beta0)

   # Extract the orginal response variable mean and standard deviation:

   meany <- fitObject$meany
   sdy <- fitObject$sdy

   # Obtain the intercept estimate for the case of family being "binomial":

   if (family=="binomial")
   {
      if (method=="MCMC") 
         beta0hat <- mean(fitObject$MCMC$beta0)
         
      if (method=="MFVB") 
         beta0hat <- fitObject$MFVB$beta0[1]
   }

   # Set up a "lenPred" by ("dLinear"+"dGeneral") matrix for storing the term-wise
   # predictions, starting with zeroes in all positions:

   predMat <- matrix(0,lenPred,(dLinear+dGeneral))

   # Obtain effect type and parameter and estimates depending
   # on the method type:

   effectTypesHat <- fitObject$effectTypesHat

   # Determine the indices of the predictors which have a linear 
   # effect estimate:

   indsEstLinEff <- (1:length(effectTypesHat))[effectTypesHat=="linear"]
   numLinSelected <- length(indsEstLinEff)

   # Determine the indices of the predictors which have a 
   # non-linear effect estimate:

   indsEstNonlinEff <- (1:length(effectTypesHat))[effectTypesHat=="nonlinear"]
   numNonlinSelected <- length(indsEstNonlinEff)

   if (method=="MCMC")
   {
      # Obtain the full "betaMCMC" matrix:

      betaTildeMCMC <- as.matrix(fitObject$MCMC$betaTilde)
      gammaBetaMCMC <- as.matrix(fitObject$MCMC$gammaBeta)
      betaMCMC <- gammaBetaMCMC*betaTildeMCMC
   }

   if (method=="MFVB")
   {
      # Obtain the full "mu.q.beta" vector:

      mu.q.betaTilde <- fitObject$MFVB$betaTilde$mu.q.betaTilde
      mu.q.gamma.beta <- fitObject$MFVB$gammaBeta
      mu.q.beta <- mu.q.gamma.beta*mu.q.betaTilde
   }

   if (numLinSelected>0)
   {
      # Process the linear effect contributions:

      for (jLin in 1:numLinSelected)
      {
         xNewCurr <- XNew[,indsEstLinEff[jLin]]
         if (method=="MCMC")
         {
            betaMCMCcurr <- betaMCMC[,indsEstLinEff[jLin]]
            termMCMCcurr <- crossprod(t(betaMCMCcurr),xNewCurr)
            predMat[,indsEstLinEff[jLin]] <- apply(termMCMCcurr,2,mean)
         }
         if (method=="MFVB")
         {
            mu.q.betaCurr <- mu.q.beta[indsEstLinEff[jLin]]
            termMFVBcurr <- as.vector(crossprod(mu.q.betaCurr,xNewCurr))
            predMat[,indsEstLinEff[jLin]] <- termMFVBcurr
         }
      }
   }
 
   if (numNonlinSelected>0)
   {
      # Extract spline basis function parameters:

      rangexList <- fitObject$rangex
      intKnotsList <- fitObject$intKnots
      OStoDRmatList <- fitObject$OStoDRmat
      truncateBasis <- fitObject$truncateBasis
      numBasis <- fitObject$numBasis

      for (jNonlin in 1:numNonlinSelected)
      {
         # Initialise the current term for the full data:

         termCurrFull <- rep(0,nrow(XgeneralNew))

         # Add the current linear effect contributions to the "etaPredMCMC" matrix:

         xNewCurr <- XNew[,indsEstNonlinEff[jNonlin]]
          
         if (method=="MCMC")
         {
            betaMCMCcurr <- as.matrix(betaMCMC[,indsEstNonlinEff[jNonlin]-dLinear])
            termMCMCcurr <- crossprod(t(betaMCMCcurr),xNewCurr)
            termCurrFull <- termCurrFull + apply(termMCMCcurr,2,mean)
         }
         if (method=="MFVB")
         {
            mu.q.betaCurr <- mu.q.beta[indsEstNonlinEff[jNonlin]-dLinear]
            termMFVBcurr <- as.vector(crossprod(mu.q.betaCurr,xNewCurr))
            termCurrFull <- termCurrFull + termMFVBcurr 
         }

         # Determine spline basis for the current original data:
          
         rangexCurr <- rangexList[[indsEstNonlinEff[jNonlin]-dLinear]]
         intKnotsCurr <- intKnotsList[[indsEstNonlinEff[jNonlin]-dLinear]]
         OStoDRmatCurr <- OStoDRmatList[[indsEstNonlinEff[jNonlin]-dLinear]]
         ZOScurr <- ZOSull(xNewCurr,range.x=rangexCurr,intKnots=intKnotsCurr)
         COScurr <- cbind(1,xNewCurr,ZOScurr)
         CcDRcurr <- tcrossprod(COScurr,t(OStoDRmatCurr))
         Zcurr <- CcDRcurr[,-c(1,2)]
         if ((truncateBasis)&(ncol(Zcurr)>=numBasis))
            Zcurr <- Zcurr[,1:numBasis]

         # Check whether any of the new values are outside the spline basis
         # range and, if so, issue a warning:

         outsideRange <- FALSE
         if (min(xNewCurr)<rangexCurr[1]) outsideRange <- TRUE
         if (max(xNewCurr)>rangexCurr[2]) outsideRange <- TRUE

         if (outsideRange) 
         {
            warnStr1 <- paste("For predictor ",jNonlin," among those having a non-linear effect",sep="")
            warnStr2 <- "some new data values are beyond the spline basis range."
            warnStr3 <- "Beyond range predictions are subject to poor quality."
            warning(paste(warnStr1,"\n  ",warnStr2,"\n  ",warnStr3,"\n",sep=""),
                 immediate.=TRUE) 
         }
         
         # Add on the spline contribution for the current "newdata" vector:

         if (method=="MCMC") 
         {
            # Extract and obtain required coefficients:

            uTildeMCMCcurr <- fitObject$MCMC$uTilde[[indsEstNonlinEff[jNonlin]-dLinear]]
            gammaUMCMCcurr <- fitObject$MCMC$gammaU[[indsEstNonlinEff[jNonlin]-dLinear]]
            uMCMCcurr <- gammaUMCMCcurr*uTildeMCMCcurr

            # Obtain the estimated component over the full data:

            termMCMCcurr <- tcrossprod(uMCMCcurr[,1:ncol(Zcurr)],Zcurr)

            # Add on the new prediction term:

            termCurrFull <- termCurrFull + apply(termMCMCcurr,2,mean)
         }

         if (method=="MFVB") 
         {
           # Extract and obtain required coefficients:

            mu.q.uTildeCurr <- fitObject$MFVB$uTilde[[1]][[indsEstNonlinEff[jNonlin]-dLinear]]
            mu.q.gammaCurr <- fitObject$MFVB$gammaU[[indsEstNonlinEff[jNonlin]-dLinear]]
            mu.q.uCurr <- mu.q.gammaCurr*mu.q.uTildeCurr
             
            # Obtain the estimated component over the full data:

            termMFVBcurr <- as.vector(crossprod(t(Zcurr),mu.q.uCurr[1:ncol(Zcurr)]))
             
            # Add on the new prediction term:

            termCurrFull <- termCurrFull + termMFVBcurr
         }

         # Insert the current term into "predMat":

         predMat[,indsEstNonlinEff[jNonlin]] <-  termCurrFull 
      }
   }

   # Scale "predMat" to correspond to the original response variable units:

   predMat <- sdy*predMat

   # If "type" is "terms" then return the term-wise predictors matrix:

   if (type=="terms")
   {
      # Convert "predMat" to data frame form and add on the appropriate names:
      
      predMat <- as.data.frame(predMat)
      if (is.null(XlinearNew)) predMatNames <- names(XgeneralNew)
      if (is.null(XgeneralNew)) predMatNames <- names(XlinearNew)
      if ((!is.null(XlinearNew))&(!is.null(XgeneralNew))) 
         predMatNames <- c(names(XlinearNew),names(XgeneralNew))

      names(predMat) <- predMatNames

      # Add on the intercept prediction as an attribute:

      if (family=="gaussian")
         attr(predMat,"intercept") <- meany

      if (family=="binomial")
         attr(predMat,"intercept") <- beta0hat
      
      # Return the attribute-adorned term-wise prediction matrix:

      return(predMat)
   }

   if (type!="terms")
   {
      # If "type" is not "terms" then combine across terms into a vector
      # and return either on the link or response scale, depending on
      # whether "type" is "link" or "response":

      predVec <- as.vector(apply(predMat,1,sum))

      if (family=="gaussian")
         predVec <- predVec + meany

      if (family=="binomial")
         predVec <- predVec + beta0hat
 
      if ((type=="response")&(family=="binomial"))
         predVec <- pnorm(predVec) 
    
      return(predVec)
   }
}

############ End of predict.gamselBayes ############
