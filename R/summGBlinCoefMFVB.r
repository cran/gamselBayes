########## R function: summGBlinCoefMFVB ##########

# For obtaining the linear coefficient summary table and
# output vector for summary.gamselBayes() when the method 
# is "MFVB".

# Last changed: 06 OCT 2021

summGBlinCoefMFVB <- function(object,indsLinEffect,dLinear,dGeneral,
                              linCoefScaFacs,credLev,sigFigs,nMC)
{
   # Obtain the units-corrected q*(beta) parameters for the 
   # "Xlinear" predictors:
   
   if (dLinear==0)
   { 
      mu.q.gammaBetaLin <- NULL
      mu.q.betaTildeLin <- NULL
      sigma.q.betaTildeLin <- NULL
   }
   if (dLinear>0)
   {
      # Extract q-density parameters corresponding to the standardised data:
   
      mu.q.gammaBetaLin <- object$MFVB$gammaBeta[1:dLinear]
      mu.q.betaTildeLin <- object$MFVB$betaTilde[[1]][1:dLinear]
      sigma.q.betaTildeLin <- sqrt(diag(object$MFVB$betaTilde[[2]])[1:dLinear])
   } 
   
   # Obtain the units-corrected q*(beta) parameters for the "Xgeneral" predictors:
   
   if (dGeneral==0)
   {
      mu.q.gammaBetaNonlin <- NULL
      mu.q.betaTildeNonlin <- NULL
      sigma.q.betaTildeNonlin <- NULL
   }
   if (dGeneral>0)
   {
      # Extract q-density parameters corresponding to the standardised data:
   
      mu.q.gammaBetaNonlin <- object$MFVB$gammaBeta[(dLinear+1):(dLinear+dGeneral)]
      mu.q.betaTildeNonlin <- object$MFVB$betaTilde[[1]][(dLinear+1):(dLinear+dGeneral)]
      sigma.q.betaTildeNonlin <- sqrt(diag(object$MFVB$betaTilde[[2]])[(dLinear+1):(dLinear+dGeneral)])
   }
   
   # Combine to obtain the full units-corrected q-density parameter vectors:
   
   mu.q.betaTilde <- c(mu.q.betaTildeLin,mu.q.betaTildeNonlin)
   sigma.q.betaTilde <- c(sigma.q.betaTildeLin,sigma.q.betaTildeNonlin)
   mu.q.gammaBeta <- c(mu.q.gammaBetaLin,mu.q.gammaBetaNonlin)     

   # Reduce the  units-corrected q-density parameter vectors to correspond 
   # to the predictors for which the effect is linear:
   
   mu.q.gammaBeta <- mu.q.gammaBeta[indsLinEffect]
   mu.q.betaTilde <- mu.q.betaTilde[indsLinEffect]
   sigma.q.betaTilde <- sigma.q.betaTilde[indsLinEffect]

   # Obtain the vectors of posterior means and credible interval limits:

   numSumms <- length(mu.q.gammaBeta)
   meanVec <- mu.q.gammaBeta*mu.q.betaTilde
   lowVec <- rep(NA,numSumms) ; uppVec <- rep(NA,numSumms)
   for (j in 1:numSumms)
   {   
      lowVec[j] <- quantile(rNormalZero(nMC,mu.q.betaTilde[j],sigma.q.betaTilde[j],
                            mu.q.gammaBeta[j]),(1-credLev)/2) 
      uppVec[j] <- quantile(rNormalZero(nMC,mu.q.betaTilde[j],sigma.q.betaTilde[j],
                            mu.q.gammaBeta[j]),(1+credLev)/2)
   }

   # Transform to correspond to the original units:

   meanVec <- signif(meanVec/linCoefScaFacs,sigFigs)
   lowVec <-  signif(lowVec/linCoefScaFacs,sigFigs)
   uppVec <-  signif(uppVec/linCoefScaFacs,sigFigs)

   # Return estimation and credible interval vectors:

   return(list(meanVec=meanVec,lowVec=lowVec,uppVec=uppVec))
}

############ summGBlinCoefMFVB ############