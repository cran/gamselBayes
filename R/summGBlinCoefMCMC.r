########## R function: summGBlinCoefMCMC ##########

# For obtaining the linear coefficient summary table and
# output vector for summary.gamselBayes() when the 
# method is "MCMC".

# Last changed: 08 NOV 2021

summGBlinCoefMCMC <- function(object,indsLinEffect,dLinear,dGeneral,
                              linCoefScaFacs,credLev,sigFigs)
{
   # Obtain the betaMCMC matrix for the "Xlinear" predictors:
   
   if (dLinear==0) betaMCMClin <- NULL
   if (dLinear>0)
   {
      gammaBetaCurr <- object$MCMC$gammaBeta[,1:dLinear]
      betaTildeCurr <- object$MCMC$betaTilde[,1:dLinear]
      betaMCMClin <- gammaBetaCurr*betaTildeCurr
   }
   
   # Obtain the betaMCMC matrix for the "Xlinear" predictors:
   
   if (dGeneral==0) betaMCMCnonlin <- NULL
   if (dGeneral>0)
   {
      gammaBetaCurr <- object$MCMC$gammaBeta[,(dLinear+1):(dLinear+dGeneral)]
      betaTildeCurr <- object$MCMC$betaTilde[,(dLinear+1):(dLinear+dGeneral)]
      betaMCMCnonlin <- gammaBetaCurr*betaTildeCurr
   }
   
   # Combine to obtain the full betaMCMC matrix:
   
   if (is.null(betaMCMClin)) betaMCMC <- as.matrix(betaMCMCnonlin)
   if (is.null(betaMCMCnonlin)) betaMCMC <- as.matrix(betaMCMClin)
   if ((!is.null(betaMCMClin))&(!is.null(betaMCMCnonlin)))
      betaMCMC <- as.matrix(cbind(betaMCMClin,betaMCMCnonlin))
   
   # Reduce the "betaMCMC" matrix to the predictors for which
   # the effect is linear:
   
   betaMCMC <- as.matrix(betaMCMC[,indsLinEffect])
   
   # Obtain and return the vectors of posterior means and credible
   # interval limits:
   
   meanVec <- apply(betaMCMC,2,mean)
   lowVec <- apply(betaMCMC,2,quantile,(1-credLev)/2)
   uppVec <- apply(betaMCMC,2,quantile,(1+credLev)/2)

   # Transform to correspond to the original units:

   meanVec <- signif(meanVec/linCoefScaFacs,sigFigs)
   lowVec <- signif(lowVec/linCoefScaFacs,sigFigs)
   uppVec <- signif(uppVec/linCoefScaFacs,sigFigs)

   # Return estimation and credible interval vectors:

   return(list(meanVec=meanVec,lowVec=lowVec,uppVec=uppVec))
}

############ summGBlinCoefMCMC ############