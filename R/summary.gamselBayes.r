########## R function: summary.gamselBayes ##########

# For summarising a gamselBayes() fit object.

# Last changed: 08 NOV 2021

summary.gamselBayes <- function(object,credLev=0.95,sigFigs=5,nMC=10000,...)
{
   # Process the predictor data inputs:

   summGBpredProcObj <- summGBpredProc(object)
   dLinear <- summGBpredProcObj$dLinear
   dGeneral <- summGBpredProcObj$dGeneral
   indsLinEffect <- summGBpredProcObj$indsLinEffect
   namesPredsLin <- summGBpredProcObj$namesPredsLin
   namesPredsNonlin <- summGBpredProcObj$namesPredsNonlin

   # Process the other inputs:

   summGBargsProcObj <- summGBargsProc(credLev,sigFigs,nMC)
   credLev <- summGBargsProcObj$credLev
   sigFigs <- summGBargsProcObj$sigFigs
   nMC <- summGBargsProcObj$nMC

   # Determine the value of "addMFVBwarn":

   addMFVBwarn <- (object$method=="MFVB")&(object$family=="binomial")

   # Obtain the scale factors that accompany the estimated linear effects:

   linCoefScaFacs <- getLinCoefScaFacs(object)

   # Make summary table for the linear effects:

   linCoefSumm <- summGBlinCoef(object,indsLinEffect,dLinear,dGeneral,
                                linCoefScaFacs,credLev,sigFigs,nMC)
   meanVec <- linCoefSumm$meanVec
   lowVec <- linCoefSumm$lowVec
   uppVec <- linCoefSumm$uppVec
   linCoefTable <- summGBcoefTable(namesPredsLin,credLev,meanVec,lowVec,uppVec,
                                   addMFVBwarn)

   # Return the linear coefficients table:

   if (!is.null(linCoefTable)) return(linCoefTable)
   if (is.null(linCoefTable)) return(invisible())
}   
 
############ End of summary.gamselBayes ############
