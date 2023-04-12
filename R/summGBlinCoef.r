########## R function: summGBlinCoef ##########

# For obtaining the linear coefficient summary table and
# output vector for summary.gamselBayes().

# Last changed: 19 JAN 2022

summGBlinCoef <- function(object,indsLinEffect,dLinear,dGeneral,
                          linCoefScaFacs,credLev,sigFigs,nMC)
{
   if (length(indsLinEffect)==0)
   {
      message("There are no linear effects to summarise.\n")
      return(NULL)
   }

   if (length(indsLinEffect)>0)
   {
      if (object$method=="MCMC")
         summGBlinCoefObj <- summGBlinCoefMCMC(object,indsLinEffect,dLinear,dGeneral,
                                               linCoefScaFacs,credLev,sigFigs)


      if (object$method=="MFVB")
         summGBlinCoefObj <- summGBlinCoefMFVB(object,indsLinEffect,dLinear,dGeneral,
                                               linCoefScaFacs,credLev,sigFigs,nMC)
   }

   return(summGBlinCoefObj)
}

############ summGBlinCoef ############
