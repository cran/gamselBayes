########## R function: getLinCoefScaFacs ##########

# For obtaining the scale factors that accompany the estimated linear
# effects of a gamselBayes() fit object.

# Last changed: 08 OCT 2021

getLinCoefScaFacs <- function(fitObject)
{
   # Extract the original data standard deviations values:

   sdy <- fitObject$sdy
   sdXlinear <- fitObject$sdXlinear
   sdXgeneral <- fitObject$sdXgeneral

   # Extract the estimated effects type vector:

   effectTypesHat <- fitObject$effectTypesHat

   # Determine the indices of the estimated linear effects: 

   indsLinEffect <- (1:length(effectTypesHat))[effectTypesHat=="linear"]

   # Treat the case where there are no estimated linear effects:

   if (length(indsLinEffect)==0)
      return(NULL)

   # Treat the case where there are estimated linear effects:

   if (length(indsLinEffect)>0)
   {
      if (is.null(sdXlinear)) sdFull <- sdXgeneral/sdy
      if (is.null(sdXgeneral)) sdFull <- sdXlinear/sdy
      if ((!is.null(sdXlinear))&(!is.null(sdXgeneral)))
         sdFull <- c(sdXlinear,sdXgeneral)/sdy

      fullScaFac <- as.numeric(sdFull)

      return(fullScaFac[indsLinEffect])
   }
}

############ End of getLinCoefScaFacs ############
