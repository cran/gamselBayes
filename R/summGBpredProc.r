########## R function: summGBpredProc ##########

# For conducting predictor data processing for the 
# summary.gamselBayes() function.

# Last changed: 06 OCT 2021

summGBpredProc <- function(object)
{
   # Extract key variables:
  
   Xlinear <- object$Xlinear
   Xgeneral <- object$Xgeneral
   effectTypesHat <- object$effectTypesHat

   # Extract the full names vector:

   predNames <- c(names(Xlinear),names(Xgeneral))

   # Determine the "dLinear" and "dGeneral" values and the number of predictors:

   if (is.null(Xlinear))     dLinear <- 0
   if (!is.null(Xlinear))    dLinear <- ncol(Xlinear)
   if (is.null(Xgeneral))  dGeneral <- 0
   if (!is.null(Xgeneral)) dGeneral <- ncol(Xgeneral)
   numPred <- dLinear + dGeneral

   # Obtain the indices of those predictors determined to be linear:

   indsLinEffect <- (1:length(effectTypesHat))[effectTypesHat=="linear"]

   # Obtain the indices of those predictors determined to be non-linear:

   indsNonlinEffect <- (1:length(effectTypesHat))[effectTypesHat=="nonlinear"]

   # Obtain the names of predictors of each non-zero type:

   if (length(indsLinEffect)==0)    namesPredsLin <- NULL
   if (length(indsLinEffect)>0)     namesPredsLin <- predNames[indsLinEffect]
   if (length(indsNonlinEffect)==0) namesPredsNonlin <- NULL
   if (length(indsNonlinEffect)>0)  namesPredsNonlin <- predNames[indsNonlinEffect]

   # Return the processed inputs:

   return(list(dLinear=dLinear,dGeneral=dGeneral,indsLinEffect=indsLinEffect,
               namesPredsLin=namesPredsLin,namesPredsNonlin=namesPredsNonlin))
}

############ summGBpredProc ############