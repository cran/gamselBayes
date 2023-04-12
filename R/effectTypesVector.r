########## R function: effectTypesVector ##########

# For obtaining the effect types from a gamselBayes() fit object.

# Last changed: 06 OCT 2021

effectTypesVector <- function(fitObject)
{
   # Extract the vector of estimated effect types:

   outVec <- fitObject$effectTypesHat

   # Add the names to the vector:

   namesVec <- NULL
   if (!is.null(fitObject$Xlinear))
      namesVec <- c(namesVec,names(fitObject$Xlinear))
   if (!is.null(fitObject$Xgeneral))
      namesVec <- c(namesVec,names(fitObject$Xgeneral))
   names(outVec) <- namesVec
 
   # Return the fully adorned effect types vector:

   return(outVec)
}
   
############ End effectTypesVector ###########

