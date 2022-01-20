########## R function: effectTypes ##########

# For tabulating the estimated effect types from a 
# gamselBayes() fit.

# Last changed: 19 JAN 2022

effectTypes <- function(fitObject)
{
   # Process the predictor data:

   summGBpredProcObj <- summGBpredProc(fitObject)
   namesPredsLin <- summGBpredProcObj$namesPredsLin
   namesPredsNonlin <- summGBpredProcObj$namesPredsNonlin

   # Obtain character string summaries of the selected predictors and 
   # the nature of their effects:

   if (!is.null(namesPredsLin))
   {
      cat("\n")
      cat("   ---------------------------------------------\n")
      if (length(namesPredsLin)==1)
         cat("   Predictor selected as having a linear effect:\n",sep="")
      if (length(namesPredsLin)>1)
         cat("   Predictors selected as having linear effects:\n",sep="")
      cat("   ---------------------------------------------\n\n")
      neatCharVecPrint(sort(namesPredsLin),3,45)
   } 
   
   if (!is.null(namesPredsNonlin))
   {
      cat("\n")
      cat("   -------------------------------------------------\n")
      if (length(namesPredsNonlin)==1)
          cat("   Predictor selected as having a non-linear effect:\n")
      if (length(namesPredsNonlin)>1)
          cat("   Predictors selected as having non-linear effects:\n")
      cat("   -------------------------------------------------\n\n")
      neatCharVecPrint(sort(namesPredsNonlin),3,49)
   }  

   if (is.null(namesPredsLin)&is.null(namesPredsNonlin))
   {
      cat("\n")
      cat("   Each of the predictors have a zero estimated effect type.\n")
   } 

   cat("\n")
   return(invisible())
}

############ effectTypes ############