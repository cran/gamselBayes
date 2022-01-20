########## R function: gBdataPreProc ##########

# For conducting data pre-processing for the gamselBayes() function.

# Last changed: 08 OCT 2021

gBdataPreProc <- function(y,Xlinear,Xgeneral,family)
{
   # Standardise the predictor data depending on the value of "family"

   if (family=="gaussian")   
   {
      meany <- mean(y)
      sdy <- sd(y)
      y <- (y - meany)/sdy
   }
   if (family=="binomial")  
   {
      meany <- 0
      sdy <- 1
   }

   # Make sure that at least one of "Xlinear" or "Xgeneral" is non-null:

   if (is.null(Xlinear)&is.null(Xgeneral))
      stop("Both Xlinear and Xgeneral are NULL. Some predictor data are required.\n")

   # Make sure that each column of "Xgeneral" contains predictor data with
   # more than 5 unique values:

   if (!is.null(Xgeneral))
     if (any(apply(Xgeneral,2,vecIsLowUniq,6)))
     {
        stopStr1 <- "At least one column of \"Xgeneral\" has 5 or fewer unique values."
        stopStr2 <- "Such predictor data are not suitable for consideration as potentially"
        stopStr3 <- "having a non-linear effect. Data for such predictors should either"
        stopStr4 <- "be moved to the \"Xlinear\" argument or omitted."
        stop(paste(stopStr1,"\n",stopStr2,"\n",stopStr3,"\n",stopStr4,"\n",sep=""))
     }

   # Make sure that the numbers of rows in "Xlinear" and "Xgeneral" match the length of "y":

   if (!is.null(Xlinear))
      if (nrow(Xlinear)!=length(y)) 
         stop("The number of rows in \"Xlinear\" does not equal the response data size.\n")

   if (!is.null(Xgeneral))
      if (nrow(Xgeneral)!=length(y)) 
         stop("The number of rows in \"Xgeneral\" does not equal the response data size.\n")

   # Possibly pre-transform the "Xlinear" predictor data:

   if (is.null(Xlinear))
   {
      meanXlinear <- NULL
      sdXlinear <- NULL
   }
   if (!is.null(Xlinear))
   {
      Xlinear <- as.data.frame(Xlinear)

      # Conduct a check of non-constancy of each column:

      doConstColsChk(Xlinear,"Xlinear")

      # Store the means and standard deviations of the original "linear" predictor data:

      meanXlinear <- apply(Xlinear,2,mean)
      sdXlinear <- apply(Xlinear,2,sd)

      # Convert the "linear" predictor data to standardised form for fitting:

      Xlinear <- standardiseDataFrame(Xlinear)
   }

   # Possibly pre-transform the "Xgeneral" predictor data:

   if (is.null(Xgeneral))
   {
      meanXgeneral <- NULL
      sdXgeneral <- NULL
   }
   if (!is.null(Xgeneral))
   {
      Xgeneral <- as.data.frame(Xgeneral)

      # Conduct a check of non-constancy of each column:

      doConstColsChk(Xgeneral,"Xgeneral")

      # Store the means and standard deviations of the original "nonlinear" predictor data:

      meanXgeneral <- apply(Xgeneral,2,mean)
      sdXgeneral <- apply(Xgeneral,2,sd)
   
      # Convert the "nonlinear" predictor data to standardised form for fitting:

      Xgeneral <- standardiseDataFrame(Xgeneral)
   } 

   # Extract key dimension variables:
 
   if (is.null(Xlinear)) dLinear <- 0
   if (!is.null(Xlinear)) 
      dLinear <- ncol(Xlinear)

   if (is.null(Xgeneral)) dGeneral <- 0
   if (!is.null(Xgeneral))  
      dGeneral <- ncol(Xgeneral)
  
   # If the names  "Xlinear" are NULL or non-unique then add default names:

   if (!is.null(Xlinear))
   {
      namesCurr <- names(Xlinear)
      if (is.null(namesCurr)|(length(unique(namesCurr))<length(namesCurr)))
         names(Xlinear) <- paste("x",1:dLinear,sep="") 
   }

   # If the names of "Xgeneral" are NULL or non-unique then add default names:

   if (!is.null(Xgeneral))
   {
      namesCurr <- names(Xgeneral)
      if (is.null(namesCurr)|(length(unique(namesCurr))<length(namesCurr)))
         names(Xgeneral) <- paste("x",(dLinear+1):(dLinear+dGeneral),sep="")
   }

   # Form the X matrix:

   if (is.null(Xlinear)) X <- as.matrix(Xgeneral)
   if (is.null(Xgeneral)) X <- as.matrix(Xlinear)
   if ((!is.null(Xlinear))&(!is.null(Xgeneral))) 
      X <- as.matrix(cbind(Xlinear,Xgeneral))

   # Return the processed inputs:

   return(list(y=y,Xlinear=Xlinear,Xgeneral=Xgeneral,X=X,dLinear=dLinear,dGeneral=dGeneral,
               meany=meany,sdy=sdy,meanXlinear=meanXlinear,sdXlinear=sdXlinear,
               meanXgeneral=meanXgeneral,sdXgeneral=sdXgeneral))
}

############ End of gBdataPreProc ############
