########## R function: predictGBnewdataProc ##########

# For conducting new data processing for the predict.gamselBayes() 
# function.

# Last changed: 12 OCT 2021

predictGBnewdataProc <- function(newdata,fitObject)
{
   # Make sure that "newdata" is a two-component list with the 
   # correct component names:

   if ((!is.list(newdata))|(length(newdata)!=2))
   {
      stopStr1 <- "The \"newdata\" must be a list with two components named"
      stopStr2 <- "\"Xlinear\" and \"Xgeneral\"." 
      stop(paste(stopStr1,"\n",stopStr2,"\n",sep=""))
   }

   # Extract the components of "newdata":

   XlinearNew <- newdata[[1]]
   XgeneralNew <- newdata[[2]]

   # Make sure that the names of the "newdata" data frames match 
   # those of "Xlinear" and "Xgeneral":

   if (length(setdiff(names(fitObject$Xlinear),names(XlinearNew)))!=0)
      stop("The names of \"Xlinear\" in \"newdata\" and the fit object must be identical.\n")

   if (length(setdiff(names(fitObject$Xgeneral),names(XgeneralNew)))!=0)
      stop("The names of \"Xgeneral\" in \"newdata\" and the fit object must be identical.\n")

   # Make sure that that "XlinearNew" and "XgeneralNew" have the correct
   # NULL versus non-NULL forms with respect to the fit object:

   if ((is.null(fitObject$Xlinear))&(!is.null(XlinearNew)))
      stop("The \"Xlinear\" component of \"newdata\" must be NULL for this fit object.\n")

   if ((is.null(fitObject$Xgeneral))&(!is.null(XgeneralNew)))
      stop("The \"Xgeneral\" component of \"newdata\" must be NULL for this fit object.\n")

   # Make sure that that "XlinearNew" and "XgeneralNew" have the correct
   # numbers of colums with respect to the fit object:

   if (!is.null(XlinearNew))
      if (ncol(XlinearNew)!=ncol(fitObject$Xlinear))
      {
         stopStr1 <- "The number of columns in the \"Xlinear\" component of \"newdata\""
         stopStr2 <- "differs from its fit object counterpart."
         stop(paste(stopStr1,"\n",stopStr2,"\n",sep=""))
      }

   if (!is.null(XgeneralNew))
      if (ncol(XgeneralNew)!=ncol(fitObject$Xgeneral))
      {
         stopStr1 <- "The number of columns in the \"Xgeneral\" component of \"newdata\""
         stopStr2 <- "differs from its fit object counterpart."
         stop(paste(stopStr1,"\n",stopStr2,"\n",sep=""))
      }
 
   # Transform the "XlinearNew" values to match the transformed "Xlinear" data:

   if (!is.null(XlinearNew))
   {
      namesSave <- names(XlinearNew)

      meanXlinear <- fitObject$meanXlinear
      sdXlinear <- fitObject$sdXlinear     

      shiftMat <- matrix(rep(meanXlinear,nrow(XlinearNew)),
                         nrow(XlinearNew),ncol(XlinearNew),byrow=TRUE)
      
      XlinearNew <- as.data.frame(t(t(XlinearNew- shiftMat)/sdXlinear))

      names(XlinearNew) <- namesSave
   }
   
   # Transform the "XgeneralNew" values to match the transformed "Xgeneral" data:

   if (!is.null(fitObject$Xgeneral))
   {
      namesSave <- names(XgeneralNew)

      meanXgeneral <- fitObject$meanXgeneral
      sdXgeneral <- fitObject$sdXgeneral      

      shiftMat <- matrix(rep(meanXgeneral,nrow(XgeneralNew)),
                         nrow(XgeneralNew),ncol(XgeneralNew),byrow=TRUE)

      XgeneralNew <- as.data.frame(t(t(XgeneralNew - shiftMat)/sdXgeneral))
      
      names(XgeneralNew) <- namesSave 
   } 

   # Extract key dimension variables:
 
   if (is.null(XlinearNew)) dLinear <- 0
   if (!is.null(XlinearNew))
   { 
      dLinear <- ncol(XlinearNew)
      lenPred <- nrow(XlinearNew)
   } 

   if (is.null(XgeneralNew)) dGeneral <- 0
   if (!is.null(XgeneralNew))  
   {
      dGeneral <- ncol(XgeneralNew)
      lenPred <- nrow(XgeneralNew)
   }

   # Form the XNew matrix:

   if (is.null(XlinearNew)) XNew <- as.matrix(XgeneralNew)
   if (is.null(XgeneralNew)) XNew <- as.matrix(XlinearNew)
   if ((!is.null(XlinearNew))&(!is.null(XgeneralNew))) 
      XNew <- as.matrix(cbind(XlinearNew,XgeneralNew))

   # Return the processed inputs:

   return(list(XlinearNew=XlinearNew,XgeneralNew=XgeneralNew,
               XNew=XNew,dLinear=dLinear,dGeneral=dGeneral,lenPred=lenPred))
}

############ End of predictGBnewdataProc ############
