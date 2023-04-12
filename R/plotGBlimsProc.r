########## R function: plotGBlimsProc ##########

# For conducting limits processing for the plot.gamselBayes() 
# function:

# Last changed: 21 JUL 2021

plotGBlimsProc <- function(xlim,ylim,numPanels)
{
   if (numPanels>0)
   {
      if (!is.null(xlim))  # Check the legality of "xlim":
      {
         xlim <- as.matrix(xlim)

         # Make sure that "xlim" has the correct number of columns:

         if (ncol(xlim)!=2)
            stop("If \"xlim\" is non-NULL then it must be a two-column matrix.\n")

         # Make sure that "xlim" has the correct number of rows:

         if (nrow(xlim)!=numPanels)
         {
            stopStr1 <- "If \"xlim\" is non-NULL then it must be a two-column matrix"
            stopStr2 <- "with the number of rows equalling the number of predictors"
            stopStr3 <- "estimated as having a non-linear effect. For the current fit"
            stopStr4 <- paste("object this number is ",numPanels,".",sep="")
            stop(paste(stopStr1,"\n",stopStr2,"\n",stopStr3,"\n",stopStr4,"\n",sep=""))
         }

         # Make sure that each row of "xlim" is ordered:

         for (iPanel in 1:numPanels)
         {
            if ((!is.na(xlim[iPanel,1]))&(!is.na(xlim[iPanel,2])))
               if (xlim[iPanel,2]<xlim[iPanel,1])
                  stop(paste("Row number ",iPanel," of \"xlim\" has unordered limits.\n",sep=""))
         }
      }

      if (!is.null(ylim))  # Check the legality of "ylim":
      {
         ylim <- as.matrix(ylim)

         # Make sure that "ylim" has the correct number of columns:

         if (ncol(ylim)!=2)
            stop("If \"ylim\" is non-NULL then it must be a two-column matrix.\n")

         # Make sure that "ylim" has the correct number of rows:

         if (nrow(ylim)!=numPanels)
         {
            stopStr1 <- "If \"ylim\" is non-NULL then it must be a two-column matrix"
            stopStr2 <- "with the number of rows equalling the number of predictors"
            stopStr3 <- "estimated as having a non-linear effect. For the current fit"
            stopStr4 <- paste("object this number is ",numPanels,".",sep="")
            stop(paste(stopStr1,"\n",stopStr2,"\n",stopStr3,"\n",stopStr4,"\n",sep=""))
         }

         # Make sure that each row of "ylim" is ordered:

         for (iPanel in 1:numPanels)
         {
            if ((!is.na(ylim[iPanel,1]))&(!is.na(ylim[iPanel,2])))
               if (ylim[iPanel,2]<ylim[iPanel,1])
                   stop(paste("Row number ",iPanel," of \"ylim\" has unordered limits.\n",sep=""))
         }
      }

   }

   # Return the processed inputs:

   return(list(xlim=xlim,ylim=ylim))
}

############ End of plotGBlimsProc ############