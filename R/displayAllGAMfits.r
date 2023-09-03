########## R script: displayAllGAMfits ##########

# For displaying all fits of a generalized additive model
# for a given list of plotting vectors. The displaying
# involves having a maximum of 20 panels per screen.

# Last changed: 05 JAN 2022

displayAllGAMfits <- function(X,xgList,estgList,lowgList,uppgList,xLabsVec,
                              yLabsVec,colVec,varBandColVec,shade,shadeColVec,
                              rug,rugSampSize,rugCol,mfrow,xlimMat,ylimMat,mai,
                              pages,cex.axis,cex.lab)
{
   # Ensure that the par() setting are returned to their default values:

   oldpar <- par(no.readonly = TRUE)
   on.exit(par(oldpar))

   # Set the number of predictors:
  
   numPred <- ncol(X)

   # If "pages" is NULL then set its default value:  

   if (is.null(pages))
      pages <- ceiling(numPred/20)
  
   # Make sure that the number of pages does not exceed 
   # the number of predictors:

   if (pages>numPred)
      pages <- numPred

   # If the number of pages is too low then increase and issue a warning:

   maxPansPerPage <- ceiling(numPred/pages)
    
   if (maxPansPerPage>20)
   {
      pages <- ceiling(numPred/20)
      warnStr1 <- "The inputted number of pages is too low relative"
      warnStr2 <- "to the number of estimated non-linear functions."
      warnStr3 <- paste("The number of pages was increased to ",pages,".\n",sep="")
      warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
              immediate.=TRUE)
      maxPansPerPage <- ceiling(numPred/pages)
   }
      
   # Set the number of panels vector:

   numPanelsVec <- c(rep(maxPansPerPage ,pages-1),numPred-maxPansPerPage *(pages-1))

   # Determine the default value of "ylim":

   ylimDefault <- range(c(unlist(lowgList),unlist(uppgList)))

   # Loop through the views:

   icntLab <- 1

   iPred <- 0
   for (iPage in 1:pages)
   {
      numPanelsCurr <- numPanelsVec[iPage]
 
      if (numPanelsCurr==1) dimVec <- c(1,1)
      if (numPanelsCurr==2) dimVec <- c(1,2)

      if (any(numPanelsCurr==(3:4))) dimVec <- c(2,2)
      if (any(numPanelsCurr==(5:6))) dimVec <- c(3,2)
      if (any(numPanelsCurr==(7:9))) dimVec <- c(3,3)
      if (any(numPanelsCurr==(10:12))) dimVec <- c(4,3)
      if (any(numPanelsCurr==(13:16))) dimVec <- c(4,4)
      if (any(numPanelsCurr==(17:20))) dimVec <- c(5,4)

      # Set the "mai" values:

      if (is.null(mai))
      {
         if (all(dimVec==c(1,1))) maiVec <- c(0.85,0.85,0.3,0.1)
         if (all(dimVec==c(1,2))) maiVec <- c(0.85,0.85,0.3,0.1)
         if (all(dimVec==c(2,2))) maiVec <- c(0.65,0.35,0.05,0.1)
         if (all(dimVec==c(3,2))) maiVec <- c(0.62,0.35,0.05,0.1)
         if (all(dimVec==c(3,3))) maiVec <- c(0.55,0.35,0.05,0.1)
         if (all(dimVec==c(4,3))) maiVec <- c(0.55,0.35,0.05,0.1)
         if (all(dimVec==c(4,4))) maiVec <- c(0.55,0.35,0.05,0.1)
         if (all(dimVec==c(5,4))) maiVec <- c(0.55,0.35,0.05,0.1)
      } 
      if (!is.null(mai))
      {
          if (!is.numeric(mai)) stop("mai must be a numeric vector.\n")
          if (length(mai)!=4) stop("mai must be a vector of length 4.\n")
          if (any(mai<0)) stop("mai contains negative entries.\n")
          maiVec <- mai
      }

      if (!is.null(mfrow)) mfrowToUse <- mfrow
      if (is.null(mfrow)) mfrowToUse <- dimVec

      par(mfrow=mfrowToUse,mai=maiVec)
      for (iPlot in 1:numPanelsCurr)
      {  
         iPred <- iPred + 1

         # Set the "xlim" value for the current panel:

         xlimVal <- range(xgList[[iPred]])
         if (!is.null(xlimMat))  # Check if "xlimVal" to be overwritten by 
         {                       # a user-specified value:
            if ((!is.na(xlimMat[iPred,1]))&(!is.na(xlimMat[iPred,2])))
               xlimVal <- xlimMat[iPred,]
         }

         # Set the "ylim" value for the current panel:

         ylimVal <- ylimDefault
         if (!is.null(ylimMat))  # Check if "ylimVal" to be overwritten by 
         {                       # a user-specified value:
            if ((!is.na(ylimMat[iPred,1]))&(!is.na(ylimMat[iPred,2])))
               ylimVal <- ylimMat[iPred,]
         }

         # Extract current colour values:

         shadeColCurr <- shadeColVec[iPred]
         estColCurr <- colVec[iPred]

         # Extract current plotting vectors:

         xgCurr <- xgList[[iPred]]
         estgCurr <- estgList[[iPred]]
         lowgCurr <- lowgList[[iPred]]
         uppgCurr <- uppgList[[iPred]]

         # Extract current rug data vector:

         rugData <- X[,iPred]
       
         # Possibly carry out horizontal truncation of the plotting
         # vectors and rug data vector:

         indsToKeep <- (1:length(xgCurr))[(xgCurr>=xlimVal[1])&(xgCurr<=xlimVal[2])]
         xgCurr <- xgCurr[indsToKeep]
         estgCurr <- estgCurr[indsToKeep]
         lowgCurr <- lowgCurr[indsToKeep]
         uppgCurr <- uppgCurr[indsToKeep]

         indsToKeep <- (1:length(rugData))[(rugData>=xlimVal[1])&(rugData<=xlimVal[2])]
         rugData <- rugData[indsToKeep]

         # Possibly replace "rugData" with a random sub-sample:

         if ((length(rugData)>0)&(!is.null(rugSampSize)))
         {
            if (rugSampSize<length(rugData))
            {   
               subInds <- sample(1:length(rugData),rugSampSize,replace=FALSE)
               rugData <- rugData[subInds]
            }
         }

         # Create the plot for the current panel:

         plot(0,type="n",bty="l",xlab=xLabsVec[icntLab],ylab=yLabsVec[icntLab],xlim=xlimVal,ylim=ylimVal,
              cex.axis = cex.axis,cex.lab = cex.lab)
         icntLab <- icntLab + 1

         if (length(xgCurr)>0)
         {
            if (shade)
               polygon(c(xgCurr,rev(xgCurr)),c(lowgCurr,rev(uppgCurr)),
                       col=shadeColCurr,border=FALSE)

            if (!shade)
            {
               lines(xgCurr,lowgCurr,col=varBandColVec[iPage],lwd=2,lty=2)
               lines(xgCurr,uppgCurr,col=varBandColVec[iPage],lwd=2,lty=2)
            } 

            lines(xgCurr,estgCurr,col=estColCurr,lwd=2)
         }
         if (rug&(length(rugData)>0)) rug(rugData,quiet=TRUE,col=rugCol)
      }
   
      if (iPage<pages)
         readline("Hit Enter to continue.\n")
   }
}

############ End of displayAllGAMfits ############
