########## R-function: summChainsGamsel ##########

# Summarises chains for gamselBayes() MCMC-fit objects.

# Last changed: 06 JUL 2021

summChainsGamsel <- function(xMat,parNames,colourVersion,paletteNum)
{
   # Set the number of samples per view value:
 
   numSampsPerView <- 6

   # Define required functions:

   empty.panel <- function()
   {
      plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",
           yaxt="n",xlab="",ylab="",bty="o")
      invisible()
   }

   # Set dimension and other variables:

   numPar <- ncol(xMat)
   sampSize <- nrow(xMat)  

   if (colourVersion)
   {
      if (paletteNum==1)
         columnCols <- c("purple4","tomato","mediumblue","olivedrab4")
      
      if (paletteNum==2)
         columnCols <- c("darkmagenta", "green4","darkorange","dodgerblue")
   }
   if (!colourVersion)
   {
      columnCols <- rep(NA,numPar)
      for (i in 1:4)
         columnCols[i] <- "black"
   }

   oldpar <- par(no.readonly = TRUE)
   on.exit(par(oldpar))

   par(mfrow=c((numSampsPerView+1),4))

   par(ann=F,mar=rep(0,4),xaxt="n",yaxt="n",xpd=TRUE)
   empty.panel()
   text(0.5,0.5,"predictor",cex=2.9,col=columnCols[1])

   par(ann=F,mar=rep(0,4),xaxt="n",yaxt="n",xpd=TRUE)
   empty.panel()
   text(0.5,0.5,"trace",cex=2.9,col=columnCols[2])

   par(ann=F,mar=rep(0,4),xaxt="n",yaxt="n",xpd=TRUE)
   empty.panel()
   text(0.5,0.5,"lag 1",cex=2.9,col=columnCols[3])

   par(ann=F,mar=rep(0,4),xaxt="n",yaxt="n",xpd=TRUE)
   empty.panel()
   text(0.5,0.5,"acf",cex=2.9,col=columnCols[4])

   for (j in 1:numPar)
   {
      # Write the variable name:

      par(ann=F,mar=c(0,0,0,0),xaxt="n",yaxt="n")

      empty.panel()
      if (length(parNames[[j]])==1)
         text(0.5,0.5,parNames[[j]][1],cex=3.0,col=columnCols[1])
      if (length(parNames[[j]])==2)
      {
         text(0.5,0.7,parNames[[j]][1],cex=1.9,col=columnCols[1])
         text(0.5,0.3,parNames[[j]][2],cex=1.9,col=columnCols[1])
      }
      if (length(parNames[[j]])==3)
      {
         text(0.5,0.8,parNames[[j]][1],cex=1.55,col=columnCols[1])
         text(0.5,0.5,parNames[[j]][2],cex=1.55,col=columnCols[1])
         text(0.5,0.2,parNames[[j]][3],cex=1.55,col=columnCols[1])
      }

      # Do the trace plot:

      plot(xMat[,j],xlab="",ylab="",type="l",col=columnCols[2])
      
      # Do the lag 1 plot:

      plot(xMat[1:(sampSize-1),j],
           xMat[2:sampSize,j],xlab="",ylab="",type="n")

      points(xMat[1:(sampSize-1),j],
             xMat[2:sampSize,j],pch=1,cex=0.5,col=columnCols[3])
      
      # Do the autocorrelation function plot:

      ci.col.val <- "black"
      if (colourVersion) ci.col.val <- "blue"
      acf(xMat[,j],lag.max=20,col=columnCols[4],
          lwd=2,ci.col=ci.col.val)
      
   }

   invisible()
}

########## End of summChainsGamsel ##########
