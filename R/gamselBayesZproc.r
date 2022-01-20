########## R function: gamselBayesZproc ##########

# For conducting Z matrix processing for the gamselBayes() function:

# Last changed: 07 DEC 2021

gamselBayesZproc <- function(X,sampSize,dLinear,dGeneral,numIntKnots,
                             truncateBasis,numBasis)
{
   if (dGeneral>0)
   {
      # Form the Z matrix:

      Z <- NULL
      ncZvec <- NULL
      rangexList <- vector("list",dGeneral)
      intKnotsList <- vector("list",dGeneral)
      OStoDRmatList <- vector("list",dGeneral)
      for (j in 1:dGeneral)
      {
         xCurr <- X[,dLinear+j]
         numUniqx <- length(unique(xCurr))

         if (numUniqx<50) 
         {
            numIntKnotsCurr <- round(numUniqx/3)
            numBasisCurr <- numIntKnotsCurr + 2
         }

         if (numUniqx>=50) 
         {
            numIntKnotsCurr <- numIntKnots
            numBasisCurr <- numBasis
         }

         rangexCurr <- c(1.05*min(xCurr) - 0.05*max(xCurr),1.05*max(xCurr) - 0.05*min(xCurr))
         intKnotsCurr <-  quantile(unique(xCurr),seq(0,1,length=numIntKnotsCurr+2)
                            [-c(1,numIntKnotsCurr+2)])

         Zcurr <- ZcDR(xCurr,rangexCurr,intKnotsCurr)

         rangexList[[j]] <- rangexCurr
         intKnotsList[[j]] <- intKnotsCurr
         OStoDRmatList[[j]] <- attr(Zcurr,"OStoDRmat")

         if ((truncateBasis)&(ncol(Zcurr)>=numBasis))
            Zcurr <- Zcurr[,1:numBasisCurr]
         ncZvec <- c(ncZvec,ncol(Zcurr))

         Z <- cbind(Z,Zcurr)
      }
   }
   if (dGeneral==0)
   {
      Z <- NULL
      ncZvec <- NULL
      rangexList <- NULL
      intKnotsList <- NULL 
      OStoDRmatList <- NULL
   }

   # Return the processed inputs:

   return(list(Z=Z,ncZvec=ncZvec,rangexList=rangexList,intKnotsList=intKnotsList,
               OStoDRmatList=OStoDRmatList))
}

############ End of gamselBayesZproc ############