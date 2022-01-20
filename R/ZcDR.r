########## R-function: ZcDR ##########

# For creation of canonical Demmler-Reinsch-type Z matrices.

# Last changed: 05 OCT 2021

ZcDR <- function(x,range.x,intKnots)
{
   # Set defaults for `range.x' and `intKnots'

   if (missing(range.x))
      range.x <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))

   if (missing(intKnots))
   {
      numIntKnots <- min(length(unique(x)),25)
      intKnots <- quantile(unique(x),seq(0,1,length=
                  (numIntKnots+2))[-c(1,(numIntKnots+2))])
   }

   # Obtain the design matrix containing canonical O'Sullivan spline basis functions:

   ZOS <- ZOSull(x,range.x,intKnots)
  
   # Obtain the corresponding Demmler-Reinsch basis:

   Cmat <- cbind(1,x,ZOS)
   Dmat <- diag(c(0,0,rep(1,ncol(ZOS))))
   svdC <- svd(Cmat)
   UC <- svdC$u
   VC <- svdC$v
   dC <- svdC$d
   svdD <- svd(t(t(crossprod(t(crossprod(VC,Dmat)),VC)/dC)/dC))
   UD <- svdD$u
   dD <- svdD$d
   CDR <- crossprod(t(UC),svdD$u) 

   # Do the damping and column ordering adjusment and then obtain the 
   # final "Z" matrix containing canonical Demmler-Reinsch basis functions.

   ncCDR <- ncol(CDR)
   dampFacVec <- sqrt(dD[ncCDR-2]/dD)
   dampFacVec[c(ncCDR-1,ncCDR)] <- rep(1,2)
   CcDR <- t(t(CDR)*dampFacVec)
   dampFacVecAdj <- dampFacVec
   dampFacVecAdj[c(ncCDR-1,ncCDR)] <- rep(1.1,2)
   CcDR <- CcDR[,rev(order(dampFacVecAdj))]

   Z <- CcDR[,-c(1,2)]

   # Obtain the O'Sullivan to Demmler-Reinsch transformation matrix:

   OStoDRmat <- t(crossprod(UD,t(VC)/dC)*dampFacVec)
   OStoDRmat <- OStoDRmat[,rev(order(dampFacVecAdj))]

   # Add the `range.x' and 'intKnots' as attributes
   # of the return object:

   attr(Z,"range.x") <- range.x
   attr(Z,"intKnots") <- intKnots
   attr(Z,"OStoDRmat") <- OStoDRmat

   # Return attribute-adorned "Z" matrix:

   return(Z)
}

############ End of ZcDR ############
