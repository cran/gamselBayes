########## R function: gBsuffStats ##########

# For computing sufficient statistic quantities for the 
# gamselBayes() function:

# Last changed: 08 OCT 2021

gBsuffStats <- function(y,X,Z)
{
   # Compute the sufficient statistic quantities:

   XTy <- as.vector(crossprod(X,y))
   XTX <- crossprod(X)

   if (!is.null(Z))
   {
      ZTy <- as.vector(crossprod(Z,y))
      ZTX <- crossprod(Z,X)
      ZTZ <- crossprod(Z)
   }
   if (is.null(Z))
   {
      ZTy <- NA
      ZTX <- NA
      ZTZ <- NA
   }

   # Return the sufficient statistic quantities:

   return(list(XTy=XTy,XTX=XTX,ZTy=ZTy,ZTX=ZTX,ZTZ=ZTZ))
}

############ End of gBsuffStats ############
