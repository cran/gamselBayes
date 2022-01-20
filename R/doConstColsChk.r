########## R-function: doConstColsChk ##########

# Do a check of whether any constants of a matrix,
# for situations where such a matrix is supposed to
# have all columns being non-constant, and report 
# violations.

# Last changed: 19 JAN 2022

doConstColsChk <- function(A,Aname)
{
  A <- as.matrix(A)
  whichColsConstant <- apply(A,2,vectorIsConstant)
  if (any(whichColsConstant)) 
  {
     constColNos <- (1:ncol(A))[whichColsConstant]
     for (j in 1:length(constColNos))
        message("The entries of column number ",constColNos[j]," of ",Aname," are identical.\n",sep="")
     stop(paste(Aname," is illegal due to having columns with all entries identical.\n",sep=""))
  }
  invisible()
}

########## End of doConstColsChk ############

