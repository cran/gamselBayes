########## R-function: standardiseDataFrame ##########

# Standardise a matrix, which entails linearly transforming 
# the matrix so that each of its columns have zero mean
# and unit standard deviation.

# Last changed: 17 SEP 2021

standardiseDataFrame <- function(A)
{
   # Ensure that A is a data frame:

   if (!is.data.frame(A))
      A <- as.data.frame(A)

   # Extract the names of A:

   Anames <- names(A)

   # Standardise each column of A:

   A <- scale(A)
   dimA <- dim(A)
   attributes(A) <- NULL
   A <- as.data.frame(matrix(A,dimA[1],dimA[2]))

   # Ensure that the column names of A are retained:

   names(A) <- Anames

   return(A)
}

########## End of standardiseDataFrame ############

