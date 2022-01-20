########## R-function: vectorIsConstant ##########

# Determine whether or not a vector has all entries
# equal to a particular constant value.

# Last changed: 23 JUN 2021

vectorIsConstant <- function(x)
{
   if (length(x)==1) 
      return(TRUE)

   if (length(x)>1)
      return(!any(x[1]!=x[-1]))
}

########## End of vectorIsConstant ############

