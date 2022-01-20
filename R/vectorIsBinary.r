########## R-function: vectorIsBinary ##########

# Determine whether or not a vector contains
# only binary data, defined to be that all 
# entries are either 0 or 1.

# Last changed: 21 JUN 2021

vectorIsBinary <- function(x)
   return(length(setdiff(x,c(0,1)))==0)

########## End of vectorIsBinary ############

