########## R-function: isColour ##########

# Determine whether or not an object is a 
# an R colour character string.

# Last changed: 09 JUL 2021

isColour <- function(x)
   return(length(setdiff(x,colors()))==0)

########## End of isColour ############

