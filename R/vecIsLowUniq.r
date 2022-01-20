########## R-function: vecIsLowUniq ##########

# Determine whether or not a the number of unique
# values in a vector is below a threshold.

# Last changed: 12 JUL 2021

vecIsLowUniq <- function(x,threshold)
   return(length(unique(x))<threshold)

########## End of vecIsLowUniq ############

