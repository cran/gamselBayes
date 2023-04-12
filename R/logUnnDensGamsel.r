########## R function: logUnnDensGamsel #########

# For carrying out slice sampling for Bayesian
# versions of gamsel-type models.

# Last changed: 27 JAN 2021

logUnnDensGamsel <- function(x,parm1,parm2,parm3,parm4,bFun)
  return(parm1*x - parm2*x*x - sum(bFun(parm3 + parm4*x))) 

############ End of logUnnDensGamsel ###########

