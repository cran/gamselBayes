########## R function: plotGBinputProc ##########

# For conducting input processing for the plot.gamselBayes() 
# function:

# Last changed: 01 NOV 2021

plotGBinputProc <- function(credLev,gridSize,nMC,varBand,yscale,shade,rug,rugSampSize,
                            estCol,varBandCol,rugCol,mai,pages,sampSizeVal)
{
   # Check the legality of the inputted "credLev" value:

   if (credLev<=0) 
   {
      warnStr1 <- "The inputted credible level is zero or negative."
      warnStr2 <- "The default value of 0.95 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      credLev <- 0.95
   }
   if (credLev>=1) 
   {
      warnStr1 <- "The inputted credible level is 1 or higher."
      warnStr2 <- "The default value of 0.95 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      credLev <- 0.95
   }

   # Check the legality of the inputted "gridSize" value:

   gridSize <- round(gridSize)
   if ((gridSize<11)|(gridSize>100000))
   {
      warnStr1 <- "The inputted plotting grid size is too low or too high."
      warnStr2 <- "The default value of 251 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      gridSize <- 251
   }

   # Check the legality of the inputted "nMC" value:

   nMC <- round(nMC)
   if (nMC<51)
   {
      warnStr1 <- "The inputted number of Monte Carlo samples for variational."
      warnStr2 <- "inference is too low. The default value of 1000 was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      nMC <- 1000
   }

   # Check the legality of the inputted "varBand" value:

   if (!(any(varBand==c(FALSE,TRUE))))
   {
      warnStr1 <- "The inputted flag for inclusion of a variability band is neither"
      warnStr2 <- "FALSE nor TRUE. The default value of TRUE was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      varBand <- TRUE
   }

   # Check the legality of the inputted "yscale" value:

   if (!(any(yscale==c("link","response"))))
   {
      warnStr1 <- "The inputted flag for the vertical scale used in display of nonlinear effects is"
      warnStr2 <- "neither \"link\" nor \"response\". The default value of \"response\" was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      yscale <- "response"
   }

   # Check the legality of the inputted "shade" value:

   if (!(any(shade==c(FALSE,TRUE))))
   {
      warnStr1 <- "The inputted flag for use of shading for variability band display is"
      warnStr2 <- "neither FALSE nor TRUE. The default value of TRUE was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      shade <- TRUE
   }

   # Check the legality of the inputted "estCol" value:

   if (!isColour(estCol))
   {
      warnStr1 <- "The inputted estimated function colour is not a"
      warnStr2 <- "character string matching an element of colors()."
      warnStr3 <- "The default value was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
             immediate.=TRUE)
      estCol <- "darkgreen"
   }

   # Check the legality and set default of the inputted 
   # "varBandCol" value:

   if (!is.null(varBandCol))
   {
      if (!isColour(varBandCol))
      {
         warnStr1 <- "The inputted variability band colour is not a"
         warnStr2 <- "character string matching an element of colors()."
         warnStr3 <- "The default value was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
                 immediate.=TRUE)
         if (!shade) varBandCol <- estCol
         if (shade) varBandCol <- "palegreen"
      }
   }
   if (is.null(varBandCol))
   {
      if (!shade) varBandCol <- estCol
      if (shade) varBandCol <- "palegreen"
   }

   # Check the legality of the inputted "rug" value:

   if (!(any(rug==c(FALSE,TRUE))))
   {
      warnStr1 <- "The inputted flag for rug-type display of predictor data is"
      warnStr2 <- "neither FALSE nor TRUE. The default value of TRUE was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",sep=""),immediate.=TRUE)
      rug <- TRUE
   }

   # Check the legality of the inputted "rugSampSize" value:

   if (!is.null(rugSampSize))
   {
      if ((rugSampSize<1)|(rugSampSize>sampSizeVal))
      {
         warnStr1 <- "The inputted sample size to be used in rug-type display of"
         warnStr2 <- "predictor data is less than 1 or higher than the sample size."
         warnStr3 <- "The default value was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),immediate.=TRUE)
         rugSampSize <- NULL
      }
   } 

   # Check the legality of the inputted "rugCol" value:

   if (!isColour(rugCol))
   {
      warnStr1 <- "The inputted rug colour is not a character"
      warnStr2 <- "string matching an element of colors()."
      warnStr3 <- "The default value was used instead."
      warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
             immediate.=TRUE)
      rugCol <- "dodgerblue"
   }

   # Check the legality of the inputted "mai" value:

   if (!is.null(mai))
   {
      if ((!is.numeric(mai))|(length(mai)!=4)|(any(mai<0)))
      {
         warnStr1 <- "The inputted inner margin value is not a vector"
         warnStr2 <- "of length 4 with all entries non-negative."
         warnStr3 <- "The default value was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
             immediate.=TRUE)
         mai <- NULL
      }
   }

   # Check the legality of the inputted "pages" value:

   if (!is.null(pages))
   {
      if ((!is.numeric(pages))|(pages<1))
      {
         pages <- round(pages)
         warnStr1 <- "The inputted number of pages for displaying non-linear"
         warnStr2 <- "function estimates is not a positive integer."
         warnStr3 <- "The default value was used instead."
         warning(paste(warnStr1,"\n",warnStr2,"\n",warnStr3,"\n",sep=""),
                immediate.=TRUE)
         pages <- NULL
      }
   }

   # Return the processed inputs:

   return(list(credLev=credLev,gridSize=gridSize,nMC=nMC,varBand=varBand,
               yscale=yscale,shade=shade,rug=rug,rugSampSize=rugSampSize,
               estCol=estCol,varBandCol=varBandCol,rugCol=rugCol,mai=mai,
               pages=pages))
}

############ End of plotGBinputProc ############