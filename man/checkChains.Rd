\name{checkChains}
\alias{checkChains}
\title{Check Markov chain Monte Carlo samples}
\description{Facilitates a graphical check of the Markov chain Monte Carlo samples ("chains") corresponding to key quantities for the predictors selected as having an effect.}
\usage{
checkChains(fitObject,colourVersion = TRUE,paletteNum = 1)
}
\arguments{
\item{fitObject}{\code{gamselBayes()} fit object.}
\item{colourVersion}{Boolean flag:\cr
                TRUE = produce a colour version of the graphical check (the default)\cr
                FALSE = produce a black-and-white version of the graphical check.}
\item{paletteNum}{The palette of colours for the graphical display when colourVersion is TRUE. Two palettes, numbered 1 and 2, are available. The value of \code{paletteNum} specifies the palette to use. The default value is 1.}
}
\details{A graphic is produced that summarises the Markov chain Monte Carlo samples ("chains") corresponding to key quantities for the predictors selected as having an effect. If the predictor is found to have a linear effect then the chain for its coefficient is graphically checked. It the predictor is found to have a non-linear effect then the chain for the vertical slice of the penalized spline fit at the median of the predictor sample is graphically checked. The columns of the graphic are the following summaries of each chain: (1) trace (time series) plot, (2) lag-1 plot in which each chain value is plotted against its previous value and (3) sample autocorrelation function plot as produced by the R function acf(). A rudimentary graphical assessment of convergence involveschecking that the trace plots have flat-lined, rather than having any noticeable trends. If the latter occurs that a longer warmup is recommended.}
\value{Nothing is returned.}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{ 
library(gamselBayes) 

# Generate some simple regression-type data:

set.seed(1) ; n <- 1000 ; x1 <- rbinom(n,1,0.5) ; 
x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)

# Obtain a gamselBayes() fit for the data:

fitMCMC <- gamselBayes(y,Xlinear,Xgeneral)
summary (fitMCMC)

# Obtain a graphic for checking the chains:

checkChains(fitMCMC)
\donttest{
if (require("Ecdat"))
{
   # Obtain a gamselBayes() fit for data on schools in California, U.S.A.:

   Caschool$log.avginc <- log(Caschool$avginc)
   mathScore <- Caschool$mathscr
   Xgeneral <- Caschool[,c("mealpct","elpct","calwpct","compstu","log.avginc")]
   fitMCMC <- gamselBayes(y = mathScore,Xgeneral = Xgeneral)
   summary(fitMCMC)

   # Obtain a graphic for checking the chains:

   checkChains(fitMCMC)
}
}
}

