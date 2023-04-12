\name{plot.gamselBayes}
\alias{plot.gamselBayes}
\title{Plot components of the selected generalized additive model from a \code{gamselBayes()} fit}
\description{The estimated non-linear components of the generalized additive model selected via \code{gamselBayes} are plotted.}
\usage{
\method{plot}{gamselBayes}(x,credLev = 0.95,gridSize = 251,nMC = 5000,varBand = TRUE,
     shade = TRUE,yscale = "response",rug = TRUE,rugSampSize = NULL,estCol = "darkgreen",
     varBandCol = NULL,rugCol = "dodgerblue",mfrow = NULL,xlim = NULL,ylim = NULL,
     xlab = NULL,ylab = NULL,mai = NULL,pages = NULL,cex.axis = 1.5,cex.lab = 1.5,...)
}
\arguments{
\item{x}{A \code{gamselBayes()} fit object.}
\item{credLev}{A number between 0 and 1 such that the credible interval band has (100*credLev)\% approximate pointwise coverage. The default value is 0.95.}
\item{gridSize}{A number of grid points used to display the density estimate curve and the pointwise credible interval band. The default value is 251.}
\item{nMC}{The size of the Monte Carlo sample, a positive integer, for carrying out approximate inference from the mean field variational Bayes-approximate posterior distributions when the method is mean field variational Bayes. The default value is 5000.}
\item{varBand}{Boolean flag specifying whether or not a variability band is included:\cr
                TRUE = add a pointwise approximate (100*credLev)\% credible set variability band (the default)\cr
                FALSE = only plot the density estimate, without a variability band.}
\item{shade}{Boolean flag specifying whether or not the variability band is displayed using shading:\cr
                TRUE = display the variability band using shading (the default)\cr
                FALSE = display the variability band using dashed curves.}
\item{yscale}{Character string specifying the vertical axis scale for display of estimated non-linear functions:\cr
                "response" = display on the response scale)\cr
                "link" = display on the link scale.}
\item{rug}{Boolean flag specifying whether or not rug-type displays for predictor data are used:\cr
                TRUE = show the predictor data using rug-type displays (the default)\cr
                FALSE = do not show the predictor data.}
\item{rugSampSize}{The size of the random sample sample of each predictor to be used in rug-type displays.}
\item{estCol}{Colour of the density estimate curve. The default value is "darkgreen".}
\item{varBandCol}{Colour of the pointwise credible interval variability band. If \code{shade}=TRUE then the default value is "palegreen". If \code{shade}=FALSE then the default value is \code{estCol}.}
\item{rugCol}{Colour of rug plot that shows values of the predictor data. The default value is "dodgerblue".}
\item{mfrow}{An optional two-entry vector for specifying the layout of the nonlinear fit displays.}
\item{xlim}{An optional two-column matrix for specification of horizontal frame limits in the plotting of the estimated non-linear predictor effects. The number of rows in \code{xlim} must equal length(effectTypesVector(x)=="nonlin"). If any of the rorows of \code{xlim} contain \code{NA} then the default horizontal frame limits value for that predictor is used. Therefore, if there are several predictors selected as non-linear and horizontal frame adjustments are required for a few of them then then \code{xlim} can be a matrix with mainly \code{NA} values, and with non-\code{NA} frame specification limits in the relevant rows.} 
\item{ylim}{The same as \code{xlim}, except for vertical frame limits.}
\item{xlab}{An optional vector of character strings containing labels for the horizontal axes. The number of entries in \code{xlab} must equal length(effectTypesVector(x)=="nonlin").}
\item{ylab}{The same as \code{xlab}, except for vertical axis labels.}
\item{mai}{An optional numerical vector of length 4 for specification of inner margin dimensions of each panel, ordered clockwise from below the panel.}
\item{pages}{An optional positive integer that specifies the number of pages used to display the non-linear function estimates.}
\item{cex.axis}{An optional positive scalar value for specification of the character expansion factor for the axis values.}
\item{cex.lab}{An optional positive scalar value for specification of the character expansion factor for the labels.}
\item{...}{Place-holder for other graphical parameters.}
}
\details{The estimated non-linear components of the selected generalized additive model are plotted. Each plot corresponds to a slice of the selected generalized additive model surface with all other predictors set to their median values. Pointwise credible intervals unless \code{varBand} is \code{FALSE}.}
\value{Nothing is returned.}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{library(gamselBayes) 

# Generate some simple regression-type data:

set.seed(1) ; n <- 1000 ; x1 <- rbinom(n,1,0.5) ; 
x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)

# Obtain a gamselBayes() fit for the data and print out a summary:

fit <- gamselBayes(y,Xlinear,Xgeneral)
summary(fit)

# Plot the predictor effect(s) estimated as being non-linear:

plot(fit)

# Plot the same fit(s) but with different colours and style:

plot(fit,shade = FALSE,estCol = "darkmagenta",varBandCol = "plum",
     rugCol = "goldenrod")

\donttest{if (require("Ecdat"))
{
   # Obtain a gamselBayes() fit for data on schools in California, U.S.A.:

   Caschool$log.avginc <- log(Caschool$avginc)
   mathScore <- Caschool$mathscr
   Xgeneral <- Caschool[,c("mealpct","elpct","calwpct","compstu","log.avginc")]
   fit <- gamselBayes(y = mathScore,Xgeneral = Xgeneral)
   summary(fit)

   # Plot the predictor effect(s) estimated as being non-linear:

   plot(fit)
}}
}



