\name{summary.gamselBayes}
\alias{summary.gamselBayes}
\title{Summarise components of the selected generalized additive model from a \code{gamselBayes()} fit}
\description{Inference summaries of the estimated linear component coefficients of the generalized additive model selected via \code{gamselBayes} are tabulated.}
\usage{
\method{summary}{gamselBayes}(object,credLev = 0.95,sigFigs = 5,nMC = 10000,...)
}
\arguments{
\item{object}{A \code{gamselBayes()} fit object.}
\item{credLev}{A number between 0 and 1 such that the credible interval band has (100*credLev)\% approximate pointwise coverage. The default value is 0.95.}
\item{sigFigs}{The number of significant figures used for the entries of the summary table.}
\item{nMC}{The size of the Monte Carlo sample, a positive integer, for carrying out approximate inference from the mean field variational Bayes-approximate posterior distributions when the method is mean field variational Bayes. The default value is 10000.}
\item{...}{Place-holder for other summary parameters.}
}
\details{If the selected generalized additive model has at least one predictor having a linear effect then a data frame is returned. The columns of the data correspond to posterior means and credible interval limits of the linear effects coefficients.}
\value{A data frame containing linear effect Bayesian inferential summaries.} 
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

# Print the summary with different values of some of the arguments:

summary(fit,credLev=0.99,sigFigs=3)

\donttest{if (require("Ecdat"))
{
   # Obtain a gamselBayes() fit for data on schools in California, U.S.A.:

   Caschool$log.avginc <- log(Caschool$avginc)
   mathScore <- Caschool$mathscr
   Xgeneral <- Caschool[,c("mealpct","elpct","calwpct","compstu","log.avginc")]
   fit <- gamselBayes(y = mathScore,Xgeneral = Xgeneral)
   summary(fit)
}}
}




