\name{gamselBayes.control}
\alias{gamselBayes.control}
\title{Controlling Bayesian generalized additive model selection}
\description{Function for optional use in calls to \code{gamselBayes()} to control spline basis dimension, hyperparameter choice and other specifications for generalized additive model selection via using Bayesian inference.}
\usage{
gamselBayes.control(numIntKnots = 25,truncateBasis = TRUE,numBasis = 12,
                    sigmabeta0 = 100000,sbeta = 1000,sepsilon = 1000,
                    su= 1000,rhoBeta = 0.5,rhoU = 0.5,nWarm = 1000,
                    nKept = 1000,nThin = 1,maxIter = 1000,toler = 1e-8,
                    msgCode = 1)
}
\arguments{
\item{numIntKnots}{The number of interior knots used in construction of the Demmler-Reinsch spline basis functions. The value of \code{numIntKnots} must be an integer between 8 and 50.}
\item{truncateBasis}{Boolean flag:\cr
                TRUE = truncate the Demmler-Reinsch spline basis to the value specified by \code{numBasis} (the default)\cr
                FALSE = produce a black-and-white version of the graphical check.}
\item{numBasis}{The number of spline basis functions retained after truncation when \code{truncateBasis} is \code{TRUE}. The value of \code{numBasis} cannot exceed \code{numIntKnots}+2. The default value of \code{numBasis} is 12.} 
\item{sigmabeta0}{The standard deviation hyperparameter for the Normal prior distribution on the intercept parameter for the standardized version of the data used in Bayesian inference. The default is 100000.} 
\item{sbeta}{The standard deviation hyperparameter for the Normal prior distribution on the linear component coefficients parameters for the standardized version of the data used in Bayesian inference. The default is 1000.} 
\item{sepsilon}{The scale hyperparameter for the Half-Cauchy prior distribution on the error standard deviation parameter for the standardized version of the data used in Bayesian inference. The default is 1000.}
\item{su}{The scale hyperparameter for the Half-Cauchy prior distribution on the spline basis coefficients standard deviation parameter for the standardized version of the data used in Bayesian inference. The default is 1000.}
\item{rhoBeta}{The probability parameter for the Bernoulli prior distribution on the linear component coefficients spike-and-slab auxiliary indicator variables probability parameter. The default is 0.5.}
\item{rhoU}{The probability parameter for the Bernoulli prior distribution on the spline basis coefficients spike-and-slab auxiliary indicator variables probability parameter. The default is 0.5.}
\item{nWarm}{The size of the Markov chain Monte Carlo warmup, a positive integer, when method is Markov chain Monte Carlo. The default is 1000.}
\item{nKept}{The size of the kept Markov chain Monte Carlo samples, a positive integer, when the method is Markov chain Monte Carlo. The default is 1000.}
\item{nThin}{The thinning factor for the kept Markov chain Monte Carlo samples, a positive integer, when the method is Markov chain Monte Carlo. The default is 1.}
\item{maxIter}{The maximum number of mean field variational Bayes iterations, a positive integer, when the method is mean field variational Bayes. The default is 1000.}
\item{toler}{The convergence tolerance for mean field variational Bayes iterations, a positive number less than 0.5, when the method is mean field variational Bayes. Convergence is deemed to have occurred if the relative change in the logarithm of the approximate marginal likelihood falls below \code{toler}. The default is 1e-8.}
\item{msgCode}{Code for the nature of progress messages printed to the screen. If msgCode=0 then no progress messages are printed. If msgCode=1 then a messages printed each time approximately each 10\% of the sampling is completed. If msgCode=2 then a messages printed each time approximately each 10\% of the sampling is completed. The default is 1.}
}
\value{A list containing values of each of the seventeen control parameters, packaged to supply the \code{control} argument to \code{gamselBayes}. The values for \code{gamselBayes.control} can be specified in the call to \code{gamselBayes}.
}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\references{He, V.X. and Wand, M.P. (2021). 
Generalized additive model selection via Bayesian inference. Submitted.}
\examples{library(gamselBayes) 

# Generate some simple regression-type data:

set.seed(1) ; n <- 1000 ; x1 <- rbinom(n,1,0.5)  
x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)

# Obtain a gamselBayes() fit for the data:

fit <- gamselBayes(y,Xlinear,Xgeneral)
summary(fit) ; plot(fit) ; checkChains(fit)

# Now modify some of the control values:

\donttest{fitControlled <- gamselBayes(y,Xlinear,Xgeneral,control = gamselBayes.control(
                             numIntKnots = 35,truncateBasis = FALSE,
                             sbeta = 10000,su = 10000,nWarm = 2000,nKept = 1500))
summary(fitControlled) ; plot(fitControlled) ; checkChains(fitControlled)}
}



