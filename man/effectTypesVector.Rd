\name{effectTypesVector}
\alias{effectTypesVector}
\title{Obtain the estimated effect types from a Bayesian generalized additive model object}
\description{Extracts the vector of estimated effect types from a gamselBayes() fit object.}
\usage{effectTypesVector(fitObject)}
\arguments{
\item{fitObject}{\code{gamselBayes()} fit object.}
}
\value{A vector of character strings having the same length as the number of predictors, which conveys the estimated effect types.}
\details{The result is a vector of character strings having the same length as the total number of predictors inputted through \code{Xlinear} and \code{Xgeneral}. The character strings are one of "linear", "nonlinear" and "zero" according to whether each predictor is estimated as having a linear effect, nonlinear effect of zero effect. The ordering in the returned vector matches that of the columns of \code{Xlinear} and then the columns of \code{Xgeneral}.}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{library(gamselBayes) 

# Generate some simple regression-type data:

set.seed(1) ; n <- 1000 ; x1 <- rbinom(n,1,0.5) ; 
x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)

# Obtain a gamselBayes() fit for the data:

fit <- gamselBayes(y,Xlinear,Xgeneral)

# Obtain the vector of effect types:

effectTypesVector(fit)
}

