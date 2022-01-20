\name{effectTypes}
\alias{effectTypes}
\title{Tabulate the estimated effect types from a Bayesian generalized additive model object}
\description{Produces a tabular summary of the estimated effect types from a gamselBayes() fit object.}
\usage{effectTypes(fitObject)}
\arguments{
\item{fitObject}{\code{gamselBayes()} fit object.}
}
\details{Two tables are printed to standard output. The first table lists the names of the predictors that are estimated as having a linear effect. The second table lists the names of the predictors that are estimated as having a nonlinear effect.}
\value{Nothing is returned.}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{library(gamselBayes) 

# Generate some simple regression-type data:

set.seed(1) ; n <- 1000 ; x1 <- rbinom(n,1,0.5) 
x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)

# Obtain a gamselBayes() fit for the data:

fit <- gamselBayes(y,Xlinear,Xgeneral)

# Tabulate the estimated effect types:

effectTypes(fit)
}

