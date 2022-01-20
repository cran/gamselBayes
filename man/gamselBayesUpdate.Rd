\name{gamselBayesUpdate}
\alias{gamselBayesUpdate}
\title{Update a \code{gamselBayes()} fit object.}
\description{Facilitates updating of \code{gamselBayes} fit object when two key parameters controlling model selection are modified. Use of \code{gamselBayesUpdate()} allows for fast tweaking of such parameters without another, potentially time-consuming, call to \code{gamselBayes()}.}
\usage{
gamselBayesUpdate(fitObject,lowerMakesSparser = NULL)
}
\arguments{
\item{fitObject}{\code{gamselBayes()} fit object.}
\item{lowerMakesSparser}{A threshold parameter between 0 and 1, which is such that lower values lead to sparser fits.}
}
\details{The \code{gamselBayesUpdate()} function is applicable when a \code{gamselBayes()} fit object has been obtained for particular data inputs \code{y}, \code{Xlinear} and \code{Xgeneral} (as well as other tuning-type inputs) and the analyst is interested in changing the value of the parameter that controls model selection. This parameter is named \code{lowerMakesSparser}, and is described above. A call to \code{gamselBayesUpdate()} with a new value of \code{lowerMakesSparse} produces an updated \code{gamselBayes()} fit object with, potentially, different effect type estimates.}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\value{An object of class \code{gamselBayes} with the same components as those produced by the \code{gamselBayes()} function. See \code{help(gamselBayes)} for details.}
\examples{ 
\donttest{library(gamselBayes) 

# Generate some regression-type data:

set.seed(1) ;  n <- 5000  ; numPred <- 15
Xgeneral <- as.data.frame(matrix(runif(n*numPred),n,numPred))
names(Xgeneral) <- paste("x",1:numPred,sep="")

y <- as.vector(0.1 + 0.4*Xgeneral[,1] - 2*pnorm(3-6*Xgeneral[,2]) 
               - 0.9*Xgeneral[,4] + cos(3*pi*Xgeneral[,5]) + 2*rnorm(n))

# Obtain and assess a gamselBayes() fit:

fitOrig <- gamselBayes(y,Xgeneral = Xgeneral)
summary(fitOrig) ; plot(fitOrig)
print(fitOrig$effectTypesHat)

# Update the gamselBayes() fit object with a new value of 
# the "lowerMakesSparser" parameter:

fitUpdated <-  gamselBayesUpdate(fitOrig,lowerMakesSparser = 0.6)
summary(fitUpdated) ; plot(fitUpdated)
print(fitUpdated$effectTypesHat)
}
}


