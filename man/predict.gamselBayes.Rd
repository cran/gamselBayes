\name{predict.gamselBayes}
\alias{predict.gamselBayes}
\title{Obtain predictions from a \code{gamselBayes()} fit}
\description{The estimated non-linear components of the generalized additive model selected via \code{gamselBayes} are plotted.}
\usage{
\method{predict}{gamselBayes}(object,newdata,type = "response",...)
}
\arguments{
\item{object}{A \code{gamselBayes()} fit object.}
\item{newdata}{A two-component list the following components:\cr
A data frame containing new data on the predictors that are only permitted to have a linear or zero effect and, if not NULL, must have the same names as the \code{Xlinear} component of \code{object}.\cr 
A data frame containing new data for the predictors that are permitted to have a linear, nonlinear or zero effect and, if not NULL, must have the same names as the \code{Xgeneral} component of \code{object}.\cr
If both \code{Xlinear} and \code{Xgeneral} are not NULL then they must have the same numbers of rows.
}
\item{type}{A character string for specifying the type of prediction, with the following options:
"link", "response" or "terms", which leads to the value as described below.
}
\item{...}{A place-holder for other prediction parameters.}
}
\value{A vector or data frame depending on the value of \code{type}:\cr
If \code{type}="link" then the value is  a vector of linear predictor-scale fitted values.\cr
If \code{type}="response" and \code{family}="binomial" then the value is a vector of probability-scale fitted values.
Otherwise (i.e. \code{family}="binomial") the value is the vector of predictor-scale fitted values.\cr
If \code{type}="terms" then the value is a  a data frame with number of columns equal to the total number of predictors. Each column is the contribution to the vector of linear predictor-scale fitted values from each predictor. These contributions do not include the intercept predicted value. The intercept predicted value is included as an attribute of the returned data frame.
}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\examples{ 
library(gamselBayes) 

# Generate some simple regression-type data:

n <- 1000 ; x1 <- rbinom(n,1,0.5) ; x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)
names(Xlinear) <- c("x1") ; names(Xgeneral) <- c("x2","x3","x4")

# Obtain and summarise a gamselBayes() fit for the data:

fit <- gamselBayes(y,Xlinear,Xgeneral)
summary(fit)   

# Obtain some new data:

nNew <- 10
x1new <- rbinom(nNew,1,0.5) ; x2new <- runif(nNew) ; x3new <- runif(nNew) 
x4new <- runif(nNew)
XlinearNew <- data.frame(x1new) ; names(XlinearNew) <- "x1"
XgeneralNew <- data.frame(x2new,x3new,x4new)
names(XgeneralNew) <- c("x2","x3","x4")

newdataList <- list(XlinearNew,XgeneralNew)

# Obtain predictions at the new data:

predObjDefault <- predict(fit,newdata=newdataList)
print(predObjDefault)
}



