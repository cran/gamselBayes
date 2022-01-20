\name{gamselBayes}
\alias{gamselBayes}
\title{Bayesian generalized additive model selection including a fast variational option}
\description{Selection of predictors and the nature of their impact on the mean response (linear versus non-linear) is a fundamental problem in regression analysis. This function uses the generalized additive models framework for estimating predictors effects. An approximate Bayesian inference approach and has two options for achieving this: (1) Markov chain Monte Carlo and (2) mean field variational Bayes.}
\usage{
gamselBayes(y,Xlinear  = NULL,Xgeneral = NULL,method = "MCMC",lowerMakesSparser = NULL,   
            family = "gaussian",verbose = TRUE,control = gamselBayes.control())
}
\arguments{
\item{y}{Vector containing the response data. If 'family = "gaussian"' then the response data are modelled as being continuous with a Gaussian distribution. If 'family = "binomial"' then the response data must be binary with 0/1 coding.}
\item{Xlinear}{Data frame with number of rows equal to the length of \code{y}. Each column contains data for a predictor which potentially has a linear or zero effect, but not a nonlinear effect. Binary predictors must be inputted through this matrix.}
\item{Xgeneral}{A data frame with number of rows equal to the length of \code{y}. Each column contains data for a predictor which potentially has a linear, nonlinear or zero effect. Binary predictors cannot be inputted through this matrix.}
\item{method}{Character string for specifying the method to be used:\cr
  "MCMC" = Markov chain Monte Carlo,\cr
  "MFVB" = mean field variational Bayes.}
\item{lowerMakesSparser}{A threshold parameter between 0 and 1, which is such that lower values lead to sparser fits.}
\item{family}{Character string for specifying the response family:\cr
  "gaussian" = response assumed to be Gaussian with constant variance,\cr
  "binomial" = response assumed to be binary.\cr
The default is "gaussian".}
\item{verbose}{Boolean variable for specifying whether or not progress messages are printed to the console. The default is TRUE.}
\item{control}{Function for controlling the spline bases, Markov chain Monte Carlo sampling, mean field variational Bayes and other specifications.}
}
\details{Generalized additive model selection via approximate Bayesian inference is provided. Bayesian mixed model-based penalized splines with spike-and-slab-type coefficient prior distributions are used to facilitate fitting and selection. The approximate Bayesian inference engine options are: (1) Markov chain Monte Carlo and (2) mean field variational Bayes. Markov chain Monte Carlo has better Bayesian inferential accuracy, but requires a longer run-time. Mean field variational Bayes is faster, but less accurate. The methodology is described in He and Wand (2021) <arXiv:2201.00412>. 
}
\value{An object of class \code{gamselBayes}, which is a list with the following components:
\item{method}{the value of \code{method}.}
\item{family}{the value of \code{family}.}
 \item{Xlinear}{the inputted design matrix containing predictors that can only have linear effects.}  
\item{Xgeneral}{the inputted design matrix containing predictors that are potentially have non-linear effects.} 
\item{rangex}{the value of the control parameter \code{rangex}.} 
\item{intKnots}{the value of the control parameter \code{intKnots}.} 
\item{truncateBasis}{the value of the control parameter \code{truncateBasis}.} 
\item{numBasis}{the value of the control parameter \code{numBasis}.} 
\item{MCMC}{a list such that each component is the retained Markov chain Monte Carlo (MCMC)sample for a model parameter. The components are:\cr
   \code{beta0} = overall intercept.\cr
   \code{betaTilde} = linear component coefficients without multiplication by the \code{gammaBeta} values.\cr
   \code{gammaBeta} = linear component coefficients spike-and-slab auxiliary indicator variables.\cr
   \code{sigmaBeta} = standard deviation of the linear component coefficients.\cr
   \code{rhoBeta} = the Bernoulli distribution probability parameter of the linear component coefficients spike-and-slab auxiliary indicator variables.\cr
   \code{uTilde} = spline basis function coefficients without multiplication by the \code{gammaUMCMC} values. The MCMC samples are stored in a list. Each list component corresponds to a predictor that is treated as potentially having a non-linear effect, and is a matrix with columns corresponding to the spline basis function coefficients for that predictor and rows corresponding to the retained MCMC samples.\cr
   \code{gammaU} = spline basis coefficients spike-and-slab auxiliary indicator variables.  The MCMC samples are stored in a list. Each list component corresponds to a predictor that is treated as potentially having a non-linear effect, and is a matrix with columns corresponding to the spline basis function coefficients for that predictor and rows corresponding to the retained MCMC samples.\cr
   \code{rhoU} = the Bernoulli distribution probability parameters of the spline basis component coefficients spike-and-slab auxiliary indicator variables. The MCMC samples are stored in a matrix. Each column corresponds to a predictor that is treated as potentially having a non-linear effect. The rows of the matrix correspond to the retained MCMC samples.\cr
  \code{sigmaEps} = error standard deviation.
} 
\item{MFVB}{a list such that each component is the mean field variational Bayes approximate posterior density function, or q-density, parameters. The components are:\cr 
   \code{beta0} = a vector with 2 entries, consisting of the mean and variance of the Univariate Normal q-density of the overall intercept.\cr
   \code{betaTilde} = a two-component list containing the Multivariate Normal q-density parameters of linear component coefficients without multiplication by the means of the \code{gammaBeta} q-densities. The list components are: \code{mu.q.betaTilde}, the mean vector; \code{Sigma.q.betaTilde}, the covariance matrix.\cr
   \code{gammaBeta} = a vector containing the Bernoulli q-density means of the linear component coefficients spike-and-slab auxiliary indicator variables.\cr
   \code{sigmaBeta} = a vector with 2 entries, consisting of the Inverse Gamma q-density shape and rate parameters of the variance of the linear component coefficients.\cr
   \code{rhoBeta} = a vector with 2 entries, consisting of the Beta q-density shape parameters of the Bernoulli probability parameter of the linear component coefficients spike-and-slab auxiliary indicator variables.\cr
   \code{uTilde} = a two-component list containing the Multivariate Normal q-density parameters of the spline basis function coefficients without multiplication by the means of the \code{gammaU} q-densities. The list components are: \code{mu.q.uTilde}, the mean vectors for each predictor that is treated as potentially having a non-linear effect; \code{sigsq.q.uTilde}, the diagonal entries of the covariance matrices of each predictor that is treated as potentially having a non-linear effect.\cr
   \code{gammaU} = a list containing the q-density means of the spline basis coefficients spike-and-slab auxiliary indicator variables. Each list component corresponds to a predictor that is treated as potentially having a non-linear effect.\cr
   \code{rhoU} = a two-component list with components \code{A.q.rho.u} and \code{B.q.rho.u}. The \code{A.q.rho.u} list component is a vector of Beta q-density first (one plus the power of rho) shape parameters corresponding to the spline basis coefficients spike-and-slab auxiliary indicator variables for each predictor that is treated as potentially having a non-linear effect. The \code{B.q.rho.u} list component is a vector of Beta q-density second (one plus the power of 1-rho) shape parameters corresponding to the spline basis coefficients spike-and-slab auxiliary indicator variables for each predictor that is treated as potentially having a non-linear effect.\cr
   \code{sigmaEps} = a vector with 2 entries, consisting of the Inverse Gamma q-density shape and rate parameters of the error variance.
} 
\item{effectTypeHat}{an array of character strings, with entry either "zero", "linear" or "nonlinear", signifying the estimated effect type for each candidate predictor.} 
\item{meanXlinear}{an array containing the sample means of each column of \code{Xlinear}.} 
\item{sdXlinear}{an array containing the sample standard deviations of each column of \code{Xlinear}.} 
\item{meanXgeneral}{an array containing the sample means of each column of \code{Xgeneral}.} 
\item{sdXgeneral}{an array containing the sample standard deviations of each column of \code{Xgeneral}.} 
}
\author{Virginia X. He \email{virginia.x.he@student.uts.edu.au} and Matt P. Wand \email{matt.wand@uts.edu.au}}
\references{
Chouldechova, A. and Hastie, T. (2015). Generalized additive model selection. <arXiv:1506.03850v2>.

He, V.X. and Wand, M.P. (2021). Bayesian generalized additive model selection including a fast variational option. <arXiv:2021.PLACE-HOLDER>. 
}
\examples{library(gamselBayes) 

# Generate some simple regression-type data:

set.seed(1) ; n <- 1000 ; x1 <- rbinom(n,1,0.5) ; 
x2 <- runif(n) ; x3 <- runif(n) ; x4 <- runif(n)
y <- x1 + sin(2*pi*x2) - x3 + rnorm(n)
Xlinear <- data.frame(x1) ; Xgeneral <- data.frame(x2,x3,x4)

# Obtain a gamselBayes() fit for the data, using Markov chain Monte Carlo:

fitMCMC <- gamselBayes(y,Xlinear,Xgeneral)
summary(fitMCMC) ; plot(fitMCMC) ; checkChains(fitMCMC)

# Obtain a gamselBayes() fit for the data, using mean field variational Bayes:

fitMFVB <- gamselBayes(y,Xlinear,Xgeneral,method = "MFVB")
summary(fitMFVB) ; plot(fitMFVB)

\donttest{if (require("Ecdat"))
{
   # Obtain a gamselBayes() fit for data on schools in California, U.S.A.:

   Caschool$log.avginc <- log(Caschool$avginc)
   mathScore <- Caschool$mathscr
   Xgeneral <- Caschool[,c("mealpct","elpct","calwpct","compstu","log.avginc")]

   # Obtain a gamselBayes() fit for the data, using Markov chain Monte Carlo:

   fitMCMC <- gamselBayes(y = mathScore,Xgeneral = Xgeneral)
   summary(fitMCMC) ; plot(fitMCMC) ; checkChains(fitMCMC)

   # Obtain a gamselBayes() fit for the data, using mean field variational Bayes:

   fitMFVB <- gamselBayes(y = mathScore,Xgeneral = Xgeneral,method = "MFVB")
   summary(fitMFVB) ; plot(fitMFVB)
}}
}

