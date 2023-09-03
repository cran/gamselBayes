/********** C++ Function approxLogML **********/

/* Computes the mean field variational Bayes approximate logarithm of the
   marginal likelihood of the model corresponding to gamselBayesMFVBinnner. */

/* Last changed: 11 AUG 2023 */

#include <RcppArmadillo.h>

using namespace Rcpp;
#include "logPhi.h"

//[[Rcpp::export]] 

double approxLogML(double muqBeta0, double sigsqqBeta0, double logitBeta, arma::vec muqgammaBeta, 
                   double muqrecipSigsqBeta, arma::vec muqbBeta, arma::vec muqBetaTilde,
                   arma::mat SigmaqBetaTilde, double muqrecipaBeta, double lambdaqSigsqBeta,
                   double sBetaHYP, double lambdaqaBeta, double logitU, arma::vec muqgammaU, 
                   arma::vec muqrecipSigsqU, arma::vec muqbU,arma::mat muqUtilde,
                   arma::mat sigsqqUtilde, arma::vec muqrecipaU, arma::vec lambdaqSigsqU,
                   double sUHYP, arma::vec lambdaqaU, double muqrecipaEps, double muqrecipSigsqEps,
                   double lambdaqaEps, double lambdaqSigsqEps, double sEpsHYP, int n, arma::uvec ncZvec,
                   int familyNum, arma::vec ySign, arma::vec omega20)
{
   /* Declare non-input variables: */

   int ncX = muqbBeta.n_elem;
   int dGeneral = muqbU.n_elem;
   double answer; 
   double sumCurr;
   arma::vec eigvals;

   /* Sequentially obtain the answer: */

   answer = -0.5*(muqBeta0*muqBeta0 + sigsqqBeta0) + 0.5*log(sigsqqBeta0);

   answer = answer + logitBeta*sum(muqgammaBeta);

   for (int j=0; j<ncX; j++)
   {
      if (muqgammaBeta(j)>0.0) answer = answer - muqgammaBeta(j)*log(muqgammaBeta(j));
      if (muqgammaBeta(j)<1.0) answer = answer - (1.0-muqgammaBeta(j))*log(1.0-muqgammaBeta(j)); 
   }

   sumCurr = 0.0;
   for (int j=0; j<ncX; j++)
      sumCurr = sumCurr + muqbBeta(j)*(muqBetaTilde(j)*muqBetaTilde(j) + SigmaqBetaTilde(j,j));
 
   answer = answer - 0.5*muqrecipSigsqBeta*sumCurr;

   eigvals = eig_sym(SigmaqBetaTilde);
   answer = answer + 0.5*sum(log(eigvals));

   for (int j=0; j<ncX; j++)
      answer = answer - 0.5/muqbBeta(j);

   answer = answer - muqrecipaBeta*muqrecipSigsqBeta - 0.5*(ncX + 1.0)*log(lambdaqSigsqBeta);
   answer = answer + muqrecipSigsqBeta*lambdaqSigsqBeta - muqrecipaBeta/(sBetaHYP*sBetaHYP);
   answer = answer + lambdaqaBeta*muqrecipaBeta - log(lambdaqaBeta);

   answer = answer + logitU*sum(muqgammaU);
   
   for (int j=0; j<dGeneral; j++)
   {
     
      if (muqgammaU(j)>0.0) answer = answer - muqgammaU(j)*log(muqgammaU(j));
      if (muqgammaU(j)<1.0) answer = answer - (1.0-muqgammaU(j))*log(1.0-muqgammaU(j)); 

      answer = answer - 0.5*muqrecipSigsqU(j)*muqbU(j)*(sum(muqUtilde.col(j)%muqUtilde.col(j)) 
                                                        + sum(sigsqqUtilde.col(j)));
      for (int k=0; k<ncZvec(j); k++)
         answer = answer + 0.5*log(sigsqqUtilde(k,j));

      answer = answer - 0.5/muqbU(j) - muqrecipaU(j)*muqrecipSigsqU(j);
      answer = answer - 0.5*(ncZvec(j) + 1.0)*log(lambdaqSigsqU(j));
      answer = answer + muqrecipSigsqU(j)*lambdaqSigsqU(j) - muqrecipaU(j)/(sUHYP*sUHYP);
      answer = answer + lambdaqaU(j)*muqrecipaU(j) - log(lambdaqaU(j));      
   }   

   if (familyNum==1)
   {
      answer = answer - 0.5*(n+1)*log(lambdaqSigsqEps) - muqrecipaEps/(sEpsHYP*sEpsHYP);
      answer = answer - log(lambdaqaEps) + lambdaqaEps*muqrecipaEps;
   }

   if (familyNum==2)
   {
      for (int i=0; i<n; i++)
         answer = answer + logPhi(ySign(i)*omega20(i));
   }
   
   /* Return the approximate marginal log-likelihood: */
                                               
   return answer;
} 

/************ End of approxLogML ************/




