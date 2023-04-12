/********** C++ function: drawInvGaussVec **********/

/* Draws a vector of Inverse Gaussian variates corresponding
   to a given vector of mean values (the precision parameter
   is unity). */

/* Last changed: 22 FEB 2021 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

arma::vec drawInvGaussVec(arma::vec muVec)
{
   /* Declare all non-input variables: */

   int nSamp = muVec.n_elem;
   arma::vec z;
   arma::vec x;
   arma::vec u;
   arma::vec muSqd;
   arma::vec omega1;
   arma::vec omega2;
   arma::vec omega3;

   /* Generate vector of Inverse-Gaussian draws: */
   
   z = rnorm(nSamp);
   u = runif(nSamp);
   muSqd = muVec%muVec;
   omega1 = z%z;
   omega2 = muVec + 0.5*muSqd%omega1 
            - 0.5*muVec%sqrt(4.0*muVec%omega1 + muSqd%omega1%omega1);
    
   omega3 = arma::zeros(nSamp);
   for (int i = 0; i < nSamp; i++)
   {
      if (u(i)<=muVec(i)/(muVec(i) + omega2(i))) 
         omega3(i) = 1.0;
   }
   
   x = omega2%omega3 + (muSqd/omega2)%(1.0-omega3);

   return x;
}

/************ End of drawInvGaussVec ************/

