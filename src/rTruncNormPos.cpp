/********** C++ function: rTruncNormPos **********/
  
/* Obtains a single draw from a Truncated Normal 
   distribution with mean "mu", unit standard 
   deviation and truncation over the positive 
   half-line. */
 
/* Last chaged: 16 JUN 2021 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

double rTruncNormPos(double mu)
{
   bool finished;
   double alphaStar;
   double u1;
   double u2;
   double aux;
   double rhoVal;
   double z;
   double x;
  
   if (mu<0)
   {
      finished = FALSE;
      while (!finished)
      {
         alphaStar = 0.5*(sqrt(mu*mu+4.0)-mu);
         u1 = R::runif(0,1);
         aux = (-mu)-log(u1)/alphaStar;
         rhoVal = exp(-0.5*(aux-alphaStar)*(aux-alphaStar));
         u2 = R::runif(0,1);
         if (u2 <= rhoVal)
         {
            x = aux+mu;
            finished = TRUE;
         }
      }
   }
  
   if (mu>=0)
   {
      finished = FALSE;
      while (!finished)
      {
         z = R::rnorm(0,1);
         if (z > (-mu))
         {
            x = z+mu;
            finished = TRUE;
         }
      }
   }
  
   return (x);
}

/************ End of rTruncNormPos ************/
