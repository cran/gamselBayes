/********** C++ function: logPhi **********/

/* Computes the logarithm of the cumulative distribution function
   of the standard Normal distribution.  */

/* Last changed: 17 OCT 2021 */

#include <Rcpp.h>
#include "zetad.h"

using namespace Rcpp;

// [[Rcpp::export]]

double logPhi(double x)
{
   /* Declare non-input variables: */

   double piValue;
   double logPhiAnswer; 

   if (x>0.0)
   {  
      logPhiAnswer = log(erfc(-x/sqrt(2.0))/2.0);
   }
   else
   {
      piValue = 4.0*atan(1.0);
      logPhiAnswer = -0.5*x*x - log(zetad(x)) - 0.5*log(2.0*piValue);
   }

   return logPhiAnswer;
}

/************ End of logPhi ************/



