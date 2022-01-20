/********** C++ function: zetad **********/

/* Computes the first derivative of the zeta function.  */

/* Last changed: 06 DEC 2021 */

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

double zetad(double x)
{
   /* Declare all non-input variables: */

   double piValue;
   double rt2;
   double toler;
   double ansprv;
   double anscur;
   double ans;
   double tiny;
   double Cprv;
   double Ccur;
   double Dprv;
   double Dcur;
   double Delta;
   int j;
   double aj;

   if (x>(-3.0))
   {  
      piValue = 4.0*atan(1.0);
      rt2 = sqrt(2.0);      
      ans = 2.0*exp(-0.5*x*x)/(sqrt(2.0*piValue)*erfc(-x/rt2));
   }
   else
   {
      toler = 1.0e-10;
      tiny = 1.0e-30;
      ansprv = tiny;
      Cprv = tiny;
      Dprv = 0.0;
      Delta = 2.0 + toler;
      j = 0;
      while (abs(Delta-1.0)>toler)
      {
         j = j + 1;
         aj = j - 1;
         if (j==1)
         { 
            aj = 1.0;
         }
         Dcur = aj*Dprv - x;
         if (abs(Dcur)<tiny) 
         {
            Dcur = tiny;
         }
         Dcur = 1/Dcur;
         Ccur = (aj/Cprv) - x;
         if (abs(Ccur)<tiny)
         {
            Ccur = tiny;
         }
         Delta = Ccur*Dcur;
         anscur = ansprv*Delta;
         ansprv = anscur;
         Cprv = Ccur;
         Dprv = Dcur;
 
         ans = 1/anscur;
      }
   }

   return ans;
}

/************ End of zetad ************/



