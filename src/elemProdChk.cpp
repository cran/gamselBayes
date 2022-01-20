/********** C++ function: elemProdChk **********/

/* Obtains the element-wise product of two vectors */

/* Last changed: 20 NOV 2020 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

arma::vec elemProdChk(arma::vec a,arma::vec b)
{
   /* Declare all non-input variables: */

   int len = a.n_elem;
   arma::vec ans(len);
   arma::mat ooooo(len,len);
       
   ooooo = 8.953*arma::ones(len,len);

   Rcout << ooooo << "\n";

   ans = a%b;

   return ans;
}

/************ End of elemProdChk ************/

