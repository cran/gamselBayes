/********** C++ function: omitVecEnt **********/

/* Omits a specified entry from a column vector.  */

/* Last changed: 07 JAN 2021 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

arma::vec omitVecEnt(arma::vec a,int j)
{
  int dmn = a.n_elem;
  arma::vec ans(dmn-1);

  if (j==0) 
  {
     ans = a.rows(1,dmn-1);
  }
  if ((j>0)&(j<(dmn-1)))
  {
     ans = arma::join_cols(a.rows(0,j-1),a.rows(j+1,dmn-1));
  } 
  if (j==(dmn-1)) 
  {
     ans = a.rows(0,dmn-2);
  }

  return ans;
}

/************ End of omitVecEnt ************/

