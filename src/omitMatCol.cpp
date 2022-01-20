/********** C++ function: omitMatCol **********/

/* Omits a specified column from a matrix.  */

/* Last changed: 08 JAN 2021 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat omitMatCol(arma::mat A, int j)
{
  int numCols = A.n_cols;
  arma::mat ans(A.n_rows,numCols-1);
  
  if (j==0)
  {
     ans = A.cols(1,numCols-1);
  }
  if ((j>0)&(j<(numCols-1)))
  {
     ans = arma::join_horiz(A.cols(0,j-1),A.cols(j+1,numCols-1));
  }
  if (j==(numCols-1))
  {
     ans = A.cols(0,numCols-2);
  }
  
  return ans;
}

/************ End of omitMatCol ************/



