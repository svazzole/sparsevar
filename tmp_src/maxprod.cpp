// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// using namespace Rcpp;
// using namespace arma;

double crossprod2(arma::mat X, arma::vec y, int n, int j) {
  double val=0;
  for (int i=0;i<n;i++) val += X(i,j)*y(i);
  return val;
}

// [[Rcpp::export]]
double maxprod(arma::mat X, arma::vec y, arma::uvec v, arma::vec m) {

  // Declarations
  int n = X.n_rows;
  int p = X.n_cols;
  double z = 0;
  double zz;

  for (int j=0; j<p; j++) {
    zz = crossprod2(X, y, n, v(j)-1) / m(v(j)-1);
    if (fabs(zz) > z) {
      z = fabs(zz);
    }
  }

  return z;
}
