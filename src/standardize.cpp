// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// using namespace Rcpp;
// using namespace arma;

// [[Rcpp::export]]
Rcpp::List standardize2(arma::mat X) {

  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat XX(n, p, arma::fill::zeros);
  arma::rowvec c(p, arma::fill::zeros);
  arma::rowvec s(p, arma::fill::zeros);

  for (int j=0; j<p; j++) {
    // Center
    for (int i=0; i<n; i++) {
      c(j) += X(i,j);
    }
    c(j) = c(j) / n;
    for (int i=0; i<n; i++) {
      XX(i,j) = X(i,j) - c(j);
    }
    // Scale
    for (int i=0; i<n; i++) {
      s(j) += pow(XX(i,j), 2);
    }
    s(j) = sqrt(s(j)/n);
    for (int i=0; i<n; i++) {
      XX(i,j) = XX(i,j)/s(j);
    }
  }

  Rcpp::List l = Rcpp::List::create(Rcpp::Named("Data")   = XX,
                                    Rcpp::Named("Means")  = c,
                                    Rcpp::Named("Sd")  = s);
  return l;
}
