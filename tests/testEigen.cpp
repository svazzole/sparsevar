// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


using namespace Rcpp;
using namespace arma;

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef MSpMat::InnerIterator InIterMat;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

// [[Rcpp::export]]
void mat_loop2 (MSpMat X, int j){
  Rcout << "Standard looping over a sparse matrix" << std::endl;
  for (InIterMat i_(X, j); i_; ++i_){
    Rcout << " i,j=" << i_.index() << "," << j << " value=" << i_.value() << std::endl;
  }
}

// [[Rcpp::export]]
double matTimesVec (MSpMat X, colvec x, int j) {
  
  double s = 0;
  
  for (int j = 0; j<6; j++){
    for (InIterMat i_(X, j); i_; ++i_){
      s += i_.value() * x(i_.index());
    }
  }
  return s;  
}

/*** R
v <- rep(1,6)
library(Matrix)
M1<- new("dgCMatrix"
         , i = c(1L, 2L, 3L, 0L, 2L, 3L, 0L, 1L, 3L, 0L, 
                 1L, 2L, 4L, 5L, 3L, 5L, 3L, 4L)
         , p = c(0L, 3L, 6L, 9L, 14L, 16L, 18L)
         , Dim = c(6L, 6L)
         , Dimnames = list(c("a", "b", "c", "d", "e", "f"), 
                           c("a", "b", "c", "d", "e", "f"))
         , x = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3)
         , factors = list()) 
M1
mat_loop2( M1, 2 )
matTimesVec (M1, v, 0)
*/
