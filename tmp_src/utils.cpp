#include <Rcpp.h>
using namespace Rcpp;

// Cross product of y with jth column of X
double crossprod(NumericVector X, NumericVector y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) {
    val += X(nn+i)*y(i); 
  }
  return(val);
}

// Weighted cross product of y with jth column of x
double wcrossprod(NumericVector X, NumericVector y, NumericVector w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) {
    val += X(nn+i)*y(i)*w(i); 
  }
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(NumericVector X, NumericVector w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) {
    val += w(i) * pow(X(nn+i), 2); 
  }
  return(val);
}

// Sum of squares of jth column of X
double sqsum(NumericVector X, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) {
    val += pow(X(nn+i), 2);
  }
  return(val);
}

double sum(NumericVector x, int n) {
  double val=0;
  for (int i=0;i<n;i++) {
    val += x(i);  
  }
  return(val);
}

int checkConvergence(NumericVector beta, NumericVector beta_old, double eps, int l, int J) {
  int converged = 1;
  for (int j=0; j<J; j++) {
    if (fabs((beta(l*J+j)-beta_old(j)/beta_old(j))) > eps) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

double MCP(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(v*(1+l2-1/gamma)));
  else return(z/(v*(1+l2)));
}



double LASSO(double z, double l1, double l2, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else return(s*(fabs(z)-l1)/(v*(1+l2)));
}
