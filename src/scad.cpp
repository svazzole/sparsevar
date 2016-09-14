// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Gaussian loss
// [[Rcpp::export]]
double gLoss(NumericVector r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) {
    l = l + pow(r(i),2); 
  }
  return(l);
}

// [[Rcpp::export]]
arma::sp_mat tp(const arma::sp_mat X, const arma::sp_mat Y) {
  //int n = X.n_rows, k = X.n_cols;
  arma::sp_mat out = X * Y;
  return out;
}

// [[Rcpp::export]]
arma::sp_mat cdfit_gaussian(const arma::sp_mat X, const arma::sp_mat Y, NumericVector lambda, double eps, int max_iter, double gamma, double multiplier, double alpha, int dfmax, SEXP user) {
  int n = Y.n_rows;
  int p = X.n_cols;
  int L = lambda.length();
  
  
  arma::sp_mat out = X * Y;
  return out;
}

// // Coordinate descent for gaussian models
// extern "C" SEXP cdfit_gaussian(const arma::sp_mat X, const arma::sp_mat y, SEXP penalty, double lambda, double eps, int max_iter, double gamma, double multiplier, double alpha, int dfmax, SEXP user_) {
//   
// //SEXP cdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_) {
//   
//   Rcpp::NumericVector yy(y);
//   int n = y.n_rows;
//   
//   // Declarations
//   //int n = yy.nrow(); //length(y_);
//   int p = length(X_)/n;
//   int L = length(lambda);
//   SEXP res, beta, loss, iter;
//   PROTECT(beta = allocVector(REALSXP, L*p));
//   double *b = REAL(beta);
//   for (int j=0; j<(L*p); j++) b[j] = 0;
//   PROTECT(loss = allocVector(REALSXP, L));
//   PROTECT(iter = allocVector(INTSXP, L));
//   for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
//   double *a = Calloc(p, double); // Beta from previous iteration
//   for (int j=0; j<p; j++) a[j]=0;
//   double *X = REAL(X_);
//   double *y = REAL(y_);
//   const char *penalty = CHAR(STRING_ELT(penalty_, 0));
//   double *lam = REAL(lambda);
//   double eps = REAL(eps_)[0];
//   int max_iter = INTEGER(max_iter_)[0];
//   double gamma = REAL(gamma_)[0];
//   double *m = REAL(multiplier);
//   double alpha = REAL(alpha_)[0];
//   int dfmax = INTEGER(dfmax_)[0];
//   int user = INTEGER(user_)[0];
//   double *r = Calloc(n, double);
//   for (int i=0; i<n; i++) r[i] = y[i];
//   double *z = Calloc(p, double);
//   for (int j=0; j<p; j++) z[j] = crossprod(X, r, n, j)/n;
//   int *e1 = Calloc(p, int);
//   for (int j=0; j<p; j++) e1[j] = 0;
//   int *e2 = Calloc(p, int);
//   for (int j=0; j<p; j++) e2[j] = 0;
//   double cutoff, l1, l2;
//   int converged, lstart;
//   
//   // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
//   if (user) {
//     lstart = 0;
//   } else {
//     REAL(loss)[0] = gLoss(r,n);
//     lstart = 1;
//   }
//   
//   // Path
//   for (int l=lstart;l<L;l++) {
//     R_CheckUserInterrupt();
//     if (l != 0) {
//       // Assign a
//       for (int j=0;j<p;j++) a[j] = b[(l-1)*p+j];
//       
//       // Check dfmax
//       int nv = 0;
//       for (int j=0; j<p; j++) {
//         if (a[j] != 0) nv++;
//       }
//       if (nv > dfmax) {
//         for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
//         res = cleanupG(a, r, e1, e2, z, beta, loss, iter);
//         return(res);
//       }
//       
//       // Determine eligible set
//       if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lam[l-1];
//       if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]);
//       if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1]);
//       for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
//     } else {
//       // Determine eligible set
//       double lmax = 0;
//       for (int j=0; j<p; j++) if (fabs(z[j]) > lmax) lmax = fabs(z[j]);
//       if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lmax;
//       if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax);
//       if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax);
//       for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
//     }
//     
//     while (INTEGER(iter)[l] < max_iter) {
//       while (INTEGER(iter)[l] < max_iter) {
//         while (INTEGER(iter)[l] < max_iter) {
//           // Solve over the active set
//           INTEGER(iter)[l]++;
//           for (int j=0; j<p; j++) {
//             if (e1[j]) {
//               z[j] = crossprod(X, r, n, j)/n + a[j];
//               
//               // Update beta_j
//               l1 = lam[l] * m[j] * alpha;
//               l2 = lam[l] * m[j] * (1-alpha);
//               if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
//               if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
//               if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);
//               
//               // Update r
//               double shift = b[l*p+j] - a[j];
//               if (shift !=0) for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
//             }
//           }
//           
//           // Check for convergence
//           converged = checkConvergence(b, a, eps, l, p);
//           for (int j=0; j<p; j++) a[j] = b[l*p+j];
//           if (converged) break;
//         }
//         
//         // Scan for violations in strong set
//         int violations = 0;
//         for (int j=0; j<p; j++) {
//           if (e1[j]==0 & e2[j]==1) {
//             
//             z[j] = crossprod(X, r, n, j)/n;
//             
//             // Update beta_j
//             l1 = lam[l] * m[j] * alpha;
//             l2 = lam[l] * m[j] * (1-alpha);
//             if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
//             if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
//             if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);
//             
//             // If something enters the eligible set, update eligible set & residuals
//             if (b[l*p+j] !=0) {
//               e1[j] = e2[j] = 1;
//               for (int i=0; i<n; i++) r[i] -= b[l*p+j]*X[j*n+i];
//               a[j] = b[l*p+j];
//               violations++;
//             }
//           }
//         }
//         if (violations==0) break;
//       }
//       
//       // Scan for violations in rest
//       int violations = 0;
//       for (int j=0; j<p; j++) {
//         if (e2[j]==0) {
//           
//           z[j] = crossprod(X, r, n, j)/n;
//           
//           // Update beta_j
//           l1 = lam[l] * m[j] * alpha;
//           l2 = lam[l] * m[j] * (1-alpha);
//           if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
//           if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
//           if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);
//           
//           // If something enters the eligible set, update eligible set & residuals
//           if (b[l*p+j] !=0) {
//             e1[j] = e2[j] = 1;
//             for (int i=0; i<n; i++) r[i] -= b[l*p+j]*X[j*n+i];
//             a[j] = b[l*p+j];
//             violations++;
//           }
//         }
//       }
//       
//       if (violations==0) {
//         REAL(loss)[l] = gLoss(r, n);
//         break;
//       }
//     }
//   }
//   res = cleanupG(a, r, e1, e2, z, beta, loss, iter);
//   return(res);
// }
