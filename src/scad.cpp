// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Gaussian loss
// [[Rcpp::export]]
double gLoss(arma::colvec r) {
  int n = r.n_rows;
  double l = 0;
  for (int i=0;i<n;i++) {
    l = l + pow(r(i),2); 
  }
  return(l);
}

// Cross product of y with jth column of X
// [[Rcpp::export]]
arma::colvec crossprod(arma::sp_mat X, arma::colvec y) {
  //int nn = n*j;
  int n = y.n_rows;
  // double val=0;
  // for (int i=0;i<n;i++) {
  //   val += X(i,j)*y(i); 
  // }
  // return(val);
  arma::colvec ret = (arma::trans(X) * y)/n;
  return ret;
}

double SCAD(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(v*(1+l2)));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)));
  else return(z/(v*(1+l2)));
}


// [[Rcpp::export]]
arma::sp_mat cdfit_gaussian(arma::sp_mat X, arma::colvec Y, arma::colvec lambda, double eps, int max_iter, double gamma, arma::colvec multiplier, double alpha, int dfmax, int user) {
//Rcpp::NumericMatrix cdfit_gaussian() {
  
  int n = Y.n_rows;
  int p = X.n_cols;
  int L = lambda.n_rows;
  
  arma::sp_mat beta(L,p);
  arma::colvec loss(L, arma::fill::zeros);
  arma::colvec iter(L, arma::fill::zeros);
  arma::colvec a(p, arma::fill::zeros);
  arma::colvec r = Y;
  arma::colvec z(p, arma::fill::zeros);
  arma::colvec e1(p, arma::fill::zeros);
  arma::colvec e2(p, arma::fill::zeros);
  double cutoff, l1, l2;
  int converged, lstart;
  
  z = crossprod(arma::trans(X),r);
  // for (int j=0; j<p; j++){
  //   z(j) = crossprod(X, r, n, j)/n; 
  // }
  // 
  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    loss(0) = gLoss(r);
    lstart = 1;
  }
  
  // Path
  for (int l=lstart;l<L;l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a
      for (int j=0;j<p;j++) {
        a(j) = beta(l,j); 
      }
      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
        if (a(j) != 0) nv++;
      }
      if (nv > dfmax) {
        for (int ll=l; ll<L; ll++) {
          iter(ll) = 0; 
        }
        return beta;
        // res = cleanupG(a, r, e1, e2, z, beta, loss, iter);
        // return(res);
      }
      // Determine eligible set
      // if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lam[l-1];
      // if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]);
      if (1 == 1) {
        cutoff = lambda(l) + gamma/(gamma-2)*(lambda(l) - lambda(l-1));
      } 
      for (int j=0; j<p; j++) {
        if (fabs(z(j)) > (cutoff * alpha * multiplier(j))) {
          e2(j) = 1;
        }
      }
    } else {
        // Determine eligible set
        double lmax = 0;
        for (int j=0; j<p; j++) {
          if (fabs(z[j]) > lmax) {
            lmax = fabs(z[j]);
          }
        }
        cutoff = lambda(l) + gamma/(gamma-2)*(lambda(l) - lmax);
        for (int j=0; j<p; j++) {
          if (fabs(z(j)) > (cutoff * alpha * multiplier(j))) { 
            e2(j) = 1;
          }
        }
    }
    
    while (iter(l) < max_iter) {
      while (iter(l) < max_iter) {
        while (iter(l) < max_iter) {
          // Solve over the active set
          iter(l)++;
          z = crossprod(X,r) + a;
          
          for (int j=0; j<p; j++) {
            if (e1(j)) {
              //z(j) = crossprod(X,r) + a(j);
              
              // Update beta_j
              l1 = lambda(l) * multiplier(j) * alpha;
              l2 = lambda(l) * multiplier(j) * (1-alpha);
              beta(l,j) = SCAD(z(j), l1, l2, gamma, 1);
              
              // Update r
              double shift = beta(l,j) - a(j);
              if (shift !=0) for (int i=0;i<n;i++) r(i) -= shift*X[j*n+i];
            }
          }
        }
      }
    }
  }
  return(beta);
}

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
