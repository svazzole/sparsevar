// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>

// [[Rcpp::export]]
double gLoss(arma::colvec r) {
  // Gaussian loss
  double l = arma::sum(r % r);
  return l;
}

// [[Rcpp::export]]
arma::colvec crossprod(arma::sp_mat X, arma::colvec y) {
  // Cross product of y with jth column of X
  int n = y.n_rows;
  arma::mat ret = (arma::conv_to<arma::rowvec>::from(y)* X)/n;
  return arma::conv_to< arma::colvec >::from(ret);
}

// Sum of squares of jth column of X
// double sqsum(arma::mat X, int j) {
//     return arma::accu(arma::pow(X.col(j),2));
// }

double SCAD(double z, double l1, double l2, double gamma, double v) {

  double s = 0;

  if (z > 0){
    s = 1;
  } else {
    s = -1;
  }

  if (fabs(z) <= l1) {
    return 0;
  } else if (fabs(z) <= (l1*(1+l2)+l1)){
    return s*(fabs(z)-l1)/(v*(1+l2));
  } else if (fabs(z) <= gamma*l1*(1+l2)) {
    return s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2));
  } else {
    return z/(v*(1+l2));
  }

}

int checkConvergence(arma::rowvec b, arma::rowvec beta_old, double eps, int p) {

  for (int j=0; j<p; j++) {
    if (fabs((b(j)-beta_old(j))/beta_old(j)) > eps) {
      return 0;
    }
  }
  return 1;
}

//[[Rcpp::export]]
Rcpp::List cdfit_gaussianTEST(arma::sp_mat X, arma::colvec Y, arma::colvec lambda, double eps, int max_iter, double gamma, arma::colvec multiplier, double alpha, int dfmax, int user) {

  int n = Y.n_rows;
  int p = X.n_cols;
  int L = lambda.n_elem;

  //arma::sp_mat Xt = arma::trans(X);
  arma::mat beta(L, p, arma::fill::zeros);
  arma::colvec loss(L, arma::fill::zeros);
  arma::uvec iter(L, arma::fill::zeros);
  arma::rowvec a(p, arma::fill::zeros);
  arma::colvec r = Y;
  arma::colvec z = crossprod(X,r);

  arma::uvec e1(p, arma::fill::zeros);
  arma::uvec e2(p, arma::fill::zeros);
  double cutoff, l1, l2;
  int converged, lstart;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    loss(0) = gLoss(r);
    lstart = 1;
  }

  // Path
  for (int l=lstart;l<L;l++) {
    Rcpp::checkUserInterrupt();
    if (l != 0) {
      // Assign a
      a = beta.row(l-1); // Check l or l-1?

      // Check dfmax
      int nv = 0;
      arma::uvec f = arma::find(arma::abs(a) > 0);
      nv = f.n_elem;
      if (nv > dfmax) {
        for (int ll=l; ll<L; ll++) {
          iter(ll) = arma::datum::nan;
        }
        Rcpp::List lst = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                            Rcpp::Named("loss") = loss,
                                            Rcpp::Named("iter") = iter);
        return lst;
      }
      // Determine eligible set
      // if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lam[l-1];
      // if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]);
      // if (1 == 1) {

      cutoff = lambda(l) + gamma/(gamma-2)*(lambda(l) - lambda(l-1));
      //}

      e2 = (arma::abs(z) > cutoff*alpha*multiplier);

    } else {
        // Determine eligible set
        double lmax = arma::max(arma::abs(z));
        cutoff = lambda(l) + gamma/(gamma-2)*(lambda(l) - lmax);
        e2 = (arma::abs(z) > cutoff*alpha*multiplier);
    }

    while (iter(l) < max_iter) {
      while (iter(l) < max_iter) {
        while (iter(l) < max_iter) {
          // Solve over the active set
          iter(l)++;

          z = crossprod(X,r) + arma::trans(a);
          
          for (int j=0; j<p; j++) {
            if (e1(j)) {
              // Update beta_j
              l1 = lambda(l) * multiplier(j) * alpha;
              l2 = lambda(l) * multiplier(j) * (1-alpha);
              beta(l,j) = SCAD(z(j), l1, l2, gamma, 1);
              // Update r
              double shift = beta(l,j) - a(j);
              if (shift !=0) {
                for (int i=0;i<n;i++) {
                  r(i) -= shift*X(i,j);
                }
              }
            }
          }

          // Check for convergence
          converged = checkConvergence(beta.row(l), a, eps, p);
          a = beta.row(l);
          if (converged) {
            break;
          }
        }
        
        // Scan for violations in strong set
        int violations = 0;
        z = crossprod(X,r);
        arma::uvec t = (e1 == 0) && (e2 == 1);
        arma::uvec ix = arma::find(t);

        for (int j=0; j<ix.n_elem; j++) {
          // Update beta_j
          l1 = lambda(l) * multiplier(ix(j)) * alpha;
          l2 = lambda(l) * multiplier(ix(j)) * (1-alpha);
          beta(l,ix(j)) = SCAD(z(ix(j)), l1, l2, gamma, 1);

          // If something enters the eligible set, update eligible set & residuals
          if (beta(l,ix(j)) !=0) {
            e1(ix(j)) = e2(ix(j)) = 1;
            for (int i=0; i<n; i++) {
              r(i) -= beta(l,ix(j))*X(i,ix(j));
            }
            a(ix(j)) = beta(l,ix(j));
            violations++;
          }
        }
        if (violations==0) {
          break;
        }

      }
      // Scan for violations in rest
      int violations = 0;
      z = crossprod(X,r);
      arma::uvec ix = arma::find(e2 == 0);

      for (int j=0; j<ix.n_elem; j++) {
        // Update beta_j
        l1 = lambda(l) * multiplier(ix(j)) * alpha;
        l2 = lambda(l) * multiplier(ix(j)) * (1-alpha);
        beta(l,ix(j)) = SCAD(z(ix(j)), l1, l2, gamma, 1);
        // If something enters the eligible set, update eligible set & residuals
        if (beta(l,ix(j))!=0) {
          e1(ix(j)) = e2(ix(j)) = 1;
          for (int i=0; i<n; i++) {
            r(i) -= beta(l,ix(j))*X(i,ix(j));
          }
          a(ix(j)) = beta(l,ix(j));
          violations++;
        }

      }

      if (violations==0) {
        loss(l) = gLoss(r);
        break;
      }

    }
  }
  Rcpp::List lst = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                      Rcpp::Named("loss") = loss,
                                      Rcpp::Named("iter") = iter);
  return lst;
}

