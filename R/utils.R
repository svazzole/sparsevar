#' @title L2 matrix norm
#'
#' @description Compute the L2 matrix norm of M
#' @usage l2norm(m)
#' @param m the matrix (real or complex valued)
#'
#' @export
l2norm <- function(m) {
  sqrt(spectral_radius(t(m) %*% m))
}

#' @title L1 matrix norm
#'
#' @description Compute the L1 matrix norm of m
#' @usage l1norm(m)
#' @param m the matrix (real or complex valued)
#'
#' @export
l1norm <- function(m) {
  max(colSums(Mod(m)))
}

#' @title L-infinity matrix norm
#'
#' @description Compute the L-infinity matrix norm of m
#' @usage l_infty_norm(m)
#' @param m the matrix (real or complex valued)
#'
#' @export
l_infty_norm <- function(m) {
  max(rowSums(Mod(m)))
}

#' @title Max-norm of a matrix
#'
#' @description Compute the max-norm of m
#' @usage max_norm(m)
#' @param m the matrix (real or complex valued)
#'
#' @export
max_norm <- function(m) {
  max(abs(m))
}

#' @title Frobenius norm of a matrix
#'
#' @description Compute the Frobenius norm of m
#' @usage frob_norm(M)
#' @param m the matrix (real or complex valued)
#'
#' @export
frob_norm <- function(m) {
  a <- (t(m) %*% m)
  a <- a * diag(nrow(a))
  sqrt(sum(a))
}

#' @title Spectral radius
#'
#' @description Compute the spectral radius of m
#' @usage spectral_radius(m)
#' @param m the matrix (real or complex valued)
#'
#' @export
spectral_radius <- function(m) {
  e <- eigen(m)
  max_eig <- max(Mod(e$values))
  max_eig
}

#' @title Spectral norm
#'
#' @description Compute the spectral norm of m
#' @usage spectral_norm(m)
#' @param m the matrix (real or complex valued)
#'
#' @export
spectral_norm <- function(m) {
  sqrt(spectral_radius(t(m) %*% m))
}

#' @title Accuracy metric
#'
#' @description Compute the accuracy of a fit
#' @param reference_m the matrix to use as reference
#' @param a the matrix obtained from a fit
#'
#' @usage accuracy(reference_m, a)
#'
#' @export
accuracy <- function(reference_m, a) {
  n <- ncol(a)
  l <- a
  l[l != 0] <- 1
  l[l == 0] <- 0

  gen_l <- reference_m
  gen_l[gen_l != 0] <- 1
  gen_l[gen_l == 0] <- 0

  acc <- 1 - sum(abs(l - gen_l)) / n^2
  acc
}

#' @title True Negative Rate
#'
#' @description Computes the True Negative Rate (Specificity) between a
#'   reference matrix and an estimated matrix. TNR = TN / (TN + FP), where
#'   negatives are zero entries in the reference matrix.
#'
#' @param reference_m the reference (ground truth) matrix
#' @param a the estimated matrix to compare against the reference
#'
#' @return The true negative rate (between 0 and 1), or NA if there are no
#'   actual negatives in the reference matrix
#'
#' @export
true_negative_rate <- function(reference_m, a) {
  actual_negatives <- reference_m == 0
  num_negatives <- sum(actual_negatives)

  if (num_negatives == 0) {
    return(NA_real_)
  }

  true_negatives <- sum(actual_negatives & a == 0)
  true_negatives / num_negatives
}

#' @title True Positive Rate
#'
#' @description Computes the True Positive Rate (Sensitivity/Recall) between a
#'   reference matrix and an estimated matrix. TPR = TP / (TP + FN), where
#'   positives are non-zero entries in the reference matrix.
#'
#' @param reference_m the reference (ground truth) matrix
#' @param a the estimated matrix to compare against the reference
#'
#' @return The true positive rate (between 0 and 1), or NA if there are no
#'   actual positives in the reference matrix
#'
#' @export
true_positive_rate <- function(reference_m, a) {
  actual_positives <- reference_m != 0
  num_positives <- sum(actual_positives)

  if (num_positives == 0) {
    return(NA_real_)
  }

  true_positives <- sum(actual_positives & a != 0)
  true_positives / num_positives
}

#' @title Check is var
#'
#' @description Check if the input is a var object
#' @param v the object to test
#'
#' @usage check_is_var(v)
#'
#' @export
check_is_var <- function(v) {
  if (!is.null(attr(v, "class"))) {
    ifelse(attr(v, "class") == "var" || attr(v, "class") == "varx",
           return(TRUE), return(FALSE))
  }
  FALSE
}
