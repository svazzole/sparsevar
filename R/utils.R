#' @title L2 matrix norm
#'
#' @description Compute the L2 matrix norm of M
#' @usage l2norm(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
l2norm <- function(M) {
  s <- sqrt(spectralRadius(t(M) %*% M))
  return(s)
}

#' @title L1 matrix norm
#'
#' @description Compute the L1 matrix norm of M
#' @usage l1norm(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
l1norm <- function(M) {
  c <- max(colSums(Mod(M)))
  return(c)
}

#' @title L-infinity matrix norm
#'
#' @description Compute the L-infinity matrix norm of M
#' @usage lInftyNorm(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
lInftyNorm <- function(M) {
  c <- max(rowSums(Mod(M)))
  return(c)
}

#' @title Max-norm of a matrix
#'
#' @description Compute the max-norm of M
#' @usage maxNorm(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
maxNorm <- function(M) {
  return(max(abs(M)))
}

#' @title Froebenius norm of a matrix
#'
#' @description Compute the Froebenius norm of M
#' @usage frobNorm(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
frobNorm <- function(M) {
  A <- (t(M) %*% M)
  A <- A * diag(nrow(A))
  return(sqrt(sum(A)))
}

#' @title Spectral radius
#'
#' @description Compute the spectral radius of M
#' @usage spectralRadius(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
spectralRadius <- function(M) {
  e <- eigen(M)
  maxEig <- max(Mod(e$values))
  return(maxEig)
}

#' @title Spectral norm
#'
#' @description Compute the spectral norm of M
#' @usage spectralNorm(M)
#' @param M the matrix (real or complex valued)
#'
#' @export
spectralNorm <- function(M) {
  return(sqrt(spectralRadius(t(M) %*% M)))
}

#' @title Accuracy metric
#'
#' @description Compute the accuracy of a fit
#' @param referenceM the matrix to use as reference
#' @param A the matrix obtained from a fit
#'
#' @usage accuracy(referenceM, A)
#'
#' @export
accuracy <- function(referenceM, A) {
  N <- ncol(A)
  L <- A
  L[L != 0] <- 1
  L[L == 0] <- 0

  genL <- referenceM
  genL[genL != 0] <- 1
  genL[genL == 0] <- 0

  acc <- 1 - sum(abs(L - genL)) / N^2 # accuracy    -(1 - sum(genL)/N^2)
  return(acc)
}

#' @title Check is var
#'
#' @description Check if the input is a var object
#' @param v the object to test
#'
#' @usage checkIsVar(v)
#'
#' @export
checkIsVar <- function(v) {
  if (!is.null(attr(v, "class"))) {
    ifelse(attr(v, "class") == "var" | attr(v, "class") == "varx", return(TRUE), return(FALSE))
  } else {
    return(FALSE)
  }
}
