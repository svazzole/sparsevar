#' @title Matrix norms
#' 
#' @description Some matrix norms
#' @usage 
#' \code{l2norm(M)}
#' \code{l1norm(M)}
#' \code{lInftyNorm(M)}
#' \code{maxNorm(M)}
#' \code{frobNorm(M)}
#' \code{spectralRadius(M)}
#' \code{spectralNorm(M)}
#' @param M the matrix (real or complex valued)
#' 
#' @author Simone Vazzoler

#' @export
l2norm <- function(M) {
  s <- sqrt(spectralRadius(t(M) %*% M))
  return(s)
}

#' @export
l1norm <- function(M) {
  c <- max(colSums(Mod(M)))
  return(c)
}

#' @export
lInftyNorm <- function(M) {
  c <- max(rowSums(Mod(M)))
  return(c)
}

#' @export
maxNorm <- function(M) {
  return(max(abs(M)))
}

#' @export
frobNorm <- function(M) {
  A <- (t(M) %*% M)
  A <-  A * diag(nrow(A))
  return(sqrt(sum(A)))
}

#' @export
spectralRadius <- function(M) {
  e <- eigen(M)
  maxEig <- max(Mod(e$values))
  return(maxEig)
}

#' @export
spectralNorm <- function(M) {
  return(sqrt(spectralRadius(t(M) %*% M)))
}
