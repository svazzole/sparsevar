#' @title Matrix norms
#' 
#' @description Some matrix norms
#' @usage 
#' l2norm(M)
#' l1norm(M)
#' maxNorm(M)
#' frobNorm(M)
#' spectralRadius(M)
#' operatorNorm(M)
#' @param M
#' 
#' @author Simone Vazzoler


l2norm <- function(M) {
  
  #s <- sqrt(sum(M^2))
  s <- sqrt(spectralRadius(t(M) %*% M))
  return(s)
  
}

l1norm <- function(M) {
  
  c <- max(colSums(abs(M)))
  return(c)
  
}

maxNorm <- function(M) {
  
  return(max(abs(M)))
  
}

frobNorm <- function(M) {
  
  A <- (t(M) %*% M)
  A <-  A * diag(nrow(A))

  return(sqrt(sum(A)))

}

spectralRadius <- function(M) {
  
  e <- eigen(M)
  maxEig <- max(abs(e$values))
  
  return(maxEig)
  
}

operatorNorm <- function(M) {
  
  return(sqrt(spectralRadius(t(M) %*% M)))
  
}

snr <- function(sig, err) {
  
  nc <- ncol(sig)
  SNR <- numeric(nc)
  
  s <- rowMeans(sig)
  e <- rowMeans(err)
  
  return (var(s)/var(e))
  
  for(i in 1:nc) {
    SNR[i] <- var(sig[,i])/var(err[,i])
  }

  return(SNR)
  
  # M1 <- sig %*% t(sig)
  # M2 <- err %*% t(err)
  # m1 <- Mod((eigen(M1)$values)[1])
  # m2 <- Mod((eigen(M2)$values)[1])
  # return(sqrt(m1/m2))
  
}