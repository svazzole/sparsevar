#' @title Impulse Response Function
#' 
#' @description A function to estimate the Impulse Response Function of a given VAR.
#' @param dataVar the data in the for of a VAR
#' @param len length of the impulse response function
#' 
#' @return \code{irf} a matrix containing the impulse response function.
#' 
#' @usage impulseResponse(dataVar, len = 20)
#' 
#' @export
impulseResponse <- function(dataVar, len = 20) {
  
  A <- dataVar$A
  nr <- nrow(dataVar$A[[1]])
  
  bigA <- companionVAR(dataVar)
  
  irf <- matrix(0, ncol = len, nrow = nr^2)
  
  Atmp <- diag(nrow = nrow(bigA), ncol = ncol(bigA))
  irf[,1] <- as.vector(t(Atmp))[1:nr^2]
  
  for (k in 1:(len-1)) {
    
    Atmp <- Atmp %*% bigA
    irf[,(k+1)] <- as.vector(t(Atmp))[1:nr^2]
    
  }
  return(irf)
}

companionVAR <- function(v) {
  
  A <- v$A
  nc <- ncol(A[[1]])
  p <- length(A)
  if (p>1){
    bigA <- matrix(0, p*nc, p*nc)
    for (k in 1:p) {
      ix <- ((k-1) * nc) + (1:nc)
      bigA[1:nc, ix] <- A[[k]]
    }
    
    ixR <- (nc+1):nrow(bigA)
    ixC <- 1:((p-1)*nc)
    bigA[ixR, ixC] <- diag(1, nrow = length(ixC), ncol = length(ixC))  
  } else {
    bigA <- A[[1]]
  }
  
  return(bigA)
}