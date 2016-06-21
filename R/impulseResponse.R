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
impulseResponse <- function(dataVar, len = 20, method = "non-orthogonal") {
  
  if (method == "orthogonal") {
    A <- dataVar$A
    nr <- nrow(dataVar$A[[1]])
    bigA <- companionVAR(dataVar)
    irf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len)) #matrix(0, ncol = len, nrow = nr^2)
    Atmp <- diag(nrow = nrow(bigA), ncol = ncol(bigA))
    irf[ , , 1] <- Atmp[1:nr, 1:nr]
    
    for (k in 1:(len-1)) {
      Atmp <- Atmp %*% bigA
      irf[ , , (k+1)] <- Atmp[1:nr, 1:nr]
    }
  } else if (method == "non-orthogonal") {
    stop("To be implemented")
  } else {
    stop("Unknown method. Possible values are: orthogonal or non-orthogonal")
  }
  
#  attr(irf, "class") <- "sparsevar"
#  attr(irf, "type") <- "irf"
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

#' @title Check Impulse Zero
#' 
#' @description A function to find which entries of the impulse response function 
#' are zero.
#' @param irf irf output from impulseResponse function
#' 
#' @return a matrix containing the indices of the impulse response function that
#' are 0.
#' 
#' @usage checkImpulseZero(irf)
#' 
#' @export
checkImpulseZero <- function(irf) {
  
  nx <- dim(irf)[1]
  ny <- dim(irf)[2]
  nz <- dim(irf)[3]
  
  logicalIrf <- matrix(0, nx, ny)
  
  for (z in 1:nz) {
    
    logicalIrf <- logicalIrf + abs(irf[ , , z])
    
  }
  
  logicalIrf <- logicalIrf == 0
  
  return(which(logicalIrf == TRUE, arr.ind = TRUE))
}

errorBandsIRF <- function(d, v, irf) {
  
  r <- v$residuals
  s <- cov(r) # variance / covariance of the residuals
  
  nr <- nrow(d)
  
  Q_T <- 1/nr * t(t(d)%*%d)
  P <- kronecker(s, solve(Q_T))
  
  nz <- dim(irf)[3]
  
  for (i in 1:nz){
    
    
  }
  
  return(0)
}
