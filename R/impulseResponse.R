#' @title Impulse Response Function
#' 
#' @description A function to estimate the Impulse Response Function of a given VAR.
#' @param dataVar the data in the for of a VAR
#' @param len length of the impulse response function
#' 
#' @return \code{irf} a 3d array containing the impulse response function.
#' 
#' @usage impulseResponse(dataVar, len = 20)
#' 
#' @export
impulseResponse <- function(dataVar, len = 20) {
  
  if (!checkIsVar(dataVar)) { 
    stop("Input dataVar must be a VAR object")
  }
  
  A <- dataVar$A
  
  # Numerical problems in the estimated variance covariance
  e <- eigen(dataVar$sigma)$values
  if (!is.null(e[e<0])) {
    P <- t(chol(dataVar$sigma, pivot = TRUE))
  } else {
    P <- t(chol(dataVar$sigma))
  }
  
  nr <- nrow(dataVar$A[[1]])
  bigA <- companionVAR(dataVar)
  
  irf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len+1))
  oirf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len+1))
  
  Atmp <- diag(nrow = nrow(bigA), ncol = ncol(bigA))
  
  irf[ , , 1] <- Atmp[1:nr, 1:nr] 
  oirf[ , , 1] <- Atmp[1:nr, 1:nr] %*% P
  
  for (k in 1:len) {
    Atmp <- Atmp %*% bigA
    irf[ , , (k+1)] <- as.matrix(Atmp[1:nr, 1:nr])
    oirf[ , , (k+1)] <- as.matrix(Atmp[1:nr, 1:nr] %*% P)
  }
  
  ## TODO: add cumulative response functions
  out <- list()
  out$irf <- irf
  out$oirf <- oirf
  out$cholP <- P
  attr(out, "class") <- "irf"
  
  return(out)
}

companionVAR <- function(v) {
  
  A <- v$A
  nc <- ncol(A[[1]])
  p <- length(A)
  if (p>1){
    bigA <- Matrix::Matrix(0, nrow = p*nc, ncol = p*nc, sparse = TRUE)
    for (k in 1:p) {
      ix <- ((k-1) * nc) + (1:nc)
      bigA[1:nc, ix] <- A[[k]]
    }
    
    ixR <- (nc+1):nrow(bigA)
    ixC <- 1:((p-1)*nc)
    bigA[ixR, ixC] <- diag(1, nrow = length(ixC), ncol = length(ixC))  
  } else {
    bigA <- Matrix::Matrix(A[[1]], sparse = TRUE)
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

errorBandsIRF <- function(v, irf, alpha = 0.05) {
  
  M <- 100       # number of repeats
  lambda <- v$lambda 
  p <- length(v$A)
  nc <- ncol(v$A[[1]]) 
  len <- dim(irf$irf)[3]
  
  irfs <- array(data = rep(0,len*nc^2*M), dim = c(nc,nc,len+1, M))
  oirfs <- array(data = rep(0,len*nc^2*M), dim = c(nc,nc,len+1, M))
  
  for (k in 1:M) {
    # create Xs and Ys (temp variables)
    o <- genVARfromRes(v)
    nr <- nrow(o)
    nc <- ncol(o)
    tmpX <- o[1:(nr-1), ]
    tmpY <- o[2:(nr), ]
    
    # create the data matrix
    tmpX <- sparsevar::duplicateMatrix(tmpX, p)
    tmpY <- tmpY[p:nrow(tmpY), ]
    
    y <- as.vector(tmpY)
    
    # Hadamard product for data
    I <- Matrix::Diagonal(nc)
    X <- kronecker(I, tmpX)
    
    fit <- glmnet::glmnet(X, y, lambda = lambda)
    
    Avector <- stats::coef(fit, s = lambda)
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    M <- cbind(diag(x = 1, nrow = (nc*(p-1)), ncol = (nc*(p-1))), matrix(0, nrow = (nc*(p-1)), ncol = nc))
    bigA <- rbind(A,M)  

    tmpRes <- getIRF(v, bigA, len, irf$cholP)
    irfs[,,,k] <- tmpRes$irf
    oirfs[,,,k] <- tmpRes$oirf
  
  }
  
  irfUB <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len))
  oirfUB <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len))
  
  irfLB <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len))
  oirfLB <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len))

  irfSD <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len))
  oirfSD <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len))
  
  a <- alpha/2
  
  for (i in 1:nc) {
    for (j in 1:nc) {
      for (k in 1:len) {
        iUB <- sd(irfs[i,j,k,])
        oiUB <- sd(oirfs[i,j,k,])
        iUB <- quantile(irfs[i,j,k,], probs = (1-a))
        oiUB <- quantile(oirfs[i,j,k,], probs = (1-a))
        irfUB[i,j,k] <- iUB
        oirfUB[i,j,k] <- oiUB
        iLB <- quantile(irfs[i,j,k,], probs = a)
        oiLB <- quantile(oirfs[i,j,k,], probs = a)
        irfLB[i,j,k] <- iLB
        oirfLB[i,j,k] <- oiLB
        iSD <- sd(irfs[i,j,k,])
        oiSD <- sd(oirfs[i,j,k,])
        irfSD[i,j,k] <- iSD
        oirfSD[i,j,k] <- oiSD
        
      }
    }
  }
  
  output <- list()
  output$irfUB <- irfUB
  output$oirfUB <- oirfUB
  
  output$irfLB <- irfLB
  output$oirfLB <- oirfLB
  
  output$irfSD <- irfSD
  output$oirfSD <- oirfSD
  
  return(output)
}

genVARfromRes <- function(v) {
  
  ## This function creates the simulated time series 
  r <- v$residuals
  s <- v$series
  A <- v$A
  N <- ncol(A[[1]])
  p <- length(A)
  t <- nrow(r)
  
  zt <- matrix(0, nrow = t, ncol = N)
  
  for (t0 in (p+1):t) {
    ix <- sample(1:t, 1)
    u <- r[ix, ]
    vv <- rep(0, N) + u    
    for (i in 1:p){
      ph <- A[[i]]
      vv <- vv + ph %*% s[(t0-i), ]
    }
    zt[t0, ] <- vv
  }
  zt <- zt[(p+1):t, ]
  
  return(zt)
  
}

getIRF <- function (v, bigA, len = 20, P) {
  
  nr <- nrow(v$A[[1]])
  
  irf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len+1))
  oirf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len+1))
  
  Atmp <- diag(nrow = nrow(bigA), ncol = ncol(bigA))
  
  irf[ , , 1] <- Atmp[1:nr, 1:nr] 
  oirf[ , , 1] <- Atmp[1:nr, 1:nr] %*% P
  
  for (k in 1:len) {
    Atmp <- Atmp %*% bigA
    irf[ , , (k+1)] <- as.matrix(Atmp[1:nr, 1:nr])
    oirf[ , , (k+1)] <- as.matrix(Atmp[1:nr, 1:nr] %*% P)
  }
  
  ## TODO: add cumulative response functions
  out <- list()
  out$irf <- irf
  out$oirf <- oirf
  attr(out, "class") <- "irf"
  
  return(out)
}
