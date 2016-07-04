#' @title Transorm data  
#' 
#' @description Transform the input data
#' 
#' @usage transformData(data, p, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param opt a list containing the options
#' 
#' @export
transformData <- function(data, p, opt) {
  
  # get the number of rows and columns
  nr <- nrow(data)
  nc <- ncol(data)
  
  # make sure the data is in matrix format
  data <- as.matrix(data)
  
  # scale the matrix columns
  scale <- ifelse(is.null(opt$scale), FALSE, opt$scale)  
  # center the matrix columns (default)
  center <- ifelse(is.null(opt$center), TRUE, opt$center)
  
  if (scale == TRUE) {
    data <- apply(FUN = scale, X = data, MARGIN = 2)
  } 
  
  if (center == TRUE) {
    m <- colMeans(data)
    cm <- matrix(rep(m, nrow(data)), nrow = nrow(data), byrow = TRUE) 
    data <- data - cm  
  }
  
  # create Xs and Ys (temp variables)
  tmpX <- data[1:(nr-1), ]
  tmpY <- data[2:(nr), ]
  
  # create the data matrix
  tmpX <- duplicateMatrix(tmpX, p)
  tmpY <- tmpY[p:nrow(tmpY), ]
  
  y <- as.vector(tmpY)
  
  # Hadamard product for data
  I <- Matrix::Diagonal(nc)
  X <- kronecker(I, tmpX)
  
  output <- list()
  output$X <- X
  output$y <- y
  output$series <- data
  output$mu <- t(m)
  
  return(output)
  
}

#' @title VAR ENET  
#' 
#' @description Estimate VAR using ENET penalty
#' 
#' @usage varENET(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#' 
#' @export
varENET <- function(data, p, lambdas, opt) {
  
  ## Fit a VAR for a sequence of lambdas 
  nc <- ncol(data)
  nr <- nrow(data)
  
  # transform the dataset
  trDt <- transformData(data, p, opt)
  
  fit <- glmnet::glmnet(trDt$X, trDt$y, lambda = lambdas)
  
  return(fit)
  
}

#' @title VAR SCAD  
#' 
#' @description Estimate VAR using SCAD penalty
#' 
#' @usage varSCAD(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#' 
#' @export
varSCAD <- function(data, p, lambdas, opt) {
  
  ## Fit a VAR for a sequence of lambdas 
  nc <- ncol(data)
  nr <- nrow(data)
  
  # transform the dataset
  trDt <- transformData(data, p, opt)
  
  fit <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y, family = "gaussian", penalty = "SCAD",
                        alpha = 1, lambda = lambdas)
  
  return(fit)
  
}

#' @title VAR MCP  
#' 
#' @description Estimate VAR using MCP penalty
#' 
#' @usage varMCP(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#' 
#' @export
varMCP <- function(data, p, lambdas, opt) {
  
  ## Fit a VAR for a sequence of lambdas 
  nc <- ncol(data)
  nr <- nrow(data)
  
  # transform the dataset
  trDt <- transformData(data, p, opt)
  
  fit <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y, family = "gaussian", penalty = "MCP",
                        alpha = 1, lambda = lambdas)
  
  return(fit)
  
}

splitMatrix <- function(M, p) {
  
  nr <- nrow(M)
  A <- list()
  
  for (i in 1:p) {
    
    ix <- ((i-1) * nr) + (1:nr)
    A[[i]] <- M[1:nr, ix]  
    
  }
  
  return(A)
}

duplicateMatrix <- function(data, p) {
  
  nr <- nrow(data)
  nc <- ncol(data)
  
  outputData <- data
  
  if (p > 1) {
    for (i in 1:(p-1)) {
      
      tmpData <- matrix(0, nrow = nr, ncol = nc)
      tmpData[(i+1):nr, ] <- data[1:(nr-i), ]
      outputData <- cbind(outputData, tmpData)
      
    }
  }
  
  outputData <- outputData[p:nr, ]
  return(outputData)
  
}

computeResiduals <- function(data, A) {
  
  nr <- nrow(data)
  nc <- ncol(data)
  p <- length(A)
  
  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)
  
  for (i in 1:p) {
    
    tmpD <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data)-i), ])
    tmpF <- t(A[[i]] %*% t(tmpD))
    f <- f + tmpF  
    
  }
  
  res <- data - f
  return(res)
  
}

#' @title Companion VAR  
#' 
#' @description Build the VAR(1) representation of a VAR(p) process
#' 
#' @usage companionVAR(v)
#' 
#' @param v the VAR object as from fitVAR or simulateVAR
#'
#' @export
companionVAR <- function(v) {
  
  ## TODO: check that v is a var object
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

#' @title Bootstrap VAR  
#' 
#' @description Build the bootstrapped series from the original var
#' 
#' @usage bootstrappedVAR(v)
#' 
#' @param v the VAR object as from fitVAR or simulateVAR
#'
#' @export
bootstrappedVAR <- function(v) {

  ## This function creates the bootstrapped time series
  ## TODO: check that v is a var object
  r <- v$residuals
  s <- v$series
  A <- v$A
  N <- ncol(A[[1]])
  p <- length(A)
  t <- nrow(r)

  zt <- matrix(0, nrow = t, ncol = N)
  zt[1:p,] <- s[1:p,]

  for (t0 in (p+1):t) {
    ix <- sample(1:t, 1)
    u <- r[ix, ]
    vv <- rep(0, N)   
    for (i in 1:p){
      ph <- A[[i]]
      vv <- vv + ph %*% zt[(t0-i), ]
    }
    vv <- vv + u
    zt[t0, ] <- vv
  }

return(zt)

}

#' @title Test for Ganger Causality
#' 
#' @description This function should retain only the coefficients of the 
#' matrices of the VAR that are statistically significative (from the bootstrap) 
#' 
#' @usage testGranger(v, eb)
#' 
#' @param v the VAR object as from fitVAR or simulateVAR
#' @param eb the error bands as obtained from errorBands
#'
#' @export
testGranger <- function(v, eb) {
 
  p <- length(v$A)
  A <- list()
  for (i in 1:p) {
    L <- (eb$irfQUB[,,i+1] >= 0 & eb$irfQLB[,,i+1]<=0)
    A[[i]] <- v$A[[i]] * (1-L)
  }
 
  
  return(A)
  
}