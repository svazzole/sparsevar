#' @title Multivariate VECM estimation
#' 
#' @description A function to estimate a (possibly big) multivariate VECM time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#'  
#' @usage fitVECM(data, p, penalty, method, logScale, ...)
#' 
#' @param data the data from the time series: variables in columns and observations in 
#' rows
#' @param p order of the VECM model 
#' @param penalty the penalty function to use. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}
#' @param logScale should the function consider the \code{log} of the inputs? By default
#' this is set to \code{TRUE} 
#' @param method \code{"cv"} or \code{"timeSlice"}
#' @param ... options for the function (TODO: specify)
#' 
#' @return Pi the matrix \code{Pi} for the VECM model 
#' @return G the list (of length \code{p-1}) of the estimated matrices of the process
#' @return fit the results of the penalized LS estimation
#' @return mse the mean square error of the cross validation
#' @return time elapsed time for the estimation
#' 
#' @export
fitVECM <- function(data, p = 0, penalty = "ENET", method = "cv", logScale = TRUE, ...) {
  
  nr <- nrow(data)
  nc <- ncol(data)
  
  p <- p + 1
  
  opt <- list(...)
  opt$center <- FALSE
  
  # by default log-scale the data
  if (logScale == TRUE) {
    data <- log(data)
    data[is.na(data)] <- 0
    # data[is.infinite(data)] <- 0
  }
  
  resultsVAR <- fitVAR(data, p = p, penalty = penalty, method = method, ...)
  M <- resultsVAR$A
  I <- diag(x = 1, nrow = nc, ncol = nc)
  
  # Coint matrix
  Pi <- -(I - matrixSum(M, ix = 1))
  
  # Gamma matrices
  G <- list()
  
  if (p>1){
    for (k in 1:(p-1)) {
      G[[k]] <- - matrixSum(M, ix = k+1)
    }
  }
  
  output <- list()
  output$mu <- resultsVAR$mu
  output$Pi <- Pi
  output$G <- G
  output$A <- resultsVAR$A
  output$fit <- resultsVAR$fit
  output$mse <- resultsVAR$mse
  output$mseSD <- resultsVAR$mseSD
  output$time <- resultsVAR$time
  output$residuals <- resultsVAR$residuals
  output$lambda <- resultsVAR$lambda
  output$series <- resultsVAR$series
  
  if (is.null(opt$methodCov)) {
    output$sigma <- estimateCovariance(output$residuals)
  } else {
    output$sigma <- estimateCovariance(output$residuals, methodCovariance = opt$methodCov)
  }
  
  output$penalty <- resultsVAR$penalty
  output$method <- resultsVAR$method
  attr(output, "class") <- "vecm"
  attr(output, "type") <- "fit"
  
  return(output)
  
}

matrixSum <- function(M, ix = 1) {
  
  l <- length(M)
  nc <- ncol(M[[1]])
  
  A <- matrix(0, nrow = nc, ncol = nc)
  
  for (i in ix:l) {
    A = A + M[[i]]  
  }
  
  return(A)
  
}

#' @export
decomposePi <- function(vecm, rk) {
  
  if(attr(vecm, "class")!="vecm") {
    stop("The input is not a vecm object.")
  }
  
  nc <- ncol(vecm$Pi)
  Pi <- vecm$Pi
  colnames(Pi) <- NULL
  rownames(Pi) <- NULL
  sig <- corpcor::invcov.shrink(vecm$residuals, verbose = FALSE)
  colnames(sig) <- NULL
  rownames(sig) <- NULL
  
  if(rk >=1) {
    a <- Pi[,1:rk]
    b <- t(solve(t(a)%*%sig%*%a)%*%(t(a)%*%sig%*%Pi[,(rk+1):nc]))
    b <- rbind(diag(1, rk, rk), b)
  } else {
    a <- numeric(0,length = nc)
    b <- Pi
  }

  out <- list()
  out$alpha <- a
  out$beta <- b
  return(out)
  
}

#' @export
decomposePi2 <- function(vecm, rk) {
  
  if(attr(vecm, "class")!="vecm") {
    stop("The input is not a vecm object.")
  }
  
  nc <- ncol(vecm$Pi)
  Pi <- vecm$Pi
  colnames(Pi) <- NULL
  rownames(Pi) <- NULL
  
  if(rk >=1) {
    a <- Pi[,1:rk]
    # s <- solve(vecm$sigma)
    # b <- t(solve(t(a)%*%s%*%a)%*%(t(a)%*%s%*%vecm$Pi[,(rk+1):nc]))
    # b <- rbind(diag(1, rk, rk), b)
    A <- kronecker(diag(1,nc,nc), a)
    B <- as.numeric(Pi)
    b <- matrix(qr.solve(A,B), ncol = rk, nrow = nc, byrow = TRUE)
    bT <- matrix(qr.solve(A,B), ncol = rk, nrow = nc, byrow = FALSE)
  } else {
    a <- numeric(0,length = nc)
    b <- Pi
    bT <- t(Pi)
  }
  
  out <- list()
  out$alpha <- a
  out$beta <- b
  out$betaT <- bT
  return(out)
  
}
