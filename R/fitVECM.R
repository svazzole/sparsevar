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
fitVECM <- function(data, p = 1, penalty = "ENET", method = "cv", logScale = TRUE, ...) {
  
  nr <- nrow(data)
  nc <- ncol(data)
  
  p <- p + 1
  
  opt <- list(...)
  
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
  
  for (k in 1:(p-1)) {
    G[[k]] <- - matrixSum(M, ix = k+1)
  }
  
  output <- list()
  output$mu <- resultsVAR$mu
  output$Pi <- Pi
  output$G <- G
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