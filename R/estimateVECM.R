#' @title Multivariate VECM estimation
#' 
#' @description A function to estimate a (possibly big) multivariate VECM time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#' @param \code{data} the data from the time series: variables in columns and observations in 
#' rows
#' @param \code{p} order of the VECM model 
#' @param \code{penalty} the penalty function to use. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}
#' @param \code{logScale} should the function consider the \code{log} of the inputs? By default
#' this is set to \code{TRUE} 
#' @param \code{options} options for the function (TODO: specify)
#' 
#' @return \code{Pi} the matrix \code{Pi} for the VECM model 
#' @return \code{G} the list (of length \code{p-1}) of the estimated matrices of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#'
#' @author Simone Vazzoler
#'
#' @export
#' 
estimateVECM <- function(data, p = 2, penalty = "ENET", logScale = TRUE, 
                         options = NULL) {
  
  nr <- nrow(data)
  nc <- ncol(data)
  
  # add one to the order of VAR
  # p <- p + 1
  
  # by default log-scale the data
  if (logScale == TRUE) {
    # for (i in 1:nc) {
    #   data[, i] <- log(data[, i])
    # }
    data <- log(data)
    data[is.na(data)] <- 0
    # data[is.infinite(data)] <- 0
  }
  
  resultsVAR <- estimateVAR(data, p = p, penalty = penalty, options = options)
  
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
  output$Pi <- Pi
  output$G <- G
  output$fit <- resultsVAR$fit
  output$mse <- resultsVAR$mse
  output$time <- resultsVAR$time
  
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