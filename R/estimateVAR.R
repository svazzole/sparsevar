#' @title Multivariate VAR estimation
#' 
#' @description A wrapper to estimate a (possibly big) multivariate VAR time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#' @param rets the data from the time series: variables in columns and observations in 
#' rows
#' @param penalty the penalty function to use. Possible values are \code{"ENET"}, \code{"SCAD"}
#' or \code{"MCP"}
#' @param options options for the function
#' 
#' @author Simone Vazzoler
#'
#' @export
estimateVAR <- function(rets, penalty = "ENET", options = NULL){

  # get the number of rows and columns
  nr <- nrow(rets)
  nc <- ncol(rets)
  
  # make sure the data is in matrix format
  rets <- as.matrix(rets)
  
  # scale the matrix columns
  for (i in 1:nc) {
    rets[, i] <- scale(rets[, i])
  }
  
  # create Xs and Ys (temp variables)
  tmpX <- rets[1:(nr-1), ]
  tmpY <- rets[2:(nr), ]
  
  # vectorization for VAR
  y <- as.vector(tmpY)
  
  if (penalty == "ENET") {

    # Hadamard product for data
    I <- Diagonal(nc)
    X <- kronecker(I, tmpX)

    # fit the ENET model
    t <- Sys.time()
    fit <- varENET(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(options$lambda), "lambda.min", options$lambda)

    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = lambda)
    A <- matrix(Avector[2:nrow(Avector)], nrow = nc, ncol = nc, byrow = TRUE)

    mse <- min(fit$cvm)
    
  } else if (penalty == "SCAD") {
    
    I <- diag(nc)
    X <- as.matrix(I %x% tmpX)
    
    # fit the SCAD model
    t <- Sys.time()
    fit <- varSCAD(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc, byrow = TRUE)
    mse <- min(fit$cve)
    
  } else if (penalty == "MCP") {
    
    I <- diag(nc)
    X <- as.matrix(I %x% tmpX)
    
    # fit the MCP model
    t <- Sys.time()
    fit <- varMCP(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc, byrow = TRUE)
    mse <- min(fit$cve)
    
  } else {
    
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
  
  }
  
  output = list()
  output$A <- A
  output$fit <- fit
  output$mse <- mse
  output$time <- elapsed
  
  return(output)
  
}

varENET <- function(X,y, options = NULL) {
  
  a  <- ifelse(is.null(options$alpha), 1, options$alpha)
  nl <- ifelse(is.null(options$nlambda), 100, options$nlambda)
  tm <- ifelse(is.null(options$type.measure), "mse", options$type.measure)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)

  if(parall == TRUE) {
    cl <- registerDoMC(ncores)
  }

  if(ncores < 1) {
    stop("The number of cores must be > 1")
  }
    
  cvfit = cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = parall)
  
  return(cvfit)
  
}

varSCAD <- function(X, y, options = NULL) {
  
  e <- ifelse(is.null(options$eps), 0.01, options$eps)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  
  cvfit = cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e)
  
  return(cvfit)
  
}

varMCP <- function(X, y, options = NULL) {
  
  e <- ifelse(is.null(options$eps), 0.01, options$eps)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)

  cvfit = cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e)
  
  return(cvfit)
  
}
