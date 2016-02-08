#' @title Multivariate VAR estimation
#' 
#' @description A wrapper to estimate a (possibly big) multivariate VAR time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#' @param rets the data from the time series: variables in columns and observations in 
#' rows
#' @param penalty the penalty function to use. Possible values are \code{"ENET"}, \code{"SCAD"}
#' or \code{"MCP"}
#' @param p order of the VAR model (only for \code{"ENET"} penalty)
#' @param options options for the function
#' 
#' @author Simone Vazzoler
#'
#' @export
#' 
estimateVAR <- function(rets, penalty = "ENET", p = 1, options = NULL){

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
  
  if (penalty == "ENET") {
    
    if (p == 1) {
      # vectorization for VAR
      y <- as.vector(tmpY)
      # Hadamard product for data
      I <- Diagonal(nc)
      X <- kronecker(I, tmpX)
    } else {
      # create the data matrix
      #tmpX <- createDataMatrix(tmpX, p)
      tmpX <- duplicateMatrix(tmpX, p)
      #tmpY <- createDataMatrix(tmpY, p)
      tmpY <- tmpY[(p+1):nrow(tmpY), ]
      
      # replicate tmpX p times
      y <- as.vector(tmpY)
      I <- Diagonal(nc)
      X <- kronecker(I, tmpX)
    }

    # fit the ENET model
    t <- Sys.time()
    fit <- varENET(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(options$lambda), "lambda.min", options$lambda)

    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = lambda)
    if (p == 1) {
      A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc, byrow = TRUE)
    } else {
      # A <- Avector
      A <- matrix(Avector[2:length(Avector)], ncol = nc*p, byrow = TRUE) 
      A <- splitMatrix(A, p)
    }
    
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

# createDataMatrix <- function(data, p) {
#   
#   nr <- nrow(data)
#   nc <- ncol(data)
# 
#   tmpX <- matrix(0, nrow = (nr-1) * p, ncol = nc)
#   ix <- matrix(1:(nr*p), ncol = p, nrow = nr, byrow = TRUE)
#   
#   for (i in p:nr) {
#     
#     tmpIx <- ix[i-p+1, ]
#     tmpX[tmpIx, ] <- data[((i-p+1):i), ]
#     
#   }
#   
#   return(tmpX)
#   
# }

splitMatrix <- function(M, p) {
  
  nr <- nrow(M)
  
  A <- list()
  
  for (i in 1:p) {
    
    ix <- ((i-1) * nr) + (1:nr)
    A[[p-(i-1)]] <- M[1:nr, ix]  
    
  }

  return(A)
}

duplicateMatrix <- function(M, p) {
  
  nr <- nrow(M)
  nc <- ncol(M)
  
  finalM <- M
  
  for (i in 1:(p-1)) {
    
    tmpM <- matrix(0, nrow = nr, ncol = nc)
    tmpM[(i+1):nr, ] <- M[1:(nr-i), ]
    finalM <- cbind(finalM, tmpM)
    
  }
  
  finalM <- finalM[(p+1):nr, ]
  return(finalM)
  
}
