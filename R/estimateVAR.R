#' @title Multivariate VAR estimation
#' 
#' @description A function to estimate a (possibly big) multivariate VAR time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#' @param \code{data} the data from the time series: variables in columns and observations in 
#' rows
#' @param \code{p} order of the VAR model
#' @param \code{penalty} the penalty function to use. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}
#' @param \code{options} options for the function (TODO: specify)
#' 
#' @return \code{A} the list (of length \code{p}) of the estimated matrices of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#'
#' @author Simone Vazzoler
#'
#' @export
#' 
estimateVAR <- function(data, p = 1, penalty = "ENET", options = NULL) {

  # get the number of rows and columns
  nr <- nrow(data)
  nc <- ncol(data)
  
  # make sure the data is in matrix format
  data <- as.matrix(data)
  
  # scale the matrix columns
  # for (i in 1:nc) {
  #   data[, i] <- scale(data[, i])
  # }
  
  # create Xs and Ys (temp variables)
  tmpX <- data[1:(nr-1), ]
  tmpY <- data[2:(nr), ]
  
  # create the data matrix
  tmpX <- duplicateMatrix(tmpX, p)
  tmpY <- tmpY[p:nrow(tmpY), ]

  y <- as.vector(tmpY)
  
  # Hadamard product for data
  I <- Diagonal(nc)
  X <- kronecker(I, tmpX)
  
  if (penalty == "ENET") {
    
    # fit the ENET model
    t <- Sys.time()
    fit <- varENET(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract what is needed
    lambda <- ifelse(is.null(options$lambda), "lambda.min", options$lambda)

    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = lambda)

    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    
    if (!is.null(options$threshold)){
      if (options$threshold == TRUE) {
        tr <- 1 / sqrt(p*nc)
        L <- abs(A) >= tr
        A <- A * L
      }
    }
    
    A <- splitMatrix(A, p)
    mse <- min(fit$cvm)
    
  } else if (penalty == "SCAD") {
    
    # convert from sparse matrix to std matrix
    X <- as.matrix(X)
    
    # fit the SCAD model
    t <- Sys.time()
    fit <- varSCAD(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    A <- splitMatrix(A, p)
    
    mse <- min(fit$cve)
    
  } else if (penalty == "MCP") {
    
    # convert from sparse matrix to std matrix
    X <- as.matrix(X)
    
    # fit the MCP model
    t <- Sys.time()
    fit <- varMCP(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    A <- splitMatrix(A, p)
    
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
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      #cl <- registerDoMC(ncores)
      cl <- parallel::makeCluster(ncores)
      cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = TRUE)
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = FALSE)
  }
  
  return(cvfit)
  
}

varSCAD <- function(X, y, options = NULL) {
  
  e <- ifelse(is.null(options$eps), 0.01, options$eps)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)
  
  if(parall == TRUE) {
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e, cluster = cl)
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e)
  }

  return(cvfit)
  
}

varMCP <- function(X, y, options = NULL) {
  
  e <- ifelse(is.null(options$eps), 0.01, options$eps)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)
  
  if(parall == TRUE) {
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e, cluster = cl)
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e)
  }
  
  return(cvfit)
  
}


splitMatrix <- function(M, p) {
  
  nr <- nrow(M)
  
  A <- list()
  
  for (i in 1:p) {
    
    ix <- ((i-1) * nr) + (1:nr)
    # A[[p-(i-1)]] <- M[1:nr, ix]  
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
