#' @title Multivariate VAR estimation
#' 
#' @description A function to estimate a (possibly high-dimensional) multivariate VAR time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#' @param data the data from the time series: variables in columns and observations in 
#' rows
#' @param p order of the VAR model
#' @param penalty the penalty function to use. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}
#' @param options a list containing the options for the estimation. Global options are:
#' \code{threshold}: if \code{TRUE} all the entries smaller than the oracle threshold are set to zero; 
#' \code{scale}: scale the data (default = FALSE)?
#' \code{nfolds}: the number of folds used for cross validation (default = 10);
#' \code{parallel}: if \code{TRUE} use multicore backend (default = FALSE);
#' \code{ncores}: if \code{parallel} is \code{TRUE}, specify the number of cores to use
#' for parallel evaluation. Options for ENET estimation: 
#' \code{alpha}: the value of alpha to use in elastic net (0 is Ridge regression, 1 is LASSO (default));
#' \code{type.measure}: the measure to use for error evaluation (\code{"mse"} or \code{"mae"});
#' \code{nlambda}: the number of lambdas to use in the cross validation (default = 100);
#' \code{repeatedCV}: use repeated cross validation (default = FALSE);
#' \code{nRepeats}: the number of repeats in the repeated cross validation (default = 3);
#' \code{timeSlice}: use time slice validation (default = FALSE);
#' \code{leaveOutLast}: in the time slice validation leave out the last \code{leaveOutLast} observations
#' (default = 15);
#' \code{horizon}: the horizon to use for estimating mse/mae (default = 1);
#' \code{fixedWindow}: use fixed window (default) or expanding window (FALSE).
#' 
#' @return \code{A} the list (of length \code{p}) of the estimated matrices of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#' @return \code{residuals} the time series of the residuals 
#' 
#' @usage estimateVAR(data, p = 1, penalty = "ENET", options = NULL)
#' 
#' @export
estimateVAR <- function(data, p = 1, penalty = "ENET", options = NULL) {

  # get the number of rows and columns
  nr <- nrow(data)
  nc <- ncol(data)
  
  # make sure the data is in matrix format
  data <- as.matrix(data)

  # scale the matrix columns
  scale <- ifelse(is.null(options$scale), FALSE, options$scale)  
  if (scale == TRUE) {
    data <- apply(FUN = scale, X = data, MARGIN = 2)
  } else {
    # only center the data
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
  
  if (penalty == "ENET") {
    
    # By default repeatedCV = FALSE
    options$repeatedCV <- ifelse(is.null(options$repeatedCV), FALSE, TRUE)
    options$timeSlice <- ifelse(is.null(options$timeSlice), FALSE, TRUE)
    
    # fit the ENET model
    t <- Sys.time()
    fit <- varENET(X, y, options)
    elapsed <- Sys.time() - t
    
    if ((options$repeatedCV == FALSE) && (options$timeSlice == FALSE)){
  
      # extract what is needed
      lambda <- ifelse(is.null(options$lambda), "lambda.min", options$lambda)
      
      # extract the coefficients and reshape the matrix
      Avector <- stats::coef(fit, s = lambda)
      A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    
    } else {
      
      Avector <- fit$Avector
      A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    
    }
    
    mse <- min(fit$cvm)
    
  } else if (penalty == "SCAD") {
    
    # convert from sparse matrix to std matrix (SCAD does not work with sparse matrices)
    X <- as.matrix(X)
    
    # fit the SCAD model
    t <- Sys.time()
    fit <- varSCAD(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    mse <- min(fit$cve)
    
  } else if (penalty == "MCP") {
    
    # convert from sparse matrix to std matrix (MCP does not work with sparse matrices)
    X <- as.matrix(X)
    
    # fit the MCP model
    t <- Sys.time()
    fit <- varMCP(X, y, options)
    elapsed <- Sys.time() - t
    
    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    mse <- min(fit$cve)
    
  } else {
    
    # Unknown penalty error
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
    
  }
  
  # If threshold = TRUE then set to zero all the entries that are small
  if (!is.null(options$threshold)) {
    if (options$threshold == TRUE) {
      tr <- 1 / sqrt(p*nc*log(nr))
      L <- abs(A) >= tr
      A <- A * L
    }
  }
  
  # Get back the list of VAR matrices (of length p)
  A <- splitMatrix(A, p)
  
  # Now that we have the matrices compute the residuals
  res <- computeResiduals(data, A)
  
  # Create the output
  output = list()
  output$mu <- t(m)
  output$A <- A
  
  # Do you want the fit?
  if (!is.null(options$returnFit)) {
    if (options$returnFit == TRUE) {
      output$fit <- fit
    }
  }
  
  # If ENET is used, return the lambda 
  if (penalty == "ENET") {
    output$lambda <- fit$lambda.min
  }
  
  output$mse <- mse
  output$time <- elapsed
  output$series <- data
  output$residuals <- res
  output$sigma <- cov(res)
  attr(output, "class") <- "var"
  attr(output, "type") <- "estimate"
  return(output)
  
}

varENET <- function(X,y, options = NULL) {
  
  a  <- ifelse(is.null(options$alpha), 1, options$alpha)
  nl <- ifelse(is.null(options$nlambda), 100, options$nlambda)
  tm <- ifelse(is.null(options$type.measure), "mse", options$type.measure)
  nf <- ifelse(is.null(options$nfolds), 10, options$nfolds)
  parall <- ifelse(is.null(options$parallel), FALSE, options$parallel)
  ncores <- ifelse(is.null(options$ncores), 1, options$ncores)
  repeatedCV <- options$repeatedCV # ifelse(is.null(options$repeatedCV), FALSE, TRUE)
  nRepeats <- ifelse(is.null(options$nRepeats), 3, options$nRepeats)
  timeSlice <- options$timeSlice
  initialWindow <- ifelse(is.null(options$leaveOutLast), nrow(X) - 15, nrow(X) - options$leaveOutLast)
  horizon <- ifelse(is.null(options$horizon), 1, options$horizon)
  fixedWindow <- ifelse(is.null(options$fixedWindow), TRUE, options$fixedWindow)
  
  if (repeatedCV == TRUE) {
    
    # Suppress warnings... Is this correct?
    options(warn = -1)
    trCtrl <- caret::trainControl(method = "repeatedcv", number = nf, repeats = nRepeats, returnData = FALSE)
    lam <- glmnet::glmnet(X, y, alpha = a)$lambda
    gr <- expand.grid(.alpha = a, .lambda = lam)
    fit <- caret::train(x = X, y = y, method = "glmnet", trControl = trCtrl, tuneGrid = gr, metric = "RMSE")
    b <- stats::coef(fit$finalModel, fit$bestTune$lambda)
    cvm <- min(fit$results$RMSE)^2 # extract the MSE
    fit <- list()
    fit$Avector <- b
    fit$cvm <- cvm  
    return(fit)
    
  } 
  
  if (timeSlice == TRUE) {
    
    # Suppress warnings... Is this correct?
    options(warn = -1)
    #inWind <- nrow(X) - 20
    trCtrl <- caret::trainControl(method = "timeslice", returnData = FALSE,
                                  initialWindow = initialWindow, horizon = horizon, fixedWindow = fixedWindow)
    lam <- glmnet::glmnet(X, y, alpha = a)$lambda
    gr <- expand.grid(.alpha = a, .lambda = lam)
    fit <- caret::train(x = X, y = y, method = "glmnet", trControl = trCtrl, tuneGrid = gr, metric = "RMSE")
    b <- stats::coef(fit$finalModel, fit$bestTune$lambda)
    cvm <- min(fit$results$RMSE)^2 # extract the MSE
    fit <- list()
    fit$Avector <- b
    fit$cvm <- cvm  
    return(fit)
    
  } 
  
  # Assign ids to the CV-folds (useful for replication of results)  
  if (is.null(options$foldsIDs)) {
    foldsIDs <- numeric(0)
  } else {
    nr <- nrow(X)
    foldsIDs <- sort(rep(seq(nf), length = nr))
    # foldsIDs <- rep(seq(nf), length = nr)
  }
  
  if(parall == TRUE) {
    if(ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
        # cl <- doMC::registerDoMC(cores = ncores) # using doMC as in glmnet vignettes
        cl <- doParallel::registerDoParallel(cores = ncores)
      if (length(foldsIDs) == 0) {
        cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = TRUE)
      } else {
        cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, foldid = foldsIDs, parallel = TRUE)
      }
    }
  } else {
    if (length(foldsIDs) == 0) {
      cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, nfolds = nf, parallel = FALSE)
    } else {
      cvfit <- glmnet::cv.glmnet(X, y, alpha = a, nlambda = nl, type.measure = tm, foldid = foldsIDs, parallel = FALSE)
    }
    
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
    A[[i]] <- M[1:nr, ix]  
    
  }

  return(A)
}

#' @export
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
