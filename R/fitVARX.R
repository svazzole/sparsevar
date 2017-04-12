#' @title Multivariate VARX estimation
#'
#' @description A function to estimate a (possibly high-dimensional) multivariate VARX time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#'
#' @usage fitVARX(data, p = 1, Xt, m = 1, penalty = "ENET", method = "cv", ...)
#'
#' @param data the data from the time series: variables in columns and observations in
#' rows
#' @param p order of the VAR model
#' @param Xt the exogenous variables
#' @param m order of the exogenous variables
#' @param penalty the penalty function to use. Possible values are \code{"ENET"},
#' \code{"SCAD"} or \code{"MCP"}
#' @param method possible values are \code{"cv"} or \code{"timeSlice"}
#' @param ... the options for the estimation. Global options are:
#' \code{threshold}: if \code{TRUE} all the entries smaller than the oracle threshold are set to zero;
#' \code{scale}: scale the data (default = FALSE)?
#' \code{nfolds}: the number of folds used for cross validation (default = 10);
#' \code{parallel}: if \code{TRUE} use multicore backend (default = FALSE);
#' \code{ncores}: if \code{parallel} is \code{TRUE}, specify the number of cores to use
#' for parallel evaluation. Options for ENET estimation:
#' \code{alpha}: the value of alpha to use in elastic net (0 is Ridge regression, 1 is LASSO (default));
#' \code{type.measure}: the measure to use for error evaluation (\code{"mse"} or \code{"mae"});
#' \code{nlambda}: the number of lambdas to use in the cross validation (default = 100);
#' \code{leaveOut}: in the time slice validation leave out the last \code{leaveOutLast} observations
#' (default = 15);
#' \code{horizon}: the horizon to use for estimating mse/mae (default = 1);
#' \code{picasso}: use picasso package for estimation (only available for \code{penalty = "SCAD"} 
#' and \code{method = "timeSlice"}).
#'
#' @return \code{A} the list (of length \code{p}) of the estimated matrices of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#' @return \code{residuals} the time series of the residuals
#'
#' @export
fitVARX <- function(data, p = 1, Xt, m = 1, penalty = "ENET", method = "cv", ...) {
  
  opt <- list(...)
  
  # convert data to matrix 
  if (!is.matrix(data)){
    data <- as.matrix(data)
  }
  
  # convert data to matrix 
  if (!is.matrix(Xt)){
    Xt <- as.matrix(Xt)
  }
  
  dataXt <- cbind(data, Xt)
  
  cnames <- colnames(data)
  cnamesX <- colnames(Xt)
  
  pX <- max(p,m)
  
  if (method == "cv") {
    # use CV to find lambda
    opt$method <- "cv"
    out <- cvVAR(dataXt, pX, penalty, opt)
  } else if (method == "timeSlice") {
    # use timeslice to find lambda
    opt$method <- "timeSlice"
    out <- timeSliceVAR(dataXt, pX, penalty, opt)
  } else {
    # error: unknown method
    stop("Unknown method. Possible values are \"cv\" or \"timeSlice\"")
  }
  
  nc <- ncol(data)
  ncX <- ncol(Xt)
  
  out <- VARtoVARX(out, p, m, nc, ncX)
  
  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in 1:length(out$A)) {
      colnames(out$A[[k]]) <- cnames
      rownames(out$A[[k]]) <- cnames
    }
  }
  
  return(out)
  
}

VARtoVARX <- function(v, p, m, nc, ncX) {
  l <- length(v$A)
  newA <- list()
  B <- list()
  for(i in 1:l) {
    newA[[i]] <- v$A[[i]][1:nc,1:nc]
    B[[i]] <- v$A[[i]][1:nc, (nc+1):ncol(v$A[[i]])]
  }
  if (p<l) {
    v$newA <- newA[-((p+1):l)]
    v$B <- B
  } else if (m<l) {
    v$newA <- newA
    v$B <- B[-((m+1):l)]
  } else {
    v$newA <- newA
    v$B <- B
  }
  return(v)
}