#' @title Multivariate VAR estimation
#'
#' @description A function to estimate a (possibly high-dimensional)
#' multivariate VAR time series using penalized least squares methods,
#' such as ENET, SCAD or MC+.
#'
#' @usage fitVAR(data, p = 1, penalty = "ENET", method = "cv", ...)
#'
#' @param data the data from the time series: variables in columns and
#' observations in rows
#' @param p order of the VAR model
#' @param penalty the penalty function to use. Possible values
#' are \code{"ENET"}, \code{"SCAD"} or \code{"MCP"}
#' @param method possible values are \code{"cv"} or \code{"timeSlice"}
#' @param ... the options for the estimation. Global options are:
#' \code{threshold}: if \code{TRUE} all the entries smaller than the oracle
#' threshold are set to zero;
#' \code{scale}: scale the data (default = FALSE)?
#' \code{nfolds}: the number of folds used for cross validation (default = 10);
#' \code{parallel}: if \code{TRUE} use multicore backend (default = FALSE);
#' \code{ncores}: if \code{parallel} is \code{TRUE}, specify the number
#' of cores to use for parallel evaluation. Options for ENET estimation:
#' \code{alpha}: the value of alpha to use in elastic net
#' (0 is Ridge regression, 1 is LASSO (default));
#' \code{type.measure}: the measure to use for error evaluation
#' (\code{"mse"} or \code{"mae"});
#' \code{nlambda}: the number of lambdas to use in the cross
#' validation (default = 100);
#' \code{leaveOut}: in the time slice validation leave out the
#' last \code{leaveOutLast} observations (default = 15);
#' \code{horizon}: the horizon to use for estimating mse/mae (default = 1);
#' \code{picasso}: use picasso package for estimation (only available
#' for \code{penalty = "SCAD"} and \code{method = "timeSlice"}).
#'
#' @return \code{A} the list (of length \code{p}) of the estimated matrices
#' of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#' @return \code{residuals} the time series of the residuals
#'
#' @export
fitVAR <- function(data, p = 1, penalty = "ENET", method = "cv", ...) {
  opt <- list(...)

  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  cnames <- colnames(data)

  if (method == "cv") {

    # use CV to find lambda
    opt$method <- "cv"
    out <- cvVAR(data, p, penalty, opt)
  } else if (method == "timeSlice") {

    # use timeslice to find lambda
    opt$method <- "timeSlice"
    out <- timeSliceVAR(data, p, penalty, opt)
  } else {

    # error: unknown method
    stop("Unknown method. Possible values are \"cv\" or \"timeSlice\"")
  }

  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in 1:length(out$A)) {
      colnames(out$A[[k]]) <- cnames
      rownames(out$A[[k]]) <- cnames
    }
  }

  return(out)
}

cvVAR <- function(data, p, penalty = "ENET", opt = NULL) {
  nc <- ncol(data)
  nr <- nrow(data)

  picasso <- ifelse(!is.null(opt$picasso), opt$picasso, FALSE)
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)

  thresholdType <- ifelse(!is.null(opt$thresholdType),
    opt$thresholdType, "soft"
  )

  returnFit <- ifelse(!is.null(opt$returnFit), opt$returnFit, FALSE)
  methodCov <- ifelse(!is.null(opt$methodCov), opt$methodCov, "tiger")

  if (picasso) {
    stop("picasso available only with timeSlice method.")
  }
  # transform the dataset
  trDt <- transformData(data, p, opt)

  if (penalty == "ENET") {

    # fit the ENET model
    t <- Sys.time()
    fit <- cvVAR_ENET(trDt$X, trDt$y, nvar = nc, opt)
    elapsed <- Sys.time() - t

    # extract what is needed
    lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)

    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = lambda)
    A <- matrix(Avector[2:length(Avector)],
      nrow = nc, ncol = nc * p,
      byrow = TRUE
    )

    mse <- min(fit$cvm)
  } else if (penalty == "SCAD") {

    # convert from sparse matrix to std matrix (SCAD does not work with sparse
    # matrices)
    trDt$X <- as.matrix(trDt$X)

    # fit the SCAD model
    t <- Sys.time()
    fit <- cvVAR_SCAD(trDt$X, trDt$y, opt)
    elapsed <- Sys.time() - t

    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)],
      nrow = nc, ncol = nc * p,
      byrow = TRUE
    )
    mse <- min(fit$cve)
  } else if (penalty == "MCP") {

    # convert from sparse matrix to std matrix (MCP does not work with sparse
    # matrices)
    trDt$X <- as.matrix(trDt$X)

    # fit the MCP model
    t <- Sys.time()
    fit <- cvVAR_SCAD(trDt$X, trDt$y, opt)
    elapsed <- Sys.time() - t

    # extract the coefficients and reshape the matrix
    Avector <- stats::coef(fit, s = "lambda.min")
    A <- matrix(Avector[2:length(Avector)],
      nrow = nc, ncol = nc * p,
      byrow = TRUE
    )
    mse <- min(fit$cve)
  } else {

    # Unknown penalty error
    stop("Unkown penalty. Available penalties are: ENET, SCAD, MCP.")
  }

  # If threshold = TRUE then set to zero all the entries that are smaller than
  # the threshold
  if (threshold == TRUE) {
    A <- applyThreshold(A, nr, nc, p, type = thresholdType)
  }

  # The full matrix A
  fullA <- A

  # Get back the list of VAR matrices (of length p)
  A <- splitMatrix(A, p)

  # Now that we have the matrices compute the residuals
  res <- computeResiduals(trDt$series, A)

  # To extract the sd of mse
  if (penalty == "ENET") {
    ix <- which(fit$cvm == min(fit$cvm))
    mseSD <- fit$cvsd[ix]
  } else {
    ix <- which(fit$cve == min(fit$cve))
    mseSD <- fit$cvse[ix]
  }

  # Create the output
  output <- list()
  output$mu <- trDt$mu
  output$A <- A

  # Do you want the fit?
  if (returnFit == TRUE) {
    output$fit <- fit
  }

  # Return the "best" lambda
  output$lambda <- fit$lambda.min

  output$mse <- mse
  output$mseSD <- mseSD
  output$time <- elapsed
  output$series <- trDt$series
  output$residuals <- res

  # Variance/Covariance estimation
  output$sigma <- estimateCovariance(res)

  output$penalty <- penalty
  output$method <- "cv"
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  return(output)
}

cvVAR_ENET <- function(X, y, nvar, opt) {
  a <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  nl <- ifelse(is.null(opt$nlambda), 100, opt$nlambda)
  tm <- ifelse(is.null(opt$type.measure), "mse", opt$type.measure)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)

  # Vector of lambdas to work on
  if (!is.null(opt$lambdas_list)) {
    lambdas_list <- opt$lambdas_list
  } else {
    lambdas_list <- c(0)
  }

  # Assign ids to the CV-folds (useful for replication of results)
  if (is.null(opt$foldsIDs)) {
    foldsIDs <- numeric(0)
  } else {
    nr <- nrow(X)
    foldsIDs <- rep(sort(rep(seq(nf), length.out = nr / nvar)), nvar)
  }

  if (parall == TRUE) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {

      # cl <- doMC::registerDoMC(cores = ncores)
      # using doMC as in glmnet vignettes
      cl <- doParallel::registerDoParallel(cores = ncores)

      if (length(foldsIDs) == 0) {
        if (length(lambdas_list) < 2) {
          cvfit <- glmnet::cv.glmnet(X, y,
            alpha = a, nlambda = nl,
            type.measure = tm, nfolds = nf,
            parallel = TRUE, standardize = FALSE
          )
        } else {
          cvfit <- glmnet::cv.glmnet(X, y,
            alpha = a, lambda = lambdas_list,
            type.measure = tm, nfolds = nf,
            parallel = TRUE, standardize = FALSE
          )
        }
      } else {
        if (length(lambdas_list) < 2) {
          cvfit <- glmnet::cv.glmnet(X, y,
            alpha = a, nlambda = nl,
            type.measure = tm, foldid = foldsIDs,
            parallel = TRUE, standardize = FALSE
          )
        } else {
          cvfit <- glmnet::cv.glmnet(X, y,
            alpha = a, lambda = lambdas_list,
            type.measure = tm, foldid = foldsIDs,
            parallel = TRUE, standardize = FALSE
          )
        }
      }
    }
  } else {
    if (length(foldsIDs) == 0) {
      if (length(lambdas_list) < 2) {
        cvfit <- glmnet::cv.glmnet(X, y,
          alpha = a, nlambda = nl,
          type.measure = tm, nfolds = nf,
          parallel = FALSE, standardize = FALSE
        )
      } else {
        cvfit <- glmnet::cv.glmnet(X, y,
          alpha = a, lambda = lambdas_list,
          type.measure = tm, nfolds = nf,
          parallel = FALSE, standardize = FALSE
        )
      }
    } else {
      if (length(lambdas_list) < 2) {
        cvfit <- glmnet::cv.glmnet(X, y,
          alpha = a, nlambda = nl,
          type.measure = tm, foldid = foldsIDs,
          parallel = FALSE, standardize = FALSE
        )
      } else {
        cvfit <- glmnet::cv.glmnet(X, y,
          alpha = a, lambda = lambdas_list,
          type.measure = tm, foldid = foldsIDs,
          parallel = FALSE, standardize = FALSE
        )
      }
    }
  }

  return(cvfit)
}

cvVAR_SCAD <- function(X, y, opt) {
  e <- ifelse(is.null(opt$eps), 0.01, opt$eps)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)
  picasso <- ifelse(is.null(opt$picasso), FALSE, TRUE)

  if (!picasso) {
    if (parall == TRUE) {
      if (ncores < 1) {
        stop("The number of cores must be > 1")
      } else {
        cl <- parallel::makeCluster(ncores)
        cvfit <- ncvreg::cv.ncvreg(X, y,
          nfolds = nf, penalty = "SCAD",
          eps = e, cluster = cl
        )
        parallel::stopCluster(cl)
      }
    } else {
      cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "SCAD", eps = e)
    }
  } else {
    cvfit <- picasso::picasso(X, y, method = "scad")
  }

  return(cvfit)
}

cvVAR_MCP <- function(X, y, opt) {
  e <- ifelse(is.null(opt$eps), 0.01, opt$eps)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)

  if (parall == TRUE) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(X, y,
        nfolds = nf, penalty = "MCP",
        eps = e, cluster = cl
      )
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(X, y, nfolds = nf, penalty = "MCP", eps = e)
  }

  return(cvfit)
}
