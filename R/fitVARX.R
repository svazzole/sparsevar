#' @title Multivariate VARX estimation
#'
#' @description A function to estimate a (possibly high-dimensional)
#' multivariate VARX time series using penalized least squares methods,
#' such as ENET, SCAD or MC+.
#'
#' @usage fit_varx(data, p = 1, xt_matrix, m = 1,
#' penalty = "ENET", method = "cv", ...)
#'
#' @param data the data from the time series: variables in columns and
#' observations in rows
#' @param p order of the VAR model
#' @param xt_matrix the exogenous variables
#' @param m order of the exogenous variables
#' @param penalty the penalty function to use. Possible values are
#' \code{"ENET"}, \code{"SCAD"} or \code{"MCP"}
#' @param method possible values are \code{"cv"} or \code{"timeSlice"}
#' @param ... the options for the estimation. Unrecognised option names are an
#' error, not silently ignored. Global options are:
#' \code{threshold}: \code{FALSE} (default), \code{TRUE}, or a single positive
#' number used directly as the cutoff;
#' \code{threshold_type}: \code{"theory"} (default), \code{"soft"} or
#' \code{"hard"};
#' \code{threshold_const}: the constant \code{C} in the \code{"theory"} cutoff
#' (default = 1). See \code{\link{fit_var}} for the formulas and for what
#' thresholding does and does not guarantee;
#' \code{scale}: scale the data (default = FALSE)?
#' \code{nfolds}: the number of folds used for cross validation (default = 10);
#' \code{parallel}: if \code{TRUE} use multicore backend (default = FALSE);
#' \code{ncores}: if \code{parallel} is \code{TRUE}, specify the number of
#' cores to use for parallel evaluation. Options for ENET estimation:
#' \code{alpha}: the value of alpha to use in elastic net (0 is Ridge
#' regression, 1 is LASSO (default));
#' \code{type.measure}: the measure to use for error evaluation (\code{"mse"}
#' or \code{"mae"});
#' \code{nlambda}: the number of lambdas to use in the cross validation
#' (default = 100);
#' \code{leaveOut}: in the time slice validation leave out the last
#' \code{leaveOutLast} observations (default = 15);
#' \code{horizon}: the forecast horizon used to score the time-slice
#' validation. Fixed at 1 in this version; other values are not yet
#' implemented and are rejected.
#'
#' @return \code{A} the list (of length \code{p}) of the estimated matrices
#' of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#' @return \code{residuals} the time series of the residuals
#'
#' @export
fit_varx <- function(data, p = 1,
                     xt_matrix, m = 1,
                     penalty = "ENET",
                     method = "cv", ...) {
  opt <- list(...)
  validate_fit_options(opt)

  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  # convert data to matrix
  if (!is.matrix(xt_matrix)) {
    xt_matrix <- as.matrix(xt_matrix)
  }

  data_xt <- cbind(data, xt_matrix)

  cnames <- colnames(data)

  px <- max(p, m)

  if (method == "cv") {
    # use CV to find lambda
    opt$method <- "cv"
    out <- cv_var(data_xt, px, penalty, opt)
  } else if (method == "timeSlice") {
    # use timeslice to find lambda
    opt$method <- "timeSlice"
    out <- time_slice_var(data_xt, px, penalty, opt)
  } else {
    # error: unknown method
    stop("Unknown method. Possible values are \"cv\" or \"timeSlice\"")
  }

  nc <- ncol(data)
  nc_x <- ncol(xt_matrix)

  out <- var_to_varx(out, p, m, nc, nc_x)

  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in seq_along(out$A)) {
      colnames(out$A[[k]]) <- cnames
      rownames(out$A[[k]]) <- cnames
    }
  }

  out
}

var_to_varx <- function(v, p, m, nc, nc_x) {
  l <- length(v$A)
  new_a <- list()
  b <- list()
  for (i in 1:l) {
    new_a[[i]] <- v$A[[i]][1:nc, 1:nc]
    b[[i]] <- v$A[[i]][1:nc, (nc + 1):ncol(v$A[[i]])]
  }
  if (p < l) {
    v$newA <- new_a[-((p + 1):l)]
    v$B <- b
  } else if (m < l) {
    v$newA <- new_a
    v$B <- b[-((m + 1):l)]
  } else {
    v$newA <- new_a
    v$B <- b
  }
  v
}
