#' @title Multivariate VAR estimation
#'
#' @description A function to estimate a (possibly high-dimensional)
#' multivariate VAR time series using penalized least squares methods,
#' such as ENET, SCAD or MC+.
#'
#' @usage fit_var(data, p = 1, penalty = "ENET", method = "cv", ...)
#'
#' @param data the data from the time series: variables in columns and
#' observations in rows
#' @param p order of the VAR model
#' @param penalty the penalty function to use. Possible values
#' are \code{"ENET"}, \code{"SCAD"} or \code{"MCP"}
#' @param method possible values are \code{"cv"} or \code{"timeSlice"}
#' @param ... the options for the estimation. Unrecognised option names are an
#' error, not silently ignored. Global options are:
#' \code{threshold}: \code{FALSE} (default, no thresholding), \code{TRUE}
#' (threshold at the cutoff named by \code{threshold_type}), or a single
#' positive number used directly as the cutoff, in which case
#' \code{threshold_type} is ignored;
#' \code{threshold_type}: \code{"theory"} (default), \code{"soft"} or
#' \code{"hard"} -- see Details;
#' \code{threshold_const}: the constant \code{C} in the \code{"theory"} cutoff
#' (default = 1), ignored by the other two;
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
#' \code{horizon}: the forecast horizon used to score the time-slice
#' validation. Fixed at 1 in this version; other values are not yet
#' implemented and are rejected.
#'
#' @details
#' \strong{Thresholding.} Penalized estimates shrink coefficients but rarely
#' set enough of them to exactly zero, so a post-hoc cutoff is applied to the
#' estimated transition matrices. \code{threshold_type} selects the
#' \emph{cutoff}; the \emph{rule} is always hard thresholding,
#' \eqn{a \cdot 1(|a| \ge \tau)}. The names \code{"soft"} and \code{"hard"} are
#' therefore a misnomer -- they are two cutoffs under one rule -- and are kept
#' only for backward compatibility.
#'
#' \itemize{
#'   \item \code{"theory"} (default):
#'     \eqn{\tau_{ij} = (\sigma_i / s_j) \, C \sqrt{\log(n^2 p) / (T - p)}},
#'     where \eqn{\sigma_i} is the residual standard deviation of equation
#'     \eqn{i} taken from the un-thresholded fit, \eqn{s_j} is the standard
#'     deviation of series \eqn{j} (shared across all \eqn{p} lags), \eqn{n} is
#'     the number of series, \eqn{T} the number of observations, and \eqn{C} is
#'     \code{threshold_const}. The \eqn{\sigma_i / s_j} factor makes the cutoff
#'     dimensionally coherent, so it is not distorted when series are measured
#'     in different units.
#'   \item \code{"soft"}: \eqn{\tau = 1 / (p n \log(T - p))}.
#'   \item \code{"hard"}: \eqn{\tau = n^{-0.49}}.
#' }
#'
#' Neither legacy cutoff depends on \eqn{T} and \eqn{p} in a way derived from
#' theory -- \code{"hard"} does not depend on \eqn{T} at all -- which is why
#' \code{"theory"} is the default from version 1.1.0.
#'
#' \strong{Choosing \code{threshold_const}.} The default \eqn{C = 1} was
#' calibrated on a simulation grid over penalties, dimensions and sample sizes.
#' The optimum is \strong{flat}: performance is essentially unchanged for
#' \eqn{C} anywhere in \eqn{[0.75, 1.5]}, so 1 should not be read as a sharply
#' identified value. There is a mild dependence on the penalty -- SCAD and MCP
#' fits do slightly better at lower \eqn{C} (around 0.5 and 0.75 respectively)
#' -- so users fitting those penalties who care about the last fraction of
#' accuracy have a reason to lower it. A single default is used because
#' lowering it globally would cost the default ENET penalty more than it gains
#' the others.
#'
#' \strong{What thresholding does and does not do.} The cutoff \strong{controls
#' over-selection}: it removes spurious non-zero entries. It does \strong{not}
#' recover the true support. The supporting theory bounds false positives only,
#' under no minimum-signal (beta-min) condition, so there is no guarantee about
#' false negatives -- and in practice the thresholded estimator tends to
#' \emph{under}-select, discarding some genuinely non-zero coefficients along
#' with the spurious ones.
#'
#' @references
#' Basu, S. and Michailidis, G. (2015). Regularized estimation in sparse
#' high-dimensional time series models. \emph{The Annals of Statistics},
#' 43(4), 1535-1567. Proposition 4.1 states a thresholded estimator for this
#' model class and motivates the rate used here.
#'
#' Zhou, S. (2009). Thresholding procedures for high dimensional variable
#' selection and statistical estimation. \emph{Advances in Neural Information
#' Processing Systems}, 22.
#'
#' Bickel, P. J., Ritov, Y. and Tsybakov, A. B. (2009). Simultaneous analysis
#' of Lasso and Dantzig selector. \emph{The Annals of Statistics}, 37(4),
#' 1705-1732.
#'
#' The \code{"theory"} cutoff follows the noise-level form of the last two
#' references, specialised to the VAR parameter count and with a constant
#' calibrated by simulation; it is not stated verbatim in any of them.
#'
#' @return \code{A} the list (of length \code{p}) of the estimated matrices
#' of the process
#' @return \code{fit} the results of the penalized LS estimation
#' @return \code{mse} the mean square error of the cross validation
#' @return \code{time} elapsed time for the estimation
#' @return \code{residuals} the time series of the residuals
#'
#' @export
fit_var <- function(data, p = 1, penalty = "ENET", method = "cv", ...) {

  checkmate::assert_choice(penalty, c("ENET", "SCAD", "MCP"))
  checkmate::assert_count(p, positive = TRUE)
  checkmate::assert_choice(method, c("cv", "timeSlice"))

  opt <- list(...)
  validate_fit_options(opt)

  # convert data to matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  cnames <- colnames(data)

  opt$method <- method
  if (method == "cv") {
    # use CV to find lambda
    out <- cv_var(data, p, penalty, opt)
  } else {
    # use timeslice to find lambda
    out <- time_slice_var(data, p, penalty, opt)
  }

  # Add the names of the variables to the matrices
  if (!is.null(cnames)) {
    for (k in seq_along(out$A)) {
      colnames(out$A[[k]]) <- cnames
      rownames(out$A[[k]]) <- cnames
    }
  }

  out
}

cv_var <- function(data, p, penalty = "ENET", opt = NULL) {
  nc <- ncol(data)

  thr <- resolve_threshold(opt)

  return_fit <- ifelse(!is.null(opt$return_fit), opt$return_fit, FALSE)

  # transform the dataset
  tr_dt <- transform_data(data, p, opt)

  if (penalty == "ENET") {
    # fit the ENET model
    t <- Sys.time()
    fit <- cv_var_enet(tr_dt$X, tr_dt$y, nvar = nc, opt)
    elapsed <- Sys.time() - t

    # extract what is needed
    lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)

    # extract the coefficients and reshape the matrix
    a_vector <- stats::coef(fit, s = lambda)
    a <- matrix(
      a_vector[2:length(a_vector)],
      nrow = nc,
      ncol = nc * p,
      byrow = TRUE
    )

    mse <- min(fit$cvm)
  } else if (penalty == "SCAD") {
    # convert from sparse matrix to std matrix (SCAD does not work with sparse
    # matrices)
    tr_dt$X <- as.matrix(tr_dt$X)

    # fit the SCAD model
    t <- Sys.time()
    fit <- cv_var_scad(tr_dt$X, tr_dt$y, opt)
    elapsed <- Sys.time() - t

    # extract the coefficients and reshape the matrix
    a_vector <- stats::coef(fit, s = "lambda.min")
    a <- matrix(
      a_vector[2:length(a_vector)],
      nrow = nc,
      ncol = nc * p,
      byrow = TRUE
    )
    mse <- min(fit$cve)
  } else if (penalty == "MCP") {
    # convert from sparse matrix to std matrix (MCP does not work with sparse
    # matrices)
    tr_dt$X <- as.matrix(tr_dt$X)

    # fit the MCP model
    t <- Sys.time()
    fit <- cv_var_mcp(tr_dt$X, tr_dt$y, opt)
    elapsed <- Sys.time() - t

    # extract the coefficients and reshape the matrix
    a_vector <- stats::coef(fit, s = "lambda.min")
    a <- matrix(
      a_vector[2:length(a_vector)],
      nrow = nc,
      ncol = nc * p,
      byrow = TRUE
    )
    mse <- min(fit$cve)
  } else {
    # Unknown penalty error
    stop("Unknown penalty. Available penalties are: ENET, SCAD, MCP.")
  }

  # Set to zero all the entries smaller than the cutoff. "theory" needs the
  # un-thresholded residuals, so this must happen before split_matrix().
  a <- threshold_coefficients(a, tr_dt$series, nc, p, thr)

  # Get back the list of VAR matrices (of length p)
  a <- split_matrix(a, p)

  # Now that we have the matrices compute the residuals
  res <- compute_residuals(tr_dt$series, a)

  # To extract the sd of mse
  if (penalty == "ENET") {
    ix <- which(fit$cvm == min(fit$cvm))
    mse_sd <- fit$cvsd[ix]
  } else {
    ix <- which(fit$cve == min(fit$cve))
    mse_sd <- fit$cvse[ix]
  }

  # Create the output
  output <- list()
  output$mu <- tr_dt$mu
  output$sd <- tr_dt$sd
  output$A <- a

  # Do you want the fit?
  if (return_fit) {
    output$fit <- fit
  }

  # Return the "best" lambda
  output$lambda <- fit$lambda.min

  output$mse <- mse
  output$mse_sd <- mse_sd
  output$time <- elapsed
  output$series <- tr_dt$series
  output$residuals <- res

  # Variance/Covariance estimation
  output$sigma <- estimate_covariance(res)

  output$penalty <- penalty
  output$method <- "cv"
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  output
}

cv_var_enet <- function(x, y, nvar, opt) {
  a <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  nl <- ifelse(is.null(opt$nlambda), 100, opt$nlambda)
  tm <- ifelse(is.null(opt$type.measure), "mse", opt$type.measure)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)

  # Vector of lambdas to work on
  lambdas_list <- c(0)
  if (!is.null(opt$lambdas_list)) {
    lambdas_list <- opt$lambdas_list
  }

  # Assign ids to the CV-folds (useful for replication of results)
  folds_ids <- numeric(0)
  if (!is.null(opt$folds_ids)) {
    nr <- nrow(x)
    folds_ids <- rep(sort(rep(seq(nf), length.out = nr / nvar)), nvar)
  }

  if (parall) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- doParallel::registerDoParallel(cores = ncores)
    }
  }
  if (length(folds_ids) == 0) {
    if (length(lambdas_list) < 2) {
      cvfit <- glmnet::cv.glmnet(
        x,
        y,
        alpha = a,
        nlambda = nl,
        type.measure = tm,
        nfolds = nf,
        parallel = parall,
        standardize = FALSE
      )
    } else {
      cvfit <- glmnet::cv.glmnet(
        x,
        y,
        alpha = a,
        lambda = lambdas_list,
        type.measure = tm,
        nfolds = nf,
        parallel = parall,
        standardize = FALSE
      )
    }
  } else {
    if (length(lambdas_list) < 2) {
      cvfit <- glmnet::cv.glmnet(
        x,
        y,
        alpha = a,
        nlambda = nl,
        type.measure = tm,
        foldid = folds_ids,
        parallel = parall,
        standardize = FALSE
      )
    } else {
      cvfit <- glmnet::cv.glmnet(
        x,
        y,
        alpha = a,
        lambda = lambdas_list,
        type.measure = tm,
        foldid = folds_ids,
        parallel = parall,
        standardize = FALSE
      )
    }
  }
  if (parall) {
    parallel::stopCluster(cl)
  }

  cvfit
}

cv_var_scad <- function(x, y, opt) {
  e <- ifelse(is.null(opt$eps), 0.01, opt$eps)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)

  if (parall) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(
        x,
        y,
        nfolds = nf,
        penalty = "SCAD",
        eps = e,
        cluster = cl
      )
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(x, y, nfolds = nf, penalty = "SCAD", eps = e)
  }
  cvfit
}

cv_var_mcp <- function(x, y, opt) {
  e <- ifelse(is.null(opt$eps), 0.01, opt$eps)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)

  if (parall) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- parallel::makeCluster(ncores)
      cvfit <- ncvreg::cv.ncvreg(
        x,
        y,
        nfolds = nf,
        penalty = "MCP",
        eps = e,
        cluster = cl
      )
      parallel::stopCluster(cl)
    }
  } else {
    cvfit <- ncvreg::cv.ncvreg(x, y, nfolds = nf, penalty = "MCP", eps = e)
  }

  cvfit
}
