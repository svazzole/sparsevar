# Every option the fit paths read from `...`. fit_var() and fit_varx() collect
# `...` into a list and previously acted on whatever happened to be spelled
# correctly, so `fit_var(d, p = 1, thersholdd_type = "hard")` ran happily and
# returned the default answer. Every misspelled option was a silent wrong
# result; this is the list that makes it loud.
known_fit_options <- c(
  # thresholding
  "threshold", "threshold_type", "threshold_const",
  # data transformation
  "scale", "center",
  # lambda selection
  "nfolds", "nlambda", "lambda", "lambdas_list", "folds_ids",
  "type.measure", "eps", "alpha",
  # time-slice validation
  "leaveOut", "horizon",
  # covariance, output and execution
  "methodCov", "return_fit", "returnFit", "parallel", "ncores", "verbose",
  # set internally by fit_var()/fit_varx(), not by the user
  "method"
)

# Reject unknown options, suggesting the nearest known name. Called once per
# entry point, before anything is fitted.
validate_fit_options <- function(opt) {
  nms <- names(opt)
  if (is.null(nms) && length(opt) > 0) {
    stop("All options passed through `...` must be named.", call. = FALSE)
  }
  if (any(nms == "")) {
    stop("All options passed through `...` must be named.", call. = FALSE)
  }

  # `horizon` has been documented since 0.0.7 but never read: the time-slice
  # loop is hard-coded to look one step ahead. Accepting any value would keep
  # promising something that does not happen.
  if (!is.null(opt$horizon) && !isTRUE(all.equal(opt$horizon, 1))) {
    stop("`horizon` is fixed at 1 in this version; got ", opt$horizon, ".",
         call. = FALSE)
  }

  unknown <- setdiff(nms, known_fit_options)
  if (length(unknown) == 0) {
    return(invisible(TRUE))
  }

  msgs <- vapply(unknown, function(u) {
    d <- utils::adist(u, known_fit_options, ignore.case = TRUE)[1, ]
    near <- known_fit_options[d == min(d)][1]
    # only suggest when it is plausibly a typo rather than a different word
    if (min(d) <= max(3, ceiling(nchar(u) / 3))) {
      paste0("unknown option \"", u, "\" - did you mean \"", near, "\"?")
    } else {
      paste0("unknown option \"", u, "\"")
    }
  }, character(1))

  stop(paste(msgs, collapse = "\n  "), call. = FALSE)
}

# Resolve the three-way `threshold` argument together with `threshold_type` and
# `threshold_const`, once, so the five call sites share one interpretation.
#
# Returns: active (is any thresholding done), value (user-supplied scalar tau,
# or NULL), type, const.
resolve_threshold <- function(opt) {
  threshold <- opt$threshold
  if (is.null(threshold)) {
    threshold <- FALSE
  }

  if (length(threshold) != 1) {
    stop("`threshold` must be length 1, got length ", length(threshold), ".",
         call. = FALSE)
  }
  if (is.na(threshold)) {
    stop("`threshold` must not be NA.", call. = FALSE)
  }

  type_given <- !is.null(opt$threshold_type)
  type <- if (type_given) opt$threshold_type else "theory"
  checkmate::assert_choice(type, c("theory", "soft", "hard"))

  const_given <- !is.null(opt$threshold_const)
  const <- if (const_given) opt$threshold_const else 1
  if (length(const) != 1 || !is.numeric(const) || !is.finite(const) ||
        const <= 0) {
    stop("`threshold_const` must be a single finite positive number.",
         call. = FALSE)
  }

  if (is.logical(threshold)) {
    if (!threshold) {
      return(list(active = FALSE, value = NULL, type = type, const = const))
    }
    if (const_given && type != "theory") {
      warning("`threshold_const` applies only to threshold_type = \"theory\"",
              "; ignoring it.", call. = FALSE)
    }
    # A silent change of cutoff is worse than a note. Removed in 1.2.0.
    if (!type_given) {
      notify_threshold_default()
    }
    return(list(active = TRUE, value = NULL, type = type, const = const))
  }

  if (is.numeric(threshold)) {
    if (!is.finite(threshold) || threshold <= 0) {
      # threshold = 0 is far more likely to be a mistake than an intent
      stop("A numeric `threshold` must be finite and > 0, got ", threshold,
           ".", call. = FALSE)
    }
    if (type_given) {
      warning("`threshold` was given a numeric cutoff, so `threshold_type` ",
              "is ignored.", call. = FALSE)
    }
    if (const_given) {
      warning("`threshold` was given a numeric cutoff, so `threshold_const` ",
              "is ignored.", call. = FALSE)
    }
    return(list(active = TRUE, value = threshold, type = type, const = const))
  }

  stop("`threshold` must be FALSE, TRUE, or a positive number, got ",
       class(threshold)[1], ".", call. = FALSE)
}

# One note per session, not per fit.
threshold_default_notified <- new.env(parent = emptyenv())
notify_threshold_default <- function() {
  if (isTRUE(threshold_default_notified$done)) {
    return(invisible(NULL))
  }
  threshold_default_notified$done <- TRUE
  message(
    "sparsevar: threshold_type now defaults to \"theory\" (was \"soft\").\n",
    "  Set threshold_type = \"soft\" for the previous behaviour. ",
    "This note is shown once per session."
  )
  invisible(NULL)
}

# Options error_bands_irf() and jackknife() read from `...`. Same purpose as
# known_fit_options: these paths do not forward `...` any further, so anything
# else is a typo.
known_irf_options <- c(
  "threshold", "threshold_type", "threshold_const", "verbose", "mode"
)

validate_irf_options <- function(opt) {
  bad <- setdiff(names(opt), known_irf_options)
  if (length(bad) > 0) {
    stop(paste0("unknown option \"", bad, "\"", collapse = "\n  "),
         call. = FALSE)
  }
  invisible(TRUE)
}

# The resampled series on the scale its coefficients were fitted on. The IRF
# bootstrap paths fit with var_scad()/var_mcp(), which standardize internally
# and return only coefficients, so the scale has to be recomputed here. Uses
# standardization() rather than transform_data() because only the series is
# wanted, not the Kronecker-expanded design.
boot_series <- function(sample_data, v, p) {
  opt <- list(method = v$method, penalty = v$penalty)
  apply_standardization(sample_data, standardization(sample_data, opt))
}

# fit_var() arguments that reproduce a resolved threshold spec, for the IRF
# paths that re-enter fit_var() rather than thresholding in place.
threshold_fit_args <- function(sample_data, p, v, thr) {
  args <- list(sample_data, p, penalty = v$penalty, method = v$method)
  if (!thr$active) {
    return(args)
  }
  if (!is.null(thr$value)) {
    args$threshold <- thr$value
    return(args)
  }
  args$threshold <- TRUE
  args$threshold_type <- thr$type
  args$threshold_const <- thr$const
  args
}

# The location and scale a fit standardizes by, as a pair of vectors.
#
# Split out of transform_data() so that time-slice validation can put its
# held-out rows on exactly the same footing as the coefficients it is scoring.
# Restating the rule at the call site is how the two drifted apart: the
# validation loop scored centered coefficients against raw levels, which made
# lambda selection depend on where the series happened to sit.
#
# Note the timeSlice branch takes the mean over the window minus leaveOut, not
# over the whole window, so no observation used for validation contributes to
# the centering constant.
standardization <- function(data, opt) {
  data <- as.matrix(data)
  nr <- nrow(data)
  nc <- ncol(data)

  # NB: named do_scale, not scale, so it cannot shadow stats::scale for
  # readers of the code below.
  do_scale <- ifelse(is.null(opt$scale), FALSE, opt$scale)
  # center the matrix columns (default)
  center <- ifelse(is.null(opt$center), TRUE, opt$center)

  if (center == TRUE) {
    if (opt$method == "timeSlice") {
      leave_out <- ifelse(is.null(opt$leaveOut), 10, opt$leaveOut)
      m <- colMeans(data[1:(nr - leave_out), ])
    } else {
      m <- colMeans(data)
    }
  } else {
    m <- rep(0, nc)
  }

  if (do_scale == TRUE) {
    sds <- apply(sweep(data, 2, m, "-"), 2, stats::sd)
    sds[sds == 0 | !is.finite(sds)] <- 1
  } else {
    sds <- rep(1, nc)
  }

  list(mu = m, sd = sds)
}

# Put a data matrix on the scale described by standardization().
apply_standardization <- function(data, std) {
  data <- as.matrix(data)
  sweep(sweep(data, 2, std$mu, "-"), 2, std$sd, "/")
}

#' @title Transform data
#'
#' @description Transform the input data
#'
#' @usage transform_data(data, p, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param opt a list containing the options
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X}: the design matrix for the VAR estimation
#'   \item \code{y}: the response vector
#'   \item \code{series}: the (possibly centered/scaled) data matrix
#'   \item \code{mu}: the column means used for centering
#'   \item \code{sd}: the column standard deviations used for scaling; all
#'     ones when \code{opt$scale} is \code{FALSE} (the default)
#' }
#' @export
transform_data <- function(data, p, opt) {

  # get the number of rows and columns
  nr <- nrow(data)
  nc <- ncol(data)

  # make sure the data is in matrix format
  data <- as.matrix(data)

  # The location and scale are recorded, not just applied: everything
  # downstream (coefficients, residuals, forecasts, IRF bands) is then on the
  # standardized scale, and without these there is no way back to the data's
  # own units.
  std <- standardization(data, opt)
  m <- std$mu
  sds <- std$sd

  data <- apply_standardization(data, std)

  # create Xs and Ys (temp variables)
  tmp_x <- data[1:(nr - 1), ]
  tmp_y <- data[2:(nr), ]

  # create the data matrix
  tmp_x <- duplicate_matrix(tmp_x, p)
  tmp_y <- tmp_y[p:nrow(tmp_y), ]

  y <- as.vector(tmp_y)

  # Hadamard product for data
  i <- Matrix::Diagonal(nc)
  x <- kronecker(i, tmp_x)

  output <- list()
  output$X <- x
  output$y <- y
  output$series <- data
  output$mu <- t(m)
  output$sd <- sds
  output
}

#' @title VAR ENET
#'
#' @description Estimate VAR using ENET penalty
#'
#' @usage var_enet(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#'
#' @return A \code{glmnet} object containing the fitted model.
#' @export
var_enet <- function(data, p, lambdas, opt) {
  # transform the dataset
  tr_dt <- transform_data(data, p, opt)
  fit <- glmnet::glmnet(tr_dt$X, tr_dt$y, lambda = lambdas)
  fit
}

#' @title VAR SCAD
#'
#' @description Estimate VAR using SCAD penalty
#'
#' @usage var_scad(data, p, lambdas, opt, penalty)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#' @param penalty a string "SCAD" or something else
#'
#' @return An \code{ncvreg} object containing the fitted model.
#' @export

var_scad <- function(data, p, lambdas, opt, penalty = "SCAD") {
  # transform the dataset
  tr_dt <- transform_data(data, p, opt)

  if (penalty == "SCAD") {
    fit <- ncvreg::ncvreg(as.matrix(tr_dt$X), tr_dt$y,
      family = "gaussian", penalty = "SCAD",
      alpha = 1, lambda = lambdas
    )
  } else {
    stop("[WIP] Only SCAD regression is supported at the moment")
  }
  fit
}

#' @title VAR MCP
#'
#' @description Estimate VAR using MCP penalty
#'
#' @usage var_mcp(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#'
#' @return An \code{ncvreg} object containing the fitted model.
#' @export
var_mcp <- function(data, p, lambdas, opt) {
  # transform the dataset
  tr_dt <- transform_data(data, p, opt)
  fit <- ncvreg::ncvreg(as.matrix(tr_dt$X), tr_dt$y,
    family = "gaussian", penalty = "MCP",
    alpha = 1, lambda = lambdas
  )
  fit
}

split_matrix <- function(m, p) {
  nr <- nrow(m)
  a <- list()
  for (i in 1:p) {
    ix <- ((i - 1) * nr) + (1:nr)
    a[[i]] <- m[1:nr, ix]
  }
  a
}

duplicate_matrix <- function(data, p) {
  nr <- nrow(data)
  nc <- ncol(data)

  output_data <- data

  if (p > 1) {
    for (i in 1:(p - 1)) {
      tmp_data <- matrix(0, nrow = nr, ncol = nc)
      tmp_data[(i + 1):nr, ] <- data[1:(nr - i), ]
      output_data <- cbind(output_data, tmp_data)
    }
  }

  output_data <- output_data[p:nr, ]
  output_data
}

compute_residuals <- function(data, a) {
  nr <- nrow(data)
  nc <- ncol(data)
  p <- length(a)

  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)

  for (i in 1:p) {
    tmp_d <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    tmp_f <- t(a[[i]] %*% t(tmp_d))
    f <- f + tmp_f
  }

  res <- data - f
  res
}

#' @title Companion VAR
#'
#' @description Build the VAR(1) representation of a VAR(p) process
#'
#' @usage companion_var(v)
#'
#' @param v the VAR object as from \code{fitVAR} or \code{simulateVAR}
#'
#' @return A sparse matrix (of class \code{dgCMatrix}) representing the
#'   companion form of the VAR(p) process.
#' @export
companion_var <- function(v) {
  if (!check_is_var(v)) {
    stop("Expected var object, got ", class(v))
  }
  a <- v$A
  nc <- ncol(a[[1]])
  p <- length(a)
  if (p > 1) {
    big_a <- Matrix::Matrix(0, nrow = p * nc, ncol = p * nc, sparse = TRUE)
    for (k in 1:p) {
      ix <- ((k - 1) * nc) + (1:nc)
      big_a[1:nc, ix] <- a[[k]]
    }
    ix_r <- (nc + 1):nrow(big_a)
    ix_c <- 1:((p - 1) * nc)
    big_a[ix_r, ix_c] <- diag(1, nrow = length(ix_c), ncol = length(ix_c))
  } else {
    big_a <- Matrix::Matrix(a[[1]], sparse = TRUE)
  }
  big_a
}

#' @title Bootstrap VAR
#'
#' @description Build the bootstrapped series from the original var
#'
#' @usage bootstrapped_var(v)
#'
#' @param v the VAR object as from fitVAR or simulateVAR
#'
#' @return A matrix containing the bootstrapped time series with the same
#'   dimensions as the original series.
#' @export
bootstrapped_var <- function(v) {

  ## This function creates the bootstrapped time series
  if (!check_is_var(v)) {
    stop("Expected var object, got ", class(v))
  }

  r <- v$residuals
  s <- v$series
  a <- v$A
  n <- ncol(a[[1]])
  p <- length(a)
  t <- nrow(r)
  r <- r - matrix(colMeans(r), ncol = n, nrow = t)

  zt <- matrix(0, nrow = t, ncol = n)
  zt[1:p, ] <- s[1:p, ]

  for (t0 in (p + 1):t) {
    ix <- sample((p + 1):t, 1)
    u <- r[ix, ]
    vv <- rep(0, n)
    for (i in 1:p) {
      ph <- a[[i]]
      vv <- vv + ph %*% zt[(t0 - i), ]
    }
    vv <- vv + u
    zt[t0, ] <- vv
  }
  zt
}

#' @title Test for Granger Causality
#'
#' @description Tests, for every ordered pair of variables \eqn{(i, j)}, the
#' null hypothesis that variable \eqn{j} does not Granger-cause variable
#' \eqn{i}, that is
#' \deqn{H_0: A_1[i, j] = A_2[i, j] = \ldots = A_p[i, j] = 0.}
#' This is a joint hypothesis across all \eqn{p} lags, and it is assessed with
#' a sup-norm bootstrap simultaneous confidence region built from the
#' resampled coefficient matrices stored in \code{eb$Aboot}. For each pair the
#' studentised statistic \eqn{\max_l |A_l[i,j]| / s_l} is compared against the
#' \eqn{1 - \alpha} quantile of its bootstrap counterpart, where \eqn{s_l} is
#' the bootstrap standard deviation of \eqn{A_l[i,j]}. Lags whose bootstrap
#' standard deviation is zero carry no information and are excluded from the
#' maximum.
#'
#' Because the null is joint over lags, a pair that fails to reject has *all*
#' of its lags set to zero, which is the interpretation Granger causality
#' requires.
#'
#' Note that the bootstrap used here resamples residuals around a penalised
#' fit. Such a bootstrap is not in general consistent for LASSO-type
#' estimators, so the resulting levels should be read as approximate.
#'
#' @usage test_granger(v, eb, alpha = 0.05)
#'
#' @param v the VAR object as from fitVAR or simulateVAR
#' @param eb the error bands as obtained from \code{error_bands_irf}; must
#' contain the \code{Aboot} array of resampled coefficient matrices
#' @param alpha significance level of the joint test (default 0.05)
#'
#' @return A list of \code{p} matrices containing only the coefficients of
#'   variable pairs for which the non-causality null was rejected; all lags of
#'   a non-causal pair are set to zero. The logical \code{n x n} causality
#'   matrix is attached as the attribute \code{"causality"}, with entry
#'   \eqn{(i, j)} \code{TRUE} when \eqn{j} Granger-causes \eqn{i}.
#' @export
test_granger <- function(v, eb, alpha = 0.05) {
  # The studentised sup-norm statistic is only calibrated when the resampling
  # scheme reproduces sampling variability. Jackknife leave-one-out
  # perturbations are far smaller than that, so the s_l below are too small and
  # the statistic is inflated: the test then over-rejects.
  resampling <- attr(eb, "resampling")
  if (!is.null(resampling) && resampling != "bootstrap") {
    stop("test_granger() requires bootstrap error bands, but `eb` was ",
         "computed with resampling = \"", resampling, "\". ",
         "Leave-one-out perturbations are much smaller than sampling ",
         "variability, so the resampled standard deviations are too small, ",
         "the studentised statistic is inflated, and the test over-rejects. ",
         "Recompute with ",
         "error_bands_irf(..., resampling = \"bootstrap\").")
  }

  if (is.null(eb$Aboot)) {
    stop("eb does not contain the resampled coefficient matrices (Aboot). ",
         "Recompute the error bands with the current version of ",
         "error_bands_irf().")
  }

  p <- length(v$A)
  nc <- ncol(v$A[[1]])
  a_boot <- eb$Aboot
  nboot <- dim(a_boot)[3]

  # point estimates stacked as [A_1 | ... | A_p], matching Aboot
  a_hat <- do.call(cbind, v$A)

  # bootstrap standard deviation of every coefficient
  s <- apply(a_boot, c(1, 2), stats::sd)
  usable <- s > 0

  causal <- matrix(FALSE, nrow = nc, ncol = nc)

  for (i in 1:nc) {
    for (j in 1:nc) {
      # column indices of entry (i, j) across the p lags
      ix <- j + nc * (0:(p - 1))
      keep <- ix[usable[i, ix]]

      if (length(keep) == 0) {
        # no bootstrap variability at any lag: nothing to test against
        next
      }

      s_ij <- s[i, keep]

      # observed studentised sup-norm statistic
      t_obs <- max(abs(a_hat[i, keep]) / s_ij)

      # bootstrap counterpart, centred at the point estimate. draws is
      # (number of usable lags) x (number of resamples), so both the centring
      # and the scaling recycle down the columns.
      draws <- matrix(a_boot[i, keep, ], nrow = length(keep), ncol = nboot)
      t_star <- apply(abs(draws - a_hat[i, keep]) / s_ij, 2, max)

      causal[i, j] <- t_obs > stats::quantile(t_star, 1 - alpha, names = FALSE)
    }
  }

  out <- list()
  for (l in 1:p) {
    out[[l]] <- v$A[[l]] * causal
  }
  attr(out, "causality") <- causal
  out
}

#' @title Computes information criteria for VARs
#'
#' @description This function computes information criteria (AIC, Schwartz and
#' Hannan-Quinn) for VARs.
#'
#' @usage inform_crit(v)
#'
#' @param v a list of VAR objects as from fitVAR.
#'
#' @return A data frame with columns \code{AIC}, \code{BIC}, and
#'   \code{HannanQuinn} containing the information criteria values for each
#'   VAR model in the input list.
#' @export
inform_crit <- function(v) {
  if (is.list(v)) {
    k <- length(v)
    r <- matrix(0, nrow = k, ncol = 3)
    for (i in 1:k) {
      if (inherits(v[[i]], c("var", "vecm"))) {
        p <- length(v[[i]]$A)
        # Compute sparsity
        s <- 0
        for (l in 1:p) {
          s <- s + sum(v[[i]]$A[[l]] != 0)
        }
        sp <- s / (p * ncol(v[[i]]$A[[1]])^2)
      } else {
        stop("List elements must be var or vecm objects.")
      }
      sigma <- v[[i]]$sigma
      nr <- nrow(v[[i]]$residuals)
      nc <- ncol(v[[i]]$residuals)
      d <- det(sigma)

      r[i, 1] <- log(d) + (2 * p * sp * nc^2) / nr # AIC
      r[i, 2] <- log(d) + (log(nr) / nr) * (p * sp * nc^2) # BIC
      r[i, 3] <- log(d) + (2 * p * sp * nc^2) / nr * log(log(nr)) # Hannan-Quinn
    }
    results <- data.frame(r)
    colnames(results) <- c("AIC", "BIC", "HannanQuinn")
  } else {
    stop("Input must be a list of var models.")
  }
  results
}

estimate_covariance <- function(res, ...) {
  nc <- ncol(res)
  s <- corpcor::cov.shrink(res, verbose = FALSE)
  sigma <- matrix(0, ncol = nc, nrow = nc)
  for (i in 1:nc) {
    for (j in 1:nc) {
      sigma[i, j] <- s[i, j]
    }
  }
  sigma
}

#' @title Computes forecasts for VARs
#'
#' @description This function computes forecasts for a given VAR.
#'
#' @usage compute_forecasts(v, num_steps)
#'
#' @param v a VAR object as from fitVAR.
#' @param num_steps the number of forecasts to produce.
#'
#' @return A matrix of dimension (number of variables) x (num_steps) containing
#'   the forecasted values for each variable at each forecast horizon.
#' @export
compute_forecasts <- function(v, num_steps = 1) {
  if (!check_is_var(v)) {
    stop("You must pass a var object.")
  } else {
    mu <- v$mu
    # all-ones for fits made without scale = TRUE, and for objects saved by
    # versions that predate the scale being recorded
    sds <- if (is.null(v$sd)) rep(1, length(mu)) else v$sd
    data <- v$series
    v <- v$A
  }

  if (!is.list(v)) {
    stop("v must be a var object or a list of matrices.")
  } else {
    nr <- nrow(data)
    nc <- ncol(v[[1]])
    p <- length(v)

    f <- matrix(0, nrow = nc, ncol = num_steps)

    tmp_data <- matrix(data = t(data[(nr - p + 1):nr, ]),
                       nrow = nc,
                       ncol = num_steps)
    nr <- ncol(tmp_data)

    for (n in 1:num_steps) {
      for (k in 1:p) {
        if (n == 1) {
          f[, n] <- f[, n] + v[[k]] %*% tmp_data[, nr - k + 1]
        } else {
          if (nr > 1) {
            tmp_data <- cbind(tmp_data[, 2:nr], f[, n - 1])
          } else {
            tmp_data <- as.matrix(f[, n - 1])
          }
          f[, n] <- f[, n] + v[[k]] %*% tmp_data[, nr - k + 1]
        }
      }
    }
  }
  # `f` is on the scale the model was fitted on. Undo the standardization
  # before adding back the means, which were taken on the raw data: mixing the
  # two returns a standardized deviation on top of a raw-scale mean.
  f <- f * sds
  f <- f + matrix(rep(mu, num_steps), nrow = length(mu), ncol = num_steps)
  f
}

# Cutoffs. `nr` is the number of design rows, T - p, not nrow(data): the
# threshold is a property of the regression, not of the raw series.
#
# All three types select a *cutoff*; the *rule* applied is always hard
# thresholding, a * 1(|a| >= tau). "soft" and "hard" are historical names for
# two cutoffs, not two rules — see ?fit_var.
apply_threshold <- function(a_mat, nr, nc, p,
                            type = "theory",
                            sigma = NULL,
                            s = NULL,
                            const = 1) {
  if (type == "soft") {
    tr <- 1 / (p * nc * log(nr))
    # another threshold: tr <- 1 / sqrt(p * nc * log(nr))
  } else if (type == "hard") {
    tr <- (nc)^(-0.49)
  } else if (type == "theory") {
    # A missing scale silently degrades this to the bare scalar cutoff, which
    # is measurably the worse estimator, so it is an error rather than a
    # default.
    if (is.null(sigma) || is.null(s)) {
      stop("threshold_type = \"theory\" needs both sigma and s.")
    }
    if (length(sigma) != nc || length(s) != nc) {
      stop("sigma and s must both have length ", nc, ".")
    }
    # tau_ij = (sigma_i / s_j) * C * sqrt(log(n^2 p) / (T - p)), an
    # nc x (nc * p) matrix aligned with a_mat. rep(s, p) tiles one sd per
    # series across all p lag blocks.
    tr <- outer(sigma, rep(s, p), "/") * const * sqrt(log(nc^2 * p) / nr)
  } else {
    stop(
      "Unknown threshold type. Possible values are: ",
      "\"theory\", \"soft\" or \"hard\""
    )
  }

  # Matrix tr broadcasts elementwise; a scalar tr still works.
  l_mat <- abs(a_mat) >= tr
  a_mat * l_mat
}

# The whole thresholding step for one fit path, so the five call sites do not
# each re-derive it. `thr` is the resolved spec from resolve_threshold();
# `series` is the transformed series the coefficients were fitted on.
threshold_coefficients <- function(a_mat, series, nc, p, thr) {
  if (!thr$active) {
    return(a_mat)
  }

  # design rows, T - p (01). Note this is *not* nrow(tr_dt$X), which is
  # nc * (T - p) because X is the Kronecker-expanded stacked design.
  nr <- nrow(series) - p

  if (!is.null(thr$value)) {
    return(a_mat * (abs(a_mat) >= thr$value))
  }

  if (thr$type != "theory") {
    # the legacy cutoffs are frozen and scale-free
    return(apply_threshold(a_mat, nr, nc, p, type = thr$type))
  }

  # sigma_i comes from the *un-thresholded* fit, so it needs its own residual
  # pass. Only "theory" pays for it (~0.08% of a fit).
  res0 <- compute_residuals(series, split_matrix(a_mat, p))
  # margin 2: residuals are T x n. compute_residuals() zero-pads the first p
  # rows, so drop them — which also agrees with the T - p denominator.
  sigma <- apply(res0[(p + 1):nrow(res0), , drop = FALSE], 2, stats::sd)
  s <- apply(series, 2, stats::sd)

  apply_threshold(a_mat, nr, nc, p, type = thr$type,
                  sigma = sigma, s = s, const = thr$const)
}
