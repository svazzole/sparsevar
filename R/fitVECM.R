#' @title Multivariate VECM estimation
#'
#' @description A function to estimate a (possibly big)
#' multivariate VECM time series using penalized least squares
#' methods, such as ENET, SCAD or MC+.
#'
#' @usage fit_vecm(data, p, penalty, method, log_scale, ...)
#'
#' @param data the data from the time series: variables
#' in columns and observations in rows
#' @param p order of the VECM model
#' @param penalty the penalty function to use.
#' Possible values are \code{"ENET"}, \code{"SCAD"} or \code{"MCP"}
#' @param log_scale should the function consider the
#' \code{log} of the inputs? By default this is set to \code{TRUE}
#' @param method \code{"cv"} or \code{"timeSlice"}
#' @param ... options for the function (TODO: specify)
#'
#' @return Pi the matrix \code{Pi} for the VECM model
#' @return G the list (of length \code{p-1}) of the
#' estimated matrices of the process
#' @return fit the results of the penalized LS estimation
#' @return mse the mean square error of the cross validation
#' @return time elapsed time for the estimation
#'
#' @export
fit_vecm <- function(
  data,
  p = 0,
  penalty = "ENET",
  method = "cv",
  log_scale = TRUE,
  ...
) {
  nc <- ncol(data)

  p <- p + 1

  opt <- list(...)
  # NB: this used to set opt$center <- FALSE, but opt is never passed to
  # fit_var() below (only `...` is), so the override never took effect and the
  # data has always been centered. Honouring it makes the estimate worse:
  # on the cointegrated example the recovered beta goes from (1, -2.083) to
  # (1, -2.160) against a true (1, -2). Centering is the correct behaviour, so
  # the dead override is removed rather than wired up.

  # by default log-scale the data
  if (log_scale == TRUE) {
    data <- log(data)
    data[is.na(data)] <- 0
  }

  results_var <- fit_var(data, p = p, penalty = penalty, method = method, ...)
  m <- results_var$A
  i <- diag(x = 1, nrow = nc, ncol = nc)

  # Coint matrix
  pi_mat <- -(i - matrix_sum(m, ix = 1))

  # Gamma matrices
  g <- list()

  if (p > 1) {
    for (k in 1:(p - 1)) {
      g[[k]] <- -matrix_sum(m, ix = k + 1)
    }
  }

  output <- list()
  output$mu <- results_var$mu
  output$Pi <- pi_mat
  output$G <- g
  output$A <- results_var$A
  output$fit <- results_var$fit
  output$mse <- results_var$mse
  output$mseSD <- results_var$mseSD
  output$time <- results_var$time
  output$residuals <- results_var$residuals
  output$lambda <- results_var$lambda
  output$series <- results_var$series

  if (is.null(opt$methodCov)) {
    output$sigma <- estimate_covariance(output$residuals)
  } else {
    output$sigma <- estimate_covariance(output$residuals,
                                        methodCovariance = opt$methodCov)
  }

  output$penalty <- results_var$penalty
  output$method <- results_var$method
  attr(output, "class") <- "vecm"
  attr(output, "type") <- "fit"

  output
}

matrix_sum <- function(m, ix = 1) {
  l <- length(m)
  nc <- ncol(m[[1]])
  a <- matrix(0, nrow = nc, ncol = nc)

  for (i in ix:l) {
    a <- a + m[[i]]
  }

  a
}

#' @title Decompose Pi VECM matrix
#'
#' @description A function to estimate a (possibly big) multivariate
#' VECM time series using penalized least squares methods, such as
#' ENET, SCAD or MC+.
#'
#' @usage decompose_pi(vecm, rk, ...)
#'
#' @param vecm the VECM object
#' @param rk rank
#' @param ... options for the function (TODO: specify)
#'
#' @return alpha
#' @return beta
#'
#' @export
decompose_pi <- function(vecm, rk, ...) {
  if (!inherits(vecm, "vecm")) {
    stop("The input is not a vecm object.")
  }

  # Different covariance methods?
  # probably in opt <- list(...)

  nc <- ncol(vecm$Pi)
  pi_mat <- vecm$Pi
  colnames(pi_mat) <- NULL
  rownames(pi_mat) <- NULL
  sig <- corpcor::invcov.shrink(vecm$residuals, verbose = FALSE)
  colnames(sig) <- NULL
  rownames(sig) <- NULL

  if (rk < nc && rk > 0) {
    a <- pi_mat[, 1:rk]
    b <- t(solve(t(a) %*% sig %*% a)
           %*% (t(a) %*% sig %*% pi_mat[, (rk + 1):nc]))
    b <- rbind(diag(1, rk, rk), b)
  } else if (rk == nc) {
    a <- pi_mat
    b <- diag(1, rk, rk)
  } else {
    a <- numeric(length = nc)
    b <- pi_mat
  }

  out <- list()
  out$alpha <- a
  out$beta <- b
  out
}

decompose_pi2 <- function(vecm, rk) {
  if (!inherits(vecm, "vecm")) {
    stop("The input is not a vecm object.")
  }

  nc <- ncol(vecm$Pi)
  pi_mat <- vecm$Pi
  colnames(pi_mat) <- NULL
  rownames(pi_mat) <- NULL

  if (rk >= 1) {
    a <- pi_mat[, 1:rk]
    a <- kronecker(diag(1, nc, nc), a)
    b <- as.numeric(pi_mat)
    b <- matrix(qr.solve(a, b), ncol = rk, nrow = nc, byrow = TRUE)
    b_t <- matrix(qr.solve(a, b), ncol = rk, nrow = nc, byrow = FALSE)
  } else {
    a <- numeric(length = nc)
    b <- pi_mat
    b_t <- t(pi_mat)
  }

  out <- list()
  out$alpha <- a
  out$beta <- b
  out$betaT <- b_t
  out
}
