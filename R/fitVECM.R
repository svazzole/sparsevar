#' @title Multivariate VECM estimation
#'
#' @description A function to estimate a (possibly big)
#' multivariate VECM time series using penalized least squares
#' methods, such as ENET, SCAD or MC+.
#'
#' @usage fitVECM(data, p, penalty, method, logScale, ...)
#'
#' @param data the data from the time series: variables
#' in columns and observations in rows
#' @param p order of the VECM model
#' @param penalty the penalty function to use.
#' Possible values are \code{"ENET"}, \code{"SCAD"} or \code{"MCP"}
#' @param logScale should the function consider the
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
  opt$center <- FALSE

  # by default log-scale the data
  if (log_scale == TRUE) {
    data <- log(data)
    data[is.na(data)] <- 0
  }

  results_var <- fit_var(data, p = p, penalty = penalty, method = method, ...)
  m <- results_var$A
  i <- diag(x = 1, nrow = nc, ncol = nc)

  # Coint matrix
  pi_mat <- -(i - matrixSum(m, ix = 1))

  # Gamma matrices
  g <- list()

  if (p > 1) {
    for (k in 1:(p - 1)) {
      g[[k]] <- -matrixSum(m, ix = k + 1)
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

matrixSum <- function(M, ix = 1) {
  l <- length(M)
  nc <- ncol(M[[1]])

  A <- matrix(0, nrow = nc, ncol = nc)

  for (i in ix:l) {
    A <- A + M[[i]]
  }

  return(A)
}

#' @title Decompose Pi VECM matrix
#'
#' @description A function to estimate a (possibly big) multivariate VECM time series
#' using penalized least squares methods, such as ENET, SCAD or MC+.
#'
#' @usage decomposePi(vecm, rk, ...)
#'
#' @param vecm the VECM object
#' @param rk rank
#' @param ... options for the function (TODO: specify)
#'
#' @return alpha
#' @return beta
#'
#' @export
decomposePi <- function(vecm, rk, ...) {
  if (attr(vecm, "class") != "vecm") {
    stop("The input is not a vecm object.")
  }

  # Different covariance methods?
  opt <- list(...)

  nc <- ncol(vecm$Pi)
  Pi <- vecm$Pi
  colnames(Pi) <- NULL
  rownames(Pi) <- NULL
  sig <- corpcor::invcov.shrink(vecm$residuals, verbose = FALSE)
  colnames(sig) <- NULL
  rownames(sig) <- NULL

  if (rk < nc & rk > 0) {
    a <- Pi[, 1:rk]
    b <- t(solve(t(a) %*% sig %*% a) %*% (t(a) %*% sig %*% Pi[, (rk + 1):nc]))
    b <- rbind(diag(1, rk, rk), b)
  } else if (rk == nc) {
    a <- Pi
    b <- diag(1, rk, rk)
  } else {
    a <- numeric(length = nc)
    b <- Pi
  }

  out <- list()
  out$alpha <- a
  out$beta <- b
  return(out)
}

decomposePi2 <- function(vecm, rk) {
  if (attr(vecm, "class") != "vecm") {
    stop("The input is not a vecm object.")
  }

  nc <- ncol(vecm$Pi)
  Pi <- vecm$Pi
  colnames(Pi) <- NULL
  rownames(Pi) <- NULL

  if (rk >= 1) {
    a <- Pi[, 1:rk]
    # s <- solve(vecm$sigma)
    # b <- t(solve(t(a)%*%s%*%a)%*%(t(a)%*%s%*%vecm$Pi[,(rk+1):nc]))
    # b <- rbind(diag(1, rk, rk), b)
    A <- kronecker(diag(1, nc, nc), a)
    B <- as.numeric(Pi)
    b <- matrix(qr.solve(A, B), ncol = rk, nrow = nc, byrow = TRUE)
    bT <- matrix(qr.solve(A, B), ncol = rk, nrow = nc, byrow = FALSE)
  } else {
    a <- numeric(length = nc)
    b <- Pi
    bT <- t(Pi)
  }

  out <- list()
  out$alpha <- a
  out$beta <- b
  out$betaT <- bT
  return(out)
}
