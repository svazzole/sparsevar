#' @title VAR simulation
#'
#' @description This function generates a simulated multivariate VAR time
#' series.
#'
#' @usage simulate_var(n, p, nobs, rho, sparsity, mu, method, covariance, ...)
#'
#' @param n dimension of the time series (default \code{n = 100}).
#' @param p number of lags of the VAR model (default \code{p = 1}).
#' @param nobs number of observations to be generated (default
#' \code{nobs = 250}).
#' @param rho base value for the covariance matrix (default \code{rho = 0.5}).
#' @param sparsity density (in percentage) of the number of nonzero elements
#' of the VAR matrices (default \code{sparsity = 0.05}).
#' @param mu a vector containing the mean of the simulated process (default
#' \code{mu = 0}).
#' @param method which method to use to generate the VAR matrix. Possible values
#' are \code{"normal"} or \code{"bimodal"} ((default \code{method = "normal"})).
#' @param covariance type of covariance matrix to use in the simulation.
#' Possible values: \code{"Toeplitz"}, \code{"block1"}, \code{"block2"},
#' \code{"Wishart"} or simply \code{"diagonal"}
#' (default \code{covariance = "Toeplitz"}).
#' @param ... the options for the simulation. These are:
#' \code{muMat}: the mean of the entries of the VAR matrices;
#' \code{sdMat}: the sd of the entries of the matrices;
#'
#' @return A a list of NxN matrices ordered by lag
#' @return data a list with two elements: \code{series} the multivariate
#' time series and \code{noises} the time series of errors
#' @return S the variance/covariance matrix of the process
#'
#' @export
simulate_var <- function(n = 100, p = 1, nobs = 250,
                         rho = 0.5, sparsity = 0.05,
                         mu = 0, method = "normal",
                         covariance = "Toeplitz", ...) {
  opt <- list(...)

  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_count(p, positive = TRUE)
  checkmate::assert_count(nobs, positive = TRUE)

  checkmate::assert_number(sparsity, lower = 0)

  checkmate::assert_choice(method, c("normal", "bimodal"))
  checkmate::assert_choice(covariance,
                           c("Toeplitz", "Wishart",
                             "block1", "block2", "diagonal"))

  # Create a var object to save the matrices (the output)
  out <- list()
  attr(out, "class") <- "var"
  attr(out, "type") <- "simulation"
  out$A <- list()

  fixed_mat <- opt$fixedMat
  if (!is.null(fixed_mat)) {
    # The user passed a list of matrices
    out$A <- fixed_mat
    if (!check_matrices(out$A)) {
      stop("The matrices you passed are incompatible.")
    }
    c_var <- as.matrix(companion_var(out))
    if (max(Mod(eigen(c_var)$values)) >= 1) {
      warning("The VAR you passed is unstable.")
    }
  } else {
    stable <- FALSE
    while (!stable) {
      for (i in 1:p) {
        out$A[[i]] <- create_sparse_matrix(sparsity = sparsity, n = n,
                                           method = method, stationary = TRUE,
                                           p = p, ...)
        l <- max(Mod(eigen(out$A[[i]])$values))
        while ((l > 1) || (l == 0)) {
          out$A[[i]] <- create_sparse_matrix(sparsity = sparsity, n = n,
                                             method = method, stationary = TRUE,
                                             p = p, ...)
          l <- max(Mod(eigen(out$A[[i]])$values))
        }
      }
      c_var <- as.matrix(companion_var(out))
      if (max(Mod(eigen(c_var)$values)) < 1) {
        stable <- TRUE
      }
    }
  }

  # Covariance Matrix: Toeplitz, Block1 or Block2
  if (covariance == "block1") {
    l <- floor(n / 2)
    i <- diag(1 - rho, nrow = n)
    r <- matrix(0, nrow = n, ncol = n)
    r[1:l, 1:l] <- rho
    r[(l + 1):n, (l + 1):n] <- diag(rho, nrow = (n - l))
    c <- i + r
  } else if (covariance == "block2") {
    l <- floor(n / 2)
    i <- diag(1 - rho, nrow = n)
    r <- matrix(0, nrow = n, ncol = n)
    r[1:l, 1:l] <- rho
    r[(l + 1):n, (l + 1):n] <- rho
    c <- i + r
  } else if (covariance == "Toeplitz") {
    r <- rho^(1:n)
    c <- Matrix::toeplitz(r)
  } else if (covariance == "Wishart") {
    r <- rho^(1:n)
    s <- Matrix::toeplitz(r)
    c <- stats::rWishart(1, 2 * n, s)
    c <- as.matrix(c[, , 1])
  } else {
    c <- diag(x = rho, nrow = n, ncol = n)
  }

  # Adjust Signal to Noise Ratio
  snr <- opt$SNR
  if (!is.null(snr)) {
    if (snr <= 0) {
      stop("Signal to Noise Ratio must be greater than 0.")
    }
    s <- max(abs(c_var)) / snr
    c <- diag(s, n, n) %*% c %*% diag(s, n, n)
  }

  # Generate the VAR process
  data <- generate_var_series(nobs = nobs, mu = mu, ar = out$A,
                              sigma = c, skip = 200)

  # Complete the output
  out$series <- data$series
  out$noises <- data$noises
  out$sigma <- c

  out
}

generate_var_series <- function(nobs, mu, ar, sigma, skip = 200) {

  # This function creates the simulated time series

  n <- nrow(sigma)
  n_t <- nobs + skip
  at <- mvtnorm::rmvnorm(n_t, rep(0, n), sigma)

  p <- length(ar)

  ist <- p + 1
  zt <- matrix(0, n_t, n)

  if (length(mu) == 0) {
    mu <- rep(0, n)
  }

  for (it in ist:n_t) {
    tmp <- matrix(at[it, ], 1, n)

    for (i in 1:p) {
      ph <- ar[[i]]
      ztm <- matrix(zt[it - i, ], 1, n)
      tmp <- tmp + ztm %*% t(ph)
    }

    zt[it, ] <- mu + tmp
  }

  # skip the first skip points to initialize the series
  zt <- zt[(1 + skip):n_t, ]
  at <- at[(1 + skip):n_t, ]

  out <- list()
  out$series <- zt
  out$noises <- at
  out
}

check_matrices <- function(a) {

  checkmate::assert_list(a)

  l <- length(a)
  if (l > 1) {
    for (i in 1:(l - 1)) {
      if (sum(1 - (dim(a[[i]]) == dim(a[[i + 1]]))) != 0) {
        return(FALSE)
      }
    }
  }
  TRUE
}
