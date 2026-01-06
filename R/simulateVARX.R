#' @title VARX simulation
#'
#' @description This function generates a simulated multivariate VAR
#' time series.
#'
#' @usage simulateVARX(N, K, p, m, nobs, rho,
#'                     sparsityA1, sparsityA2, sparsityA3,
#'                     mu, method, covariance, ...)
#'
#' @param N dimension of the time series.
#' @param K TODO
#' @param p number of lags of the VAR model.
#' @param m TODO
#' @param nobs number of observations to be generated.
#' @param rho base value for the covariance matrix.
#' @param sparsityA1 density (in percentage) of the number of nonzero elements
#' of the A1 block.
#' @param sparsityA2 density (in percentage) of the number of nonzero elements
#' of the A2 block.
#' @param sparsityA3 density (in percentage) of the number of nonzero elements
#' of the A3 block.
#' @param mu a vector containing the mean of the simulated process.
#' @param method which method to use to generate the VAR matrix. Possible values
#' are \code{"normal"} or \code{"bimodal"}.
#' @param covariance type of covariance matrix to use in the simulation.
#' Possible values: \code{"toeplitz"}, \code{"block1"}, \code{"block2"} or
#' simply \code{"diagonal"}.
#' @param ... the options for the simulation. These are:
#' \code{muMat}: the mean of the entries of the VAR matrices;
#' \code{sdMat}: the sd of the entries of the matrices;
#'
#' @return A a list of NxN matrices ordered by lag
#' @return data a list with two elements: \code{series} the multivariate time
#' series and \code{noises} the time series of errors
#' @return S the variance/covariance matrix of the process
#'
#' @export
simulate_varx <- function(n = 40, k = 10, p = 1, m = 1, nobs = 250, rho = 0.5,
                          sparsity_a1 = 0.05, sparsity_a2 = 0.5,
                          sparsity_a3 = 0.5, mu = 0, method = "normal",
                          covariance = "Toeplitz", ...) {
  opt <- list(...)
  # Opt has fixedMat as element: fixedMat <- opt$fixedMat
  snr <- opt$SNR

  # Create a var object to save the matrices (the output)
  out <- list()
  attr(out, "class") <- "varx"
  attr(out, "type") <- "simulation"

  out$A <- list()
  out$A1 <- list()
  out$A2 <- list()
  out$A3 <- list()
  out$A4 <- list()
  p_x <- max(p, m)

  # Create D matrices (null)
  for (i in 1:p_x) {
    out$A4[[i]] <- matrix(0, nrow = k, ncol = n)
  }

  stable <- FALSE

  while (!stable) {
    # Randomly select an order for C matrices in 1:pX
    s <- sample(1:p_x, 1)

    # Create random C matrices with a given sparsity
    for (i in 1:s) {
      out$A3[[i]] <- create_sparse_matrix(sparsity = sparsity_a3, N = k,
                                          method = method, stationary = TRUE,
                                          p = 1, ...)
      l <- max(Mod(eigen(out$A3[[i]])$values))
      while ((l > 1) || (l == 0)) {
        out$A3[[i]] <- create_sparse_matrix(sparsity = sparsity_a3, N = k,
                                            method = method, stationary = TRUE,
                                            p = 1, ...)
        l <- max(Mod(eigen(out$A3[[i]])$values))
      }
    }
    if (s < p_x) {
      for (i in (s + 1):p_x) {
        out$A3[[i]] <- matrix(0, nrow = k, ncol = k)
      }
    }

    # Create random A matrices with a given sparsity
    for (i in 1:p) {
      out$A1[[i]] <- create_sparse_matrix(sparsity = sparsity_a1, N = n,
                                          method = method, stationary = TRUE,
                                          p = p, ...)
      l <- max(Mod(eigen(out$A1[[i]])$values))
      while ((l > 1) || (l == 0)) {
        out$A1[[i]] <- create_sparse_matrix(sparsity = sparsity_a1, N = n,
                                            method = method, stationary = TRUE,
                                            p = p, ...)
        l <- max(Mod(eigen(out$A1[[i]])$values))
      }
    }
    if (p < p_x) {
      for (i in (p + 1):p_x) {
        out$A1[[i]] <- matrix(0, nrow = n, ncol = n)
      }
    }

    # Create random B matrices
    for (i in 1:m) {
      r <- max(k, n)
      tmp <- create_sparse_matrix(sparsity = sparsity_a2, N = r,
                                  method = method, stationary = TRUE,
                                  p = p, ...)
      out$A2[[i]] <- tmp[1:n, 1:k]
    }
    if (m < p_x) {
      for (i in (m + 1):p_x) {
        out$A2[[i]] <- matrix(0, nrow = n, ncol = k)
      }
    }

    # Now "glue" all the matrices together
    for (i in 1:p_x) {
      tmp1 <- cbind(out$A1[[i]], out$A2[[i]])
      tmp2 <- cbind(out$A4[[i]], out$A3[[i]])
      out$A[[i]] <- rbind(tmp1, tmp2)
    }

    c_var <- as.matrix(companion_var(out))
    if (max(Mod(eigen(c_var)$values)) < 1) {
      stable <- TRUE
    }
  }

  n <- n + k
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
  } else if (covariance == "diagonal") {
    c <- diag(x = rho, nrow = n, ncol = n)
  } else {
    stop("Unknown covariance matrix type. Possible choices are:
          toeplitz, block1, block2 or diagonal")
  }

  # Adjust Signal to Noise Ratio
  if (!is.null(snr)) {
    if (snr == 0) {
      stop("Signal to Noise Ratio must be greater than 0.")
    }
    s <- max(abs(c_var)) / opt$SNR
    c <- diag(s, n, n) %*% c %*% diag(s, n, n)
  }

  # Matrix for MA part: theta <- matrix(0, n, n)

  # Generate the VAR process
  data <- generate_var_series(nobs = nobs, mu, ar = out$A,
                              sigma = C, skip = 200)

  # Complete the output
  out$series <- data$series[, 1:(n - k)]
  out$Xt <- data$series[, (n - k + 1):n]
  out$noises <- data$noises
  out$sigma <- C

  out
}

generate_varx_series <- function(nobs, mu, ar, sigma, skip = 200) {

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

check_matrices_x <- function(a) {

  # This function check if all the matrices passed have the same dimensions
  if (!is.list(a)) {
    stop("The matrices must be passed in a list")
  } else {
    l <- length(a)
    if (l > 1) {
      for (i in 1:(l - 1)) {
        if (sum(1 - (dim(a[[i]]) == dim(a[[i + 1]]))) != 0) {
          FALSE
        }
      }
    }
  }
  TRUE
}
