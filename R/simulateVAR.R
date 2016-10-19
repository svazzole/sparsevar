#' @title VAR simulation
#'
#' @description This function generates a simulated multivariate VAR time series.
#'
#' @usage simulateVAR(N, p, nobs, rho, sparsity, mu, method, covariance)
#'
#' @param N dimension of the time series.
#' @param p number of lags of the VAR model.
#' @param nobs number of observations to be generated.
#' @param rho base value for the covariance matrix.
#' @param sparsity density (in percentage) of the number of nonzero elements of the VAR matrices.
#' @param mu a vector containing the mean of the simulated process.
#' @param method which method to use to generate the VAR matrix. Possible values
#' are \code{"normal"} or \code{"bimodal"}.
#' @param covariance type of covariance matrix to use in the simulation. Possible
#' values: \code{"toeplitz"}, \code{"block1"}, \code{"block2"} or simply \code{"diagonal"}.
#'
#' @return A a list of NxN matrices ordered by lag
#' @return data a list with two elements: \code{series} the multivariate time series and
#' \code{noises} the time series of errors
#' @return S the variance/covariance matrix of the process
#'
#' @export
simulateVAR <- function(N = 100, p = 1, nobs = 250, rho = 0.5, sparsity = 0.05,
                        mu = 0, method = "normal", covariance = "Toeplitz") {

  # Create the list of the VAR matrices
  A <- list()
  for (i in 1:p) {
    A[[i]] <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE, p = p)
    l <- max(Mod(eigen(A[[i]])$values))
    while ((l > 1) | (l == 0)) {
      A[[i]] <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE, p = p)
      l <- max(Mod(eigen(A[[i]])$values))
    }
    A[[i]] <- 1/sqrt(p) * A[[i]]
  }

  # Covariance Matrix: Toeplitz, Block1 or Block2
  if (covariance == "block1"){

    l <- floor(N/2)
    I <- diag(1 - rho, nrow = N)
    r <- matrix(0, nrow = N, ncol = N)
    r[1:l, 1:l] <- rho
    r[(l+1):N, (l+1):N] <- diag(rho, nrow = (N-l))
    C <- I + r

  } else if (covariance == "block2") {

    l <- floor(N/2)
    I <- diag(1 - rho, nrow = N)
    r <- matrix(0, nrow = N, ncol = N)
    r[1:l, 1:l] <- rho
    r[(l+1):N, (l+1):N] <- rho
    C <- I + r

  } else if (covariance == "Toeplitz"){

    r <- rho^(1:N)
    C <- Matrix::toeplitz(r)

  } else if (covariance == "Wishart"){

    r <- rho^(1:N)
    S <- Matrix::toeplitz(r)
    C <- stats::rWishart(1, 2*N, S)
    C <- as.matrix(C[, , 1])

  } else if (covariance == "diagonal"){

    C <- diag(x = rho, nrow = N, ncol = N)

  } else {

    stop("Unknown covariance matrix type. Possible choices are: toeplitz, block1, block2 or diagonal")

  }

  # Matrix for MA part
  theta <- matrix(0, N, N)
  # ar <- 1:p

  # Generate the VAR process
  data <- generateVARseries(nobs = nobs, mu, AR = A, sigma = C, skip = 200)

  # Output
  out <- list()
  out$A <- A
  out$series <- data$series
  out$noises <- data$noises
  out$sigma <- C

  attr(out, "class") <- "var"
  attr(out, "type") <- "simulation"
  return(out)

}

generateVARseries <- function(nobs, mu, AR, sigma, skip = 200) {

  ## This function creates the simulated time series

  N <- nrow(sigma)
  nT <- nobs + skip
  at <- mvtnorm::rmvnorm(nT, rep(0,N), sigma)

  p <- length(AR)

  ist <- p + 1
  zt <- matrix(0, nT, N)

  if(length(mu)==0) {
    mu <- rep(0,N)
  }

  for (it in ist:nT){
    tmp <- matrix(at[it,], 1, N)

    for (i in 1:p){
      ph <- AR[[i]]
      ztm <- matrix(zt[it-i, ], 1, N)
      tmp <- tmp + ztm%*%t(ph)
    }

    zt[it, ] <- mu + tmp
  }

  # skip the first skip points to initialize the series
  zt <- zt[(1+skip):nT, ]
  at <- at[(1+skip):nT, ]

  out <- list()
  out$series <- zt
  out$noises <- at
  return(out)

}
