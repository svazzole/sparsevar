#' @title VAR simulation
#'
#' @description This function generates a simulated multivariate VAR time series.
#' 
#' @param N dimension of the time series. 
#' @param p number of lags of the VAR model.
#' @param nobs number of observations to be generated.
#' @param rho base value for the covariance matrix.
#' @param sparsity density (in percentage) of the number of nonzero elements of the VAR matrices.
#' @param method which method to use to generate the VAR matrix. Possible values
#' are \code{"normal"} or \code{"bimodal"}.
#' @param covariance type of covariance matrix to use in the simulation. Possible 
#' values: \code{"toeplitz"}, \code{"block1"} or \code{"block2"}.
#' 
#' @return A list containing ...
#' 
#' @author Simone Vazzoler
#'
#' @export
#' 
simulateVAR <- function(N = 100, nobs = 250, rho = 0.5, sparsity = 0.05, p = 1, method = "normal", covariance = "toeplitz") {
  
 # Create sparse matrices for VAR
  # if (p==1) {
  # 
  #   # only 1 lag
  #   A <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
  #   while (max(Mod(eigen(A)$values)) > 1) {
  #     A <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
  #   }
  # 
  # } else {
    # p lags
    A <- list()
    cA <- matrix(0, nrow = N, ncol = N * p)
    for (i in 1:p) {
      A[[i]] <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
      while (max(Mod(eigen(A[[i]])$values)) > 1) {
        A[[i]] <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
      }
      cA[1:N, ((i-1) * N) + (1:N)] <- A[[i]]
    }
    
  # }

  # Covariance Matrix: Toeplitz, Block1 or Block2
  if (covariance == "block1"){
    
    l <- floor(N/2)
    I <- diag(1 - rho, nrow = N)
    r <- matrix(0, nrow = N, ncol = N)
    r[1:l, 1:l] <- rho
    T <- I + r
      
  } else if (covariance == "block2") {
  
    l <- floor(N/2)
    I <- diag(1 - rho, nrow = N)
    r <- matrix(0, nrow = N, ncol = N)
    r[1:l, 1:l] <- rho
    r[(l+1):N, (l+1):N] <- rho
    T <- I + r
      
  } else if (covariance == "toeplitz"){
    
    r <- rho^(1:N)
    T <- toeplitz(r) 
  
  } else {
    
    stop("Unknown covariance matrix type. Possible choices are: toeplitz, block1, block2")
    
  }
  
  # Matrix for MA part
  theta <- matrix(0, N, N)
  
  ar <- 1:p
  
  # Generate VAR(1) process 
  data <- MTS::VARMAsim(nobs = nobs, arlags = ar, malags = 0, cnst = 0, phi = cA, theta = theta, skip = 200, sigma = T)
  
  out <- list()
  out$data <- data
  out$A <- A
  out$S <- T
  
  return(out)
  
}
