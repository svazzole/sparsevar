#' @title VAR simulation
#'
#' @description This function generates a multivariate VAR time series
#' @param N
#' @param nobs
#' @param rho
#' @param sparsity
#' @param method
#' @param covariance
#' 
#' @return A list containing ...
#' 
#' @author Simone Vazzoler
#'
#' @export
simulateVAR <- function(N = 100, nobs = 250, rho = 0.5, sparsity = 0.05, method = "normal", covariance = "toeplitz"){
  
 # Create sparse matrix for VAR
  A <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
  
  while (max(Mod(eigen(A)$values)) > 1) {
    
    A <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
  
  }

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
  
  # Generate VAR(1) process 
  data <- VARMAsim(nobs = nobs, arlags = 1, malags = 0, cnst = 0, phi = A, theta = theta, skip = 200, sigma = T)
  
  out <- list()
  out$data <- data
  out$A <- A
  
  return(out)
  
}
