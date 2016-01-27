library(MTS)

simulateVAR <- function(nobs = 250, rho = 0.75, N = 100, sparsity = 0.03){
  
  source("createSparseMatrix.R")

  # Create sparse matrix for VAR
  A <- createSparseMatrix(sparsity = sparsity, N = N)
  while (max(abs(eigen(A)$values)) > 1) {
    A <- createSparseMatrix(sparsity = sparsity, N = N)
  }
  
  # Covariance Matrix: Toeplitz
  r <- rho^(1:N)
  T <- toeplitz(r) 
  
  # Matrix for MA part
  theta <- matrix(0, N, N)
  
  # Generate VAR(1) process 
  data <- VARMAsim(nobs = nobs, arlags = 1, malags = 0, cnst = 0, phi = A, theta = theta, skip = 200, sigma = T)
  
  final <- list()
  final$data <- data
  final$A <- A
  
  return(final)
  
}
