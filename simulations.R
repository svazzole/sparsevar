# Libraries
library(MTS)

# Other files needed
source("createSparseMatrix.R")

simulateVAR <- function(nobs = 250, rho = 0.75, N = 100, sparsity = 0.05, method = "normal"){
  
  # Create sparse matrix for VAR
  A <- createSparseMatrix(sparsity = sparsity, N = N, method = method)
  while (max(Mod(eigen(A)$values)) > 1) {
    A <- createSparseMatrix(sparsity = sparsity, N = N, method = method)
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
