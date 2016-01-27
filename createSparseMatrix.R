
createSparseMatrix <- function(sparsity, N) {
  
  randEig <- runif(N, min = -1, max = 1)
  D <- diag(randEig)
  P <- matrix(rnorm(N^2, mean = 0, sd = 1), N, N)
  invP <- solve(P)
  
  Atmp <- P %*% D %*% invP

  r <- rbinom(n = N^2, 1, sparsity)
  rM <- matrix(r, nrow = N, ncol = N)
  
  A <- Atmp * rM
  
  return(A)
  
}