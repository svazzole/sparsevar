
createSparseMatrix <- function(N, sparsity, method = "normal") {
  
  if (method == "normal") {
    
    randEig <- runif(N, min = -1, max = 1)
    D <- diag(randEig)
    P <- matrix(rnorm(N^2, mean = 0, sd = 1), N, N)
    invP <- solve(P)
    
    Atmp <- (P %*% D %*% invP)

    
  } else if (method == "bimodal") {
    
    t <- c(rnorm(N^2, mean = -1, sd = 1), rnorm(N^2, mean = 1, sd = 1))
    t <- sample(t, N^2)
    Atmp <- (matrix(t, nrow = N, ncol = N))
     
  } else {
    
    stop("Unknown method. Possible methods are normal or bimodal.")
    
  }
  
  r <- rbinom(n = N^2, 1, sparsity)
  rM <- matrix(r, nrow = N, ncol = N)
  A <- 1/sqrt(sparsity * N) *  Atmp * rM
    
  return(A)
  
}