###########################################################
# What about using rsparsematrix from the library Matrix? #
###########################################################
# Which are the differences between a matrix created with #
# createSparseMatrix and rsparsematrix ?                  #
###########################################################

createSparseMatrix <- function(N, sparsity, method = "normal", stationary = FALSE) {
  
  if (method == "normal") {
    
#     randEig <- runif(N, min = -1, max = 1)
#     D <- diag(randEig)
#     P <- matrix(rnorm(N^2, mean = 0, sd = 1), N, N)
#     invP <- solve(P)
#     
#     Atmp <- (P %*% D %*% invP)

    n <- floor(sparsity * N^2)
    nonZeroEntries <- rnorm(n, mean = 0, sd = 1)
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)
#    return(A)
    
  } else if (method == "bimodal") {
    
    n <- floor(sparsity * N^2)
    
#     t <- c(rnorm(N^2, mean = -1, sd = 1), rnorm(N^2, mean = 1, sd = 1))
#     t <- sample(t, N^2)
#     Atmp <- (matrix(t, nrow = N, ncol = N))
    nonZeroEntriesLeft <- rnorm(n, mean = -1, sd = 1)
    nonZeroEntriesRight <- rnorm(n, mean = 1, sd = 1)
    
    nonZeroEntries <- sample(x = c(nonZeroEntriesLeft, nonZeroEntriesRight), size = n, replace = FALSE)
    
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)
#    return(A)
    
  } else {
    
    stop("Unknown method. Possible methods are normal or bimodal.")
    
  }
  
#   r <- rbinom(n = N^2, 1, sparsity)
#   rM <- matrix(r, nrow = N, ncol = N)
#   A <- 1/sqrt(sparsity * N) *  Atmp * rM
  if (stationary == TRUE){
    return(1/sqrt(N) * A)
  } else {
    return(A)
  }
  
}
