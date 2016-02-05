#' @title Create Sparse Matrix
#' 
#' @description Creates a sparse square matrix with a given sparsity and distribution.
#' 
#' @param N the dimension of the square matrix
#' @param sparsity the density of non zero elements
#' @param method the method used to generate the entries of the matrix. Possible values are 
#' \code{"normal"} (default) or \code{"bimodal"}.
#' @param stationary should the spectral radius of the matrix be smaller than 1? 
#' Possible values are \code{TRUE} or \code{FALSE}
#' @return An NxN sparse matrix. 
#' @examples
#' M <- createSparseMatrix(N = 30, sparsity = 0.05, method = "normal", stationary = TRUE)
#' 
#' @author Simone Vazzoler

###########################################################
# What about using rsparsematrix from the library Matrix? #
###########################################################
# Which are the differences between a matrix created with #
# createSparseMatrix and rsparsematrix ?                  #
###########################################################


createSparseMatrix <- function(N, sparsity, method = "normal", stationary = FALSE) {
  
  if (method == "normal") {

#     ##########################################################
#     # Old method    
#     ##########################################################    
#     randEig <- runif(N, min = -1, max = 1)
#     D <- diag(randEig)
#     P <- matrix(rnorm(N^2, mean = 0, sd = 1), N, N)
#     invP <- solve(P)
#     
#     Atmp <- (P %*% D %*% invP)
#     ##########################################################
    
    n <- floor(sparsity * N^2)
    nonZeroEntries <- rnorm(n, mean = 0, sd = 1)
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)
    
  } else if (method == "bimodal") {
    
    n <- floor(sparsity * N^2)
    
    nonZeroEntriesLeft <- rnorm(n, mean = -1, sd = 1)
    nonZeroEntriesRight <- rnorm(n, mean = 1, sd = 1)
    
    nonZeroEntries <- sample(x = c(nonZeroEntriesLeft, nonZeroEntriesRight), size = n, replace = FALSE)
    
    entries <- sample(x = 1:N^2, size = n, replace = FALSE)
    
    Atmp <- numeric(length = N^2)
    Atmp[entries] <- nonZeroEntries
    
    A <- matrix(Atmp, nrow = N, ncol = N)

  } else {
    
    stop("Unknown method. Possible methods are normal or bimodal.")
    
  }

#   ##########################################################  
#   # Old method
#   ##########################################################
#     r <- rbinom(n = N^2, 1, sparsity)
#     rM <- matrix(r, nrow = N, ncol = N)
#     A <- 1/sqrt(sparsity * N) *  Atmp * rM
#   ##########################################################
  
  if (stationary == TRUE){
    return(1/sqrt(N) * A)
  } else {
    return(A)
  }
  
}
