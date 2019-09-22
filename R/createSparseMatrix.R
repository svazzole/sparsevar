#' @title Create Sparse Matrix
#'
#' @description Creates a sparse square matrix with a given sparsity and
#' distribution.
#'
#' @param N the dimension of the square matrix
#' @param sparsity the density of non zero elements
#' @param method the method used to generate the entries of the matrix.
#' Possible values are \code{"normal"} (default) or \code{"bimodal"}.
#' @param stationary should the spectral radius of the matrix be smaller than 1?
#' Possible values are \code{TRUE} or \code{FALSE}. Default is \code{FALSE}.
#' @param p normalization constant (used for VAR of order greater than 1,
#' default = 1)
#' @param ... other options for the matrix (you can specify the mean
#' \code{mu_mat} and the standard deviation \code{sd_mat}).
#' @return An NxN sparse matrix.
#' @examples
#' M <- createSparseMatrix(N = 30, sparsity = 0.05, method = "normal",
#' stationary = TRUE)
#' @export
createSparseMatrix <- function(N, sparsity, method = "normal",
                               stationary = FALSE, p = 1, ...) {
  opt <- list(...)
  mu <- ifelse(!is.null(opt$mu_mat), opt$mu_mat, 0)
  sd <- ifelse(!is.null(opt$sd_mat), opt$sd_mat, 1)
  n <- floor(sparsity * (N ^ 2))

  if (method == "normal") {
    # normal distributed nonzero entries
    non_zero_entries <- stats::rnorm(n, mean = mu, sd = sd)
    entries <- sample(x = 1:(N ^ 2), size = n, replace = FALSE)
    Atmp <- numeric(length = (N ^ 2))
    Atmp[entries] <- non_zero_entries
    A <- matrix(Atmp, nrow = N, ncol = N)
  } else if (method == "bimodal") {
    # bimodal (bi-normal) distributed nonzero entries
    non_zero_entries_left <- stats::rnorm(n, mean = -mu, sd = sd)
    non_zero_entries_right <- stats::rnorm(n, mean = mu, sd = sd)
    non_zero_entries <- sample(x = c(non_zero_entries_left,
                                   non_zero_entries_right),
                             size = n, replace = FALSE)
    entries <- sample(x = 1:(N ^ 2), size = n, replace = FALSE)
    Atmp <- numeric(length = (N ^ 2))
    Atmp[entries] <- non_zero_entries
    A <- matrix(Atmp, nrow = N, ncol = N)
  } else if (method == "full") {
    # full matrix: used only for tests
    e <- 0.9 ^ (1:N)
    D <- diag(e)
    P <- matrix(0, N, N)
    while (det(P) == 0) {
      P <- createSparseMatrix(N = N, sparsity = 1, method = "bimodal")
    }
    A <- solve(P) %*% D %*% P
    stationary <- FALSE
  } else {
    # invalid method
    stop("Unknown method. Possible methods are normal or bimodal.")
  }

  if (stationary == TRUE) {
    # if spectral radius < 1 is needed, return the re-normalized matrix
    K <- 1
    return(1 / (K * base::sqrt(p * sparsity * N * sd)) * A)
  } else {
    return(A)
  }
}
