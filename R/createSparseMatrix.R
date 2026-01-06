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
#' M <- createSparseMatrix(
#'   N = 30, sparsity = 0.05, method = "normal",
#'   stationary = TRUE
#' )
#' @export
create_sparse_matrix <- function(n, sparsity, method = "normal",
                                 stationary = FALSE, p = 1, ...) {
  opt <- list(...)
  mu <- ifelse(!is.null(opt$mu_mat), opt$mu_mat, 0)
  sd <- ifelse(!is.null(opt$sd_mat), opt$sd_mat, 1)
  non_zero_elms <- floor(sparsity * (n^2))

  if (method == "normal") {
    # normal distributed nonzero entries
    non_zero_entries <- stats::rnorm(non_zero_elms, mean = mu, sd = sd)
    entries <- sample(x = 1:(n^2), size = non_zero_elms, replace = FALSE)
    a_tmp <- numeric(length = (n^2))
    a_tmp[entries] <- non_zero_entries
    a_mat <- matrix(a_tmp, nrow = n, ncol = n)
  } else if (method == "bimodal") {
    # bimodal (bi-normal) distributed nonzero entries
    non_zero_entries_left <- stats::rnorm(non_zero_elms, mean = -mu, sd = sd)
    non_zero_entries_right <- stats::rnorm(non_zero_elms, mean = mu, sd = sd)
    non_zero_entries <- sample(
      x = c(
        non_zero_entries_left,
        non_zero_entries_right
      ),
      size = non_zero_elms, replace = FALSE
    )
    entries <- sample(x = 1:(n^2), size = non_zero_elms, replace = FALSE)
    a_tmp <- numeric(length = (n^2))
    a_tmp[entries] <- non_zero_entries
    a_mat <- matrix(a_tmp, nrow = n, ncol = n)
  } else if (method == "full") {
    # full matrix: used only for tests
    e <- 0.9^(1:n)
    d <- diag(e)
    p <- matrix(0, n, n)
    while (det(p) == 0) {
      p <- create_sparse_matrix(N = n, sparsity = 1, method = "bimodal")
    }
    a_mat <- solve(p) %*% d %*% p
    stationary <- FALSE
  } else {
    # invalid method
    stop("Unknown method. Possible methods are normal or bimodal.")
  }

  if (stationary == TRUE) {
    # if spectral radius < 1 is needed, return the re-normalized matrix
    k <- 1
    a_mat <- 1 / (k * base::sqrt(p * sparsity * n * sd)) * a_mat
  }
  a_mat
}
