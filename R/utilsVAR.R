#' @title Transorm data
#'
#' @description Transform the input data
#'
#' @usage transformData(data, p, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param opt a list containing the options
#'
#' @export
transform_data <- function(data, p, opt) {

  # get the number of rows and columns
  nr <- nrow(data)
  nc <- ncol(data)

  # make sure the data is in matrix format
  data <- as.matrix(data)

  # scale the matrix columns
  scale <- ifelse(is.null(opt$scale), FALSE, opt$scale)
  # center the matrix columns (default)
  center <- ifelse(is.null(opt$center), TRUE, opt$center)

  if (center == TRUE) {
    if (opt$method == "timeSlice") {
      leave_out <- ifelse(is.null(opt$leaveOut), 10, opt$leaveOut)
      m <- colMeans(data[1:(nr - leave_out), ])
    } else {
      m <- colMeans(data)
    }
    cm <- matrix(rep(m, nrow(data)), nrow = nrow(data), byrow = TRUE)
    data <- data - cm
  } else {
    m <- rep(0, nc)
  }

  if (scale == TRUE) {
    data <- apply(FUN = scale, X = data, MARGIN = 2)
  }

  # create Xs and Ys (temp variables)
  tmp_x <- data[1:(nr - 1), ]
  tmp_y <- data[2:(nr), ]

  # create the data matrix
  tmp_x <- duplicate_matrix(tmp_x, p)
  tmp_y <- tmp_y[p:nrow(tmp_y), ]

  y <- as.vector(tmp_y)

  # Hadamard product for data
  i <- Matrix::Diagonal(nc)
  x <- kronecker(i, tmp_x)

  output <- list()
  output$X <- x
  output$y <- y
  output$series <- data
  output$mu <- t(m)
  output
}

#' @title VAR ENET
#'
#' @description Estimate VAR using ENET penalty
#'
#' @usage varENET(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#'
#' @export
var_enet <- function(data, p, lambdas, opt) {
  # transform the dataset
  tr_dt <- transform_data(data, p, opt)
  fit <- glmnet::glmnet(tr_dt$X, tr_dt$y, lambda = lambdas)
  fit
}

#' @title VAR SCAD
#'
#' @description Estimate VAR using SCAD penalty
#'
#' @usage varSCAD(data, p, lambdas, opt, penalty)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#' @param penalty a string "SCAD" or something else
#'
#' @export

var_scad <- function(data, p, lambdas, opt, penalty = "SCAD") {
  # transform the dataset
  tr_dt <- transform_data(data, p, opt)

  if (penalty == "SCAD") {
    fit <- ncvreg::ncvreg(as.matrix(tr_dt$X), tr_dt$y,
      family = "gaussian", penalty = "SCAD",
      alpha = 1, lambda = lambdas
    )
  } else {
    stop("[WIP] Only SCAD regression is supported at the moment")
  }
  fit
}

#' @title VAR MCP
#'
#' @description Estimate VAR using MCP penalty
#'
#' @usage varMCP(data, p, lambdas, opt)
#'
#' @param data the data
#' @param p the order of the VAR
#' @param lambdas a vector containing the lambdas to be used in the fit
#' @param opt a list containing the options
#'
#' @export
var_mcp <- function(data, p, lambdas, opt) {
  # transform the dataset
  tr_dt <- transform_data(data, p, opt)
  fit <- ncvreg::ncvreg(as.matrix(tr_dt$X), tr_dt$y,
    family = "gaussian", penalty = "MCP",
    alpha = 1, lambda = lambdas
  )
  fit
}

split_matrix <- function(m, p) {
  nr <- nrow(m)
  a <- list()
  for (i in 1:p) {
    ix <- ((i - 1) * nr) + (1:nr)
    a[[i]] <- m[1:nr, ix]
  }
  a
}

duplicate_matrix <- function(data, p) {
  nr <- nrow(data)
  nc <- ncol(data)

  output_data <- data

  if (p > 1) {
    for (i in 1:(p - 1)) {
      tmp_data <- matrix(0, nrow = nr, ncol = nc)
      tmp_data[(i + 1):nr, ] <- data[1:(nr - i), ]
      output_data <- cbind(output_data, tmp_data)
    }
  }

  output_data <- output_data[p:nr, ]
  output_data
}

compute_residuals <- function(data, a) {
  nr <- nrow(data)
  nc <- ncol(data)
  p <- length(a)

  res <- matrix(0, ncol = nc, nrow = nr)
  f <- matrix(0, ncol = nc, nrow = nr)

  for (i in 1:p) {
    tmp_d <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    tmp_f <- t(a[[i]] %*% t(tmp_d))
    f <- f + tmp_f
  }

  res <- data - f
  res
}

#' @title Companion VAR
#'
#' @description Build the VAR(1) representation of a VAR(p) process
#'
#' @usage companionVAR(v)
#'
#' @param v the VAR object as from \code{fitVAR} or \code{simulateVAR}
#'
#' @export
companion_var <- function(v) {
  if (!checkIsVar(v)) {
    stop("v must be a var object")
  }
  a <- v$A
  nc <- ncol(a[[1]])
  p <- length(a)
  if (p > 1) {
    big_a <- Matrix::Matrix(0, nrow = p * nc, ncol = p * nc, sparse = TRUE)
    for (k in 1:p) {
      ix <- ((k - 1) * nc) + (1:nc)
      big_a[1:nc, ix] <- a[[k]]
    }
    ix_r <- (nc + 1):nrow(big_a)
    ix_c <- 1:((p - 1) * nc)
    big_a[ix_r, ix_c] <- diag(1, nrow = length(ix_c), ncol = length(ix_c))
  } else {
    big_a <- Matrix::Matrix(a[[1]], sparse = TRUE)
  }
  big_a
}

#' @title Bootstrap VAR
#'
#' @description Build the bootstrapped series from the original var
#'
#' @usage bootstrappedVAR(v)
#'
#' @param v the VAR object as from fitVAR or simulateVAR
#'
#' @export
bootstrapped_var <- function(v) {

  ## This function creates the bootstrapped time series
  if (!checkIsVar(v)) {
    stop("v must be a var object")
  }

  r <- v$residuals
  s <- v$series
  a <- v$A
  n <- ncol(a[[1]])
  p <- length(a)
  t <- nrow(r)
  r <- r - matrix(colMeans(r), ncol = n, nrow = t)

  zt <- matrix(0, nrow = t, ncol = n)
  zt[1:p, ] <- s[1:p, ]

  for (t0 in (p + 1):t) {
    ix <- sample((p + 1):t, 1)
    u <- r[ix, ]
    vv <- rep(0, n)
    for (i in 1:p) {
      ph <- a[[i]]
      vv <- vv + ph %*% zt[(t0 - i), ]
    }
    vv <- vv + u
    zt[t0, ] <- vv
  }
  zt
}

#' @title Test for Ganger Causality
#'
#' @description This function should retain only the coefficients of the
#' matrices of the VAR that are statistically significative (from the bootstrap)
#'
#' @usage testGranger(v, eb)
#'
#' @param v the VAR object as from fitVAR or simulateVAR
#' @param eb the error bands as obtained from errorBands
#'
#' @export
test_granger <- function(v, eb) {
  p <- length(v$A)
  a <- list()
  for (i in 1:p) {
    l <- (eb$irfQUB[, , i + 1] >= 0 & eb$irfQLB[, , i + 1] <= 0)
    a[[i]] <- v$A[[i]] * (1 - l)
  }
  a
}

#' @title Computes information criteria for VARs
#'
#' @description This function computes information criteria (AIC, Schwartz and
#' Hannan-Quinn) for VARs.
#'
#' @usage informCrit(v)
#'
#' @param v a list of VAR objects as from fitVAR.
#'
#' @export
inform_crit <- function(v) {
  if (is.list(v)) {
    k <- length(v)
    r <- matrix(0, nrow = k, ncol = 3)
    for (i in 1:k) {
      if (attr(v[[1]], "class") == "var" || attr(v[[1]], "class") == "vecm") {
        p <- length(v[[i]]$A)
        # Compute sparsity
        s <- 0
        for (l in 1:p) {
          s <- s + sum(v[[i]]$A[[l]] != 0)
        }
        sp <- s / (p * ncol(v[[i]]$A[[1]])^2)
      } else {
        stop("List elements must be var or vecm objects.")
      }
      sigma <- v[[i]]$sigma
      nr <- nrow(v[[i]]$residuals)
      nc <- ncol(v[[i]]$residuals)
      d <- det(sigma)

      r[i, 1] <- log(d) + (2 * p * sp * nc^2) / nr # AIC
      r[i, 2] <- log(d) + (log(nr) / nr) * (p * sp * nc^2) # BIC
      r[i, 3] <- log(d) + (2 * p * sp * nc^2) / nr * log(log(nr)) # Hannan-Quinn
    }
    results <- data.frame(r)
    colnames(results) <- c("AIC", "BIC", "HannanQuinn")
  } else {
    stop("Input must be a list of var models.")
  }
  results
}

estimate_covariance <- function(res, ...) {
  nc <- ncol(res)
  s <- corpcor::cov.shrink(res, verbose = FALSE)
  sigma <- matrix(0, ncol = nc, nrow = nc)
  for (i in 1:nc) {
    for (j in 1:nc) {
      sigma[i, j] <- s[i, j]
    }
  }
  sigma
}

#' @title Computes forecasts for VARs
#'
#' @description This function computes forecasts for a given VAR.
#'
#' @usage computeForecasts(v, num_steps)
#'
#' @param v a VAR object as from fitVAR.
#' @param num_steps the number of forecasts to produce.
#'
#' @export
compute_forecasts <- function(v, num_steps = 1) {
  if (!checkIsVar(v)) {
    stop("You must pass a var object.")
  } else {
    mu <- v$mu
    data <- v$series
    v <- v$A
  }

  if (!is.list(v)) {
    stop("v must be a var object or a list of matrices.")
  } else {
    nr <- nrow(data)
    nc <- ncol(v[[1]])
    p <- length(v)

    f <- matrix(0, nrow = nc, ncol = num_steps)

    tmp_data <- matrix(data = t(data[(nr - p + 1):nr, ]),
                       nrow = nc,
                       ncol = num_steps)
    nr <- ncol(tmp_data)

    for (n in 1:num_steps) {
      for (k in 1:p) {
        if (n == 1) {
          f[, n] <- f[, n] + v[[k]] %*% tmp_data[, nr - k + 1]
        } else {
          if (nr > 1) {
            tmp_data <- cbind(tmp_data[, 2:nr], f[, n - 1])
          } else {
            tmp_data <- as.matrix(f[, n - 1])
          }
          f[, n] <- f[, n] + v[[k]] %*% tmp_data[, nr - k + 1]
        }
      }
    }
  }
  f <- f + matrix(rep(mu, length(mu)), length(mu), num_steps)
  f
}

apply_threshold <- function(a_mat, nr, nc, p, type = "soft") {
  if (type == "soft") {
    tr <- 1 / sqrt(p * nc * log(nr))
  } else if (type == "hard") {
    tr <- (nc) ^ (-0.49)
  } else {
    stop("Unknown threshold type. Possible values are: \"soft\" or \"hard\"")
  }

  l_mat <- abs(a_mat) >= tr
  a_mat * l_mat
}
