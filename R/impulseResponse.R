#' @title Impulse Response Function
#'
#' @description A function to estimate the Impulse Response Function of a
#' given VAR.
#'
#' @usage impulse_response(v, len = 20)
#'
#' @param v the data in the for of a VAR
#' @param len length of the impulse response function
#'
#' @return \code{irf} a 3d array containing the impulse response function.
#'
#' @export
impulse_response <- function(v, len = 20) {

  # Check if v is a VAR object
  if (!check_is_var(v)) {
    stop("Input v must be a VAR object")
  }

  # Numerical problems in the estimated variance covariance
  e <- eigen(v$sigma)$values
  if (any(e <= 0)) {
    p <- t(chol(v$sigma, pivot = TRUE))
  } else {
    p <- t(chol(v$sigma))
  }

  big_a <- companion_var(v)

  out <- get_irf(v, big_a, len = len, p)
  out$cholP <- p # Add Choleski factorization to the output
  out
}

#' @title Check Impulse Zero
#'
#' @description A function to find which entries of the impulse
#' response function are zero.
#'
#' @usage check_impulse_zero(irf)
#'
#' @param irf irf output from impulseResponse function
#'
#' @return a matrix containing the indices of the impulse response function that
#' are 0.
#'
#' @export
check_impulse_zero <- function(irf) {
  nx <- dim(irf)[1]
  ny <- dim(irf)[2]
  nz <- dim(irf)[3]
  logical_irf <- matrix(0, nx, ny)
  for (z in 1:nz) {
    logical_irf <- logical_irf + abs(irf[, , z])
  }
  logical_irf <- logical_irf == 0
  which(logical_irf == TRUE, arr.ind = TRUE)
}

#' @title Error bands for IRF
#'
#' @description A function to estimate the confidence intervals for
#' irf and oirf.
#'
#' @usage error_bands_irf(v, irf, alpha, m, resampling, ...)
#'
#' @param v a var object as from fitVAR or simulateVAR
#' @param irf irf output from impulseResponse function
#' @param alpha level of confidence (default \code{alpha = 0.01})
#' @param m number of bootstrapped series (default \code{m = 100})
#' @param resampling type of resampling: \code{"bootstrap"} or
#' \code{"jackknife"}
#' @param ... some options for the estimation: \code{verbose = TRUE}
#' or \code{FALSE}, \code{mode = "fast"} or \code{"slow"},
#' \code{threshold = TRUE} or \code{FALSE}.
#'
#' @return a matrix containing the indices of the impulse response function that
#' are 0.
#'
#' @export
error_bands_irf <- function(v, irf, alpha = 0.01, m = 100,
                            resampling = "bootstrap", ...) {
  opt <- list(...)
  verbose <- ifelse(!is.null(opt$verbose), opt$verbose, TRUE)
  mode <- ifelse(!is.null(opt$mode), opt$mode, "fast")
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  threshold_type <- ifelse(!is.null(opt$thresholdType),
                           opt$thresholdType, "soft")

  if (resampling == "bootstrap") {
    lambda <- v$lambda
    p <- length(v$A)
    nr <- ncol(v$series)
    nc <- ncol(v$A[[1]])
    len <- dim(irf$irf)[3]

    irfs <- array(data = rep(0, len * nc^2 * m), dim = c(nc, nc, len + 1, m))
    oirfs <- array(data = rep(0, len * nc^2 * m), dim = c(nc, nc, len + 1, m))

    if (verbose == TRUE) {
      cat("Step 1 of 2: bootstrapping series and re-estimating VAR...\n")
      pb <- utils::txtProgressBar(min = 0, max = m, style = 3)
    }

    for (k in 1:m) {
      # create Xs and Ys (temp variables)
      o <- bootstrapped_var(v)

      if (mode == "fast") {
        if (v$penalty == "ENET") {
          # fit ENET to a specific value of lambda
          fit <- var_enet(o, p, lambda,
                          opt = list(method = v$method, penalty = v$penalty))
          a_vector <- stats::coef(fit, s = lambda)
          a <- matrix(a_vector[2:length(a_vector)], nrow = nc,
                      ncol = nc * p, byrow = TRUE)
        } else if (v$penalty == "SCAD") {
          fit <- var_scad(o, p, lambda,
                          opt = list(method = v$method, penalty = v$penalty))
          a_vector <- fit$beta[2:nrow(fit$beta), 1]
          a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)
        } else {
          fit <- var_mcp(o, p, lambda,
                         opt = list(method = v$method, penalty = v$penalty))
          a_vector <- fit$beta[2:nrow(fit$beta), 1]
          a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)
        }

        if (threshold) {
          apply_threshold(a, nr, nc, p, type = threshold_type)
        }

        m_mat <- cbind(diag(x = 1, nrow = (nc * (p - 1)),
                            ncol = (nc * (p - 1))),
                       matrix(0, nrow = (nc * (p - 1)), ncol = nc))
        big_a <- rbind(a, m_mat)
      } else {
        # fit ENET on a series of lambdas
        if (threshold) {
          fit <- fit_var(o, p, penalty = v$penalty,
                         method = v$method, threshold = TRUE)
        } else {
          fit <- fit_var(o, p, penalty = v$penalty, method = v$method)
        }
        big_a <- companion_var(fit)
      }

      tmp_res <- get_irf(v, big_a, len, irf$cholP)
      irfs[, , , k] <- tmp_res$irf
      oirfs[, , , k] <- tmp_res$oirf

      if (verbose == TRUE) {
        utils::setTxtProgressBar(pb, k)
      }
    }

    if (verbose == TRUE) {
      close(pb)
      cat("Step 2 of 2: computing quantiles...\n")
      pb <- utils::txtProgressBar(min = 0, max = (nc * nc), style = 3)
    }

    irf_ub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
    irf_lb <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
    oirf_ub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
    oirf_lb <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
    irf_qub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
    irf_qlb <- irf_qub
    oirf_qub <- irf_qub
    oirf_qlb <- irf_qub

    a <- alpha / 2
    q_lb <- stats::qnorm(a)
    q_ub <- stats::qnorm((1 - a))

    for (i in 1:nc) {
      for (j in 1:nc) {
        for (k in 1:len) {
          irf_qub[i, j, k] <- stats::quantile(irfs[i, j, k, ], probs = (1 - a),
                                              na.rm = TRUE)
          oirf_qub[i, j, k] <- stats::quantile(oirfs[i, j, k, ],
                                               probs = (1 - a), na.rm = TRUE)
          irf_qlb[i, j, k] <- stats::quantile(irfs[i, j, k, ], probs = a,
                                              na.rm = TRUE)
          oirf_qlb[i, j, k] <- stats::quantile(oirfs[i, j, k, ],
                                               probs = a, na.rm = TRUE)

          irf_ub[i, j, k] <- q_ub * stats::sd(irfs[i, j, k, ])
          oirf_ub[i, j, k] <- q_ub * stats::sd(oirfs[i, j, k, ])
          irf_lb[i, j, k] <- q_lb * stats::sd(irfs[i, j, k, ])
          oirf_lb[i, j, k] <- q_lb * stats::sd(oirfs[i, j, k, ])
        }
        if (verbose == TRUE) {
          utils::setTxtProgressBar(pb, (i - 1) * nc + j)
        }
      }
    }

    if (verbose == TRUE) {
      close(pb)
    }

    output <- list()

    output$irfUB <- irf_ub
    output$oirfUB <- oirf_ub
    output$irfLB <- irf_lb
    output$oirfLB <- oirf_lb

    output$irfQUB <- irf_qub
    output$oirfQUB <- oirf_qub
    output$irfQLB <- irf_qlb
    output$oirfQLB <- oirf_qlb

    attr(output, "class") <- "irfBands"
    attr(output, "resampling") <- "bootstrap"
  } else if (resampling == "jackknife") {
    output <- jackknife(v, irf, alpha = alpha, ...)
  } else if (resampling == "bootstrapOLS") {
    output <- bootstrap_ols(v, irf, alpha = alpha, ...)
  } else {
    stop("Unknown resampling method. Possible values are
         \"bootstrap\" or \"jackknife\"")
  }
  output
}

jackknife <- function(v, irf, mode = "fast", alpha, ...) {
  lambda <- v$lambda
  p <- length(v$A)
  nc <- ncol(v$A[[1]])
  len <- dim(irf$irf)[3]
  nr <- nrow(v$series)

  opt <- list(...)
  verbose <- ifelse(!is.null(opt$verbose), opt$verbose, TRUE)
  mode <- ifelse(!is.null(opt$mode), opt$mode, "fast")
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  threshold_type <- ifelse(!is.null(opt$thresholdType),
                           opt$thresholdType, "soft")

  irfs <- array(data = rep(0, len * nc^2 * nr), dim = c(nc, nc, len + 1, nr))
  oirfs <- array(data = rep(0, len * nc^2 * nr), dim = c(nc, nc, len + 1, nr))

  if (verbose == TRUE) {
    cat("Step 1 of 2: jack knifing series and re-estimating VAR...\n")
    pb <- utils::txtProgressBar(min = 0, max = nr, style = 3)
  }

  for (k in 1:nr) {
    # create Xs and Ys (temp variables)
    data <- v$series[-k, ]
    tr_dt <- transform_data(data, p,
                            opt = list(method = v$method, penalty = v$penalty))
    tr_dt$X <- tr_dt$X
    tr_dt$y <- tr_dt$y

    # Things TODO: data <- v$series
    # Things TODO: AtrDt <- transformData(data, p,
    # Things TODO: opt = list(method = v$method, penalty = v$penalty))
    # Things TODO: trDt$X <- trDt$X[-k, ]
    # Things TODO: trDt$y <- trDt$y[-k]

    if (mode == "fast") {
      if (v$penalty == "ENET") {
        # fit ENET to a specific value of lambda
        fit <- var_enet(data, p, lambda,
                        opt = list(method = v$method, penalty = v$penalty))
        a_vector <- stats::coef(fit, s = lambda)
        a <- matrix(a_vector[2:length(a_vector)], nrow = nc,
                    ncol = nc * p, byrow = TRUE)
      } else if (v$penalty == "SCAD") {
        fit <- var_scad(data, p, lambda,
                        opt = list(method = v$method, penalty = v$penalty))
        a_vector <- fit$beta[2:nrow(fit$beta), 1]
        a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)
      } else {
        fit <- var_mcp(data, p, lambda,
                       opt = list(method = v$method, penalty = v$penalty))
        a_vector <- fit$beta[2:nrow(fit$beta), 1]
        a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)
      }

      if (threshold) {
        a <- apply_threshold(a, nr, nc, p, type = threshold_type)
      }

      m <- cbind(diag(x = 1, nrow = (nc * (p - 1)), ncol = (nc * (p - 1))),
                 matrix(0, nrow = (nc * (p - 1)), ncol = nc))
      big_a <- rbind(a, m)
    } else {
      # fit ENET on a series of lambdas
      if (threshold) {
        fit <- fit_var(data, p, penalty = v$penalty,
                       method = v$method, threshold = TRUE)
      } else {
        fit <- fit_var(data, p, penalty = v$penalty, method = v$method)
      }
      big_a <- companion_var(fit)
    }

    tmp_res <- get_irf(v, big_a, len, irf$cholP)
    irfs[, , , k] <- tmp_res$irf
    oirfs[, , , k] <- tmp_res$oirf

    if (verbose == TRUE) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  if (verbose == TRUE) {
    close(pb)
    cat("Step 2 of 2: computing quantiles...\n")
    pb <- utils::txtProgressBar(min = 0, max = (nc * nc), style = 3)
  }

  irf_ub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
  irf_lb <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
  oirf_ub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
  oirf_lb <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))

  a <- alpha / 2
  q_lb <- stats::qnorm(a)
  q_ub <- stats::qnorm((1 - a))

  for (i in 1:nc) {
    for (j in 1:nc) {
      for (k in 1:len) {
        irf_ub[i, j, k] <- base::mean(irfs[i, j, k, ]) +
          q_ub * stats::sd(irfs[i, j, k, ])
        oirf_ub[i, j, k] <- base::mean(oirfs[i, j, k, ]) +
          q_ub * stats::sd(oirfs[i, j, k, ])
        irf_lb[i, j, k] <- base::mean(irfs[i, j, k, ]) +
          q_lb * stats::sd(irfs[i, j, k, ])
        oirf_lb[i, j, k] <- base::mean(oirfs[i, j, k, ]) +
          q_lb * stats::sd(oirfs[i, j, k, ])
      }
      if (verbose == TRUE) {
        utils::setTxtProgressBar(pb, (i - 1) * nc + j)
      }
    }
  }

  if (verbose == TRUE) {
    close(pb)
  }

  output <- list()

  output$irfUB <- irf_ub
  output$oirfUB <- oirf_ub
  output$irfLB <- irf_lb
  output$oirfLB <- oirf_lb

  attr(output, "class") <- "irfBands"
  attr(output, "resampling") <- "jackknife"
  output
}


bootstrap_ols <- function(v, irf, mode = "fast", alpha, ...) {
  # If you want lambda: lambda <- v$lambda
  p <- length(v$A)
  nc <- ncol(v$A[[1]])
  len <- dim(irf$irf)[3]
  nr <- nrow(v$series)

  opt <- list(...)
  verbose <- ifelse(!is.null(opt$verbose), opt$verbose, TRUE)
  mode <- ifelse(!is.null(opt$mode), opt$mode, "fast")
  # Threshold: threshold <- ifelse(!is.null(opt$threshold),
  #                                opt$threshold, FALSE)
  # and its type: thresholdType <- ifelse(!is.null(opt$thresholdType),
  #                                       opt$thresholdType, "soft")

  irfs <- array(data = rep(0, len * nc^2 * nr), dim = c(nc, nc, len + 1, nr))
  oirfs <- array(data = rep(0, len * nc^2 * nr), dim = c(nc, nc, len + 1, nr))

  if (verbose == TRUE) {
    cat("Step 1 of 2: bootstrappingOLS series and re-estimating VAR...\n")
    pb <- utils::txtProgressBar(min = 0, max = nr, style = 3)
  }

  for (k in 1:nr) {
    # create Xs and Ys (temp variables)
    o <- bootstrapped_var(v)

    n <- ncol(v$A[[1]])
    nobs <- nrow(v$series)

    big_a <- companion_var(v)

    # Not needed: tr_dt <- transform_data(o, p = p,
    #                                     opt = list(method = v$method,
    #                                     scale = FALSE, center = TRUE))

    non_zero_entries <- as.matrix(big_a != 0)

    ## Create matrix R
    t <- as.vector(non_zero_entries)
    n <- sum(t != 0)
    ix <- which(t != 0)
    j <- 1:n

    r <- matrix(0, ncol = n, nrow = length(t))
    for (zz in 1:n) {
      r[ix[zz], j[zz]] <- 1
    }

    # Not needed: x <- as.matrix(tr_dt$X)
    y <- as.vector(t(o[-(1:p), ]))

    # Metodo A MANO
    s <- corpcor::invcov.shrink(v$residuals, verbose = FALSE)
    g <- t(o[-nobs, ]) %*% o[-nobs, ] / nobs

    v0 <- solve(t(r) %*% (kronecker(g, s) %*% r))
    vv <- non_zero_entries
    vv[non_zero_entries] <- diag(v0)
    g1 <- solve(t(r) %*% (kronecker(t(o[-nobs, ]) %*% o[-nobs, ], s)) %*% r)
    g2 <- t(r) %*% (kronecker(t(o[-nobs, ]), s))

    g <- g1 %*% g2 # [ , (N+1):(length(y) + N)]
    ga <- g %*% y

    b1 <- vector(length = n * n)
    b1 <- r %*% ga
    a <- matrix(b1, ncol = n, byrow = FALSE)


    m <- cbind(diag(x = 1, nrow = (nc * (p - 1)), ncol = (nc * (p - 1))),
               matrix(0, nrow = (nc * (p - 1)), ncol = nc))
    big_a <- rbind(a, m)

    tmp_res <- get_irf(v0, big_a, len, irf$cholP)
    irfs[, , , k] <- tmp_res$irf
    oirfs[, , , k] <- tmp_res$oirf

    if (verbose == TRUE) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  if (verbose == TRUE) {
    close(pb)
    cat("Step 2 of 2: computing quantiles...\n")
    pb <- utils::txtProgressBar(min = 0, max = (nc * nc), style = 3)
  }

  irf_ub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
  irf_lb <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
  oirf_ub <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))
  oirf_lb <- array(data = rep(0, len * nc^2), dim = c(nc, nc, len))

  a <- alpha / 2
  q_lb <- stats::qnorm(a)
  q_ub <- stats::qnorm((1 - a))

  for (i in 1:nc) {
    for (j in 1:nc) {
      for (k in 1:len) {
        irf_ub[i, j, k] <- base::mean(irfs[i, j, k, ]) +
          q_ub * stats::sd(irfs[i, j, k, ])
        oirf_ub[i, j, k] <- base::mean(oirfs[i, j, k, ]) +
          q_ub * stats::sd(oirfs[i, j, k, ])
        irf_lb[i, j, k] <- base::mean(irfs[i, j, k, ]) +
          q_lb * stats::sd(irfs[i, j, k, ])
        oirf_lb[i, j, k] <- base::mean(oirfs[i, j, k, ]) +
          q_lb * stats::sd(oirfs[i, j, k, ])
      }
      if (verbose == TRUE) {
        utils::setTxtProgressBar(pb, (i - 1) * nc + j)
      }
    }
  }

  if (verbose == TRUE) {
    close(pb)
  }

  output <- list()

  output$irfUB <- irf_ub
  output$oirfUB <- oirf_ub
  output$irfLB <- irf_lb
  output$oirfLB <- oirf_lb

  attr(output, "class") <- "irfBands"
  attr(output, "resampling") <- "jackknife"
  output
}


get_irf <- function(v, big_a, len = 20, p) {
  nr <- nrow(v$A[[1]])

  irf <- array(data = rep(0, len * nr^2), dim = c(nr, nr, len + 1))
  oirf <- array(data = rep(0, len * nr^2), dim = c(nr, nr, len + 1))

  a_tmp <- diag(nrow = nrow(big_a), ncol = ncol(big_a))

  irf[, , 1] <- a_tmp[1:nr, 1:nr]
  oirf[, , 1] <- a_tmp[1:nr, 1:nr] %*% p

  for (k in 1:len) {
    a_tmp <- a_tmp %*% big_a
    irf[, , (k + 1)] <- as.matrix(a_tmp[1:nr, 1:nr])
    oirf[, , (k + 1)] <- as.matrix(a_tmp[1:nr, 1:nr] %*% p)
  }

  ## TODO: add cumulative response functions
  out <- list()
  out$irf <- irf
  out$oirf <- oirf
  attr(out, "class") <- "irf"

  out
}
