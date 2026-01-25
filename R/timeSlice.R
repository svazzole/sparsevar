#' @export
time_slice_var <- function(data, p = 1, penalty = "ENET", opt) {
  if (penalty == "ENET") {
    # call timeslice with ENET
    out <- time_slice_var_enet(data, p, opt)
  } else if (penalty == "SCAD" || penalty == "MCP" || penalty == "SCAD2") {
    # call timeslice with SCAD or MCP
    out <- time_slice_var_scad(data, p, opt, penalty)
  } else {
    # error
    stop("Unknown penalty. Possible values are \"ENET\", \"SCAD\" or \"MCP\".")
  }

  out$penalty <- penalty
  out
}

time_slice_var_enet <- function(data, p, opt) {
  t <- Sys.time()
  nr <- nrow(data)
  nc <- ncol(data)

  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  return_fit <- ifelse(!is.null(opt$returnFit), opt$returnFit, FALSE)
  # Cov method: method_cov <- ifelse(!is.null(opt$methodCov),
  #                                  opt$methodCov, "tiger")
  a <- ifelse(!is.null(opt$alpha), opt$alpha, 1)
  l <- ifelse(!is.null(opt$leaveOut), opt$leaveOut, 10)
  ## TODO: Add the look ahead period > 1
  # Look ahead: horizon <- 1
  win_length <- nr - l

  tr_dt <- transform_data(data[1:win_length, ], p, opt)
  lam <- glmnet::glmnet(tr_dt$X, tr_dt$y, alpha = a)$lambda

  res_ts <- matrix(0, ncol = l + 1, nrow = length(lam))
  res_ts[, 1] <- lam

  for (i in 1:l) {
    d <- data[i:(win_length + i), ]
    fit <- var_enet(d[1:(nrow(d) - 1), ], p, lam, opt)
    res_ts[, i + 1] <- compute_errors(d, p, fit)
  }

  final_res <- matrix(0, ncol = 3, nrow = length(lam))
  final_res[, 1] <- lam
  final_res[, 2] <- rowMeans(res_ts[, 2:(l + 1)])
  for (k in seq_along(lam)) {
    final_res[k, 3] <- stats::sd(res_ts[k, 2:(l + 1)])
  }

  ix <- which(final_res[, 2] == min(final_res[, 2]))[1]
  fit <- var_enet(data, p, final_res[ix, 1], opt)

  a_vector <- stats::coef(fit, s = final_res[ix, 1])
  a <- matrix(a_vector[2:length(a_vector)], nrow = nc,
              ncol = nc * p, byrow = TRUE)

  elapsed <- Sys.time() - t

  # If threshold = TRUE then set to zero all the entries that are smaller than
  # the threshold
  if (threshold) {
    a <- apply_threshold(a, nr, nc, p)
  }

  # Get back the list of VAR matrices (of length p)
  a <- split_matrix(a, p)

  # Now that we have the matrices compute the residuals
  res <- compute_residuals(tr_dt$series, a)

  # Create the output
  output <- list()
  output$mu <- tr_dt$mu
  output$A <- a

  # Do you want the fit?
  if (return_fit == TRUE) {
    output$fit <- fit
  }

  output$lambda <- final_res[ix, 1]
  output$mse <- final_res[ix, 2]
  output$mseSD <- final_res[ix, 3]
  output$time <- elapsed
  output$series <- tr_dt$series
  output$residuals <- res

  # Variance/Covariance estimation
  output$sigma <- estimate_covariance(res)

  output$penalty <- "ENET"
  output$method <- "timeSlice"
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"

  output
}

time_slice_var_scad <- function(data, p, opt, penalty) {
  t <- Sys.time()
  nr <- nrow(data)
  nc <- ncol(data)

  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  a <- ifelse(!is.null(opt$alpha), opt$alpha, 1)
  ## TODO: Add the look ahead period > 1
  l <- ifelse(!is.null(opt$leaveOut), opt$leaveOut, 10)
  win_length <- nr - l

  tr_dt <- transform_data(data[1:win_length, ], p, opt)

  if (penalty == "SCAD") {
    lam <- ncvreg::ncvreg(as.matrix(tr_dt$X), tr_dt$y,
      family = "gaussian", penalty = "SCAD",
      alpha = 1
    )$lambda
  } else if (penalty == "MCP") {
    lam <- ncvreg::ncvreg(as.matrix(tr_dt$X), tr_dt$y,
      family = "gaussian", penalty = "MCP",
      alpha = 1
    )$lambda
  } else {
    stop("[WIP] Only SCAD and MCP regression are supported.")
  }

  res_ts <- matrix(0, ncol = l + 1, nrow = length(lam))
  res_ts[, 1] <- lam

  for (i in 1:l) {
    d <- data[i:(win_length + i), ]
    if (penalty == "SCAD" || penalty == "SCAD2") {
      fit <- var_scad(d[1:(nrow(d) - 1), ], p, lam, opt, penalty)
      res_ts[, i + 1] <- compute_errors(d, p, fit, penalty = penalty)
    } else {
      fit <- var_mcp(d[1:(nrow(d) - 1), ], p, lam, opt)
      res_ts[, i + 1] <- compute_errors(d, p, fit, penalty = "MCP")
    }
  }

  final_res <- matrix(0, ncol = 3, nrow = length(lam))
  final_res[, 1] <- lam
  final_res[, 2] <- rowMeans(res_ts[, 2:(l + 1)])
  for (k in seq_along(lam)) {
    final_res[k, 3] <- stats::sd(res_ts[k, 2:(l + 1)])
  }

  ix <- which(final_res[, 2] == min(final_res[, 2]))[1]

  if (penalty == "SCAD") {
    fit <- var_scad(data, p, final_res[ix, 1], opt)
    a_vector <- fit$beta[2:nrow(fit$beta), 1]
  } else if (penalty == "MCP") {
    fit <- var_mcp(data, p, final_res[ix, 1], opt)
    a_vector <- fit$beta[2:nrow(fit$beta), 1]
  } else {
    fit <- var_scad(data, p, final_res[ix, 1], opt, penalty == "SCAD2")
    a_vector <- fit$beta[seq_along(fit$beta), 1]
  }
  a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)

  elapsed <- Sys.time() - t

  # If threshold = TRUE then set to zero all the entries that are smaller than
  # the threshold
  if (threshold) {
    a <- apply_threshold(a, nr, nc, p)
  }

  # Get back the list of VAR matrices (of length p)
  a <- split_matrix(a, p)

  # Now that we have the matrices compute the residuals
  res <- compute_residuals(data, a)

  # Create the output
  output <- list()
  output$mu <- tr_dt$mu
  output$A <- a

  # Do you want the fit?
  if (!is.null(opt$returnFit)) {
    if (opt$returnFit == TRUE) {
      output$fit <- fit
    }
  }

  output$lambda <- final_res[ix, 1]
  output$mse <- final_res[ix, 2]
  output$mseSD <- final_res[ix, 3]
  output$time <- elapsed
  output$series <- tr_dt$series
  output$residuals <- res

  # Variance/Covariance estimation
  output$sigma <- estimate_covariance(res)

  output$penalty <- penalty
  output$method <- "timeSlice"
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  output
}

compute_errors <- function(data, p, fit, penalty = "ENET") {
  nr <- nrow(data)
  nc <- ncol(data)
  l <- length(fit$lambda)

  err <- rep(0, ncol = 1, nrow = nr)

  for (i in 1:l) {
    if (penalty == "ENET") {
      a_vector <- stats::coef(fit, s = fit$lambda[i])
      a <- matrix(a_vector[2:length(a_vector)],
                  nrow = nc, ncol = nc * p, byrow = TRUE)
    } else if (penalty == "SCAD" || penalty == "MCP") {
      a_vector <- fit$beta[2:nrow(fit$beta), i]
      a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)
    } else {
      a_vector <- fit$beta[seq_along(fit$beta), i]
      a <- matrix(a_vector, nrow = nc, ncol = nc * p, byrow = TRUE)
    }

    a <- split_matrix(a, p)
    n <- data[nr, ]
    f <- rep(0, nrow = nc, ncol = 1)
    for (k in 1:p) {
      f <- f + a[[k]] %*% data[((nr - 1) - (k - 1)), ]
    }
    err[i] <- mean((f - n)^2)
  }
  err
}
