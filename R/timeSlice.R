timeSliceVAR <- function(data, p = 1, penalty = "ENET", opt) {

  if (penalty == "ENET") {
    # call timeslice with ENET
    out <- timeSliceVAR_ENET(data, p, opt)

  } else if (penalty == "SCAD" | penalty == "MCP" | penalty == "SCAD2") {
    # call timeslice with SCAD or MCP
    out <- timeSliceVAR_SCAD(data, p, opt, penalty)

  } else {
    # error
    stop("Unknown penalty. Possible values are \"ENET\", \"SCAD\" or \"MCP\".")

  }

  out$penalty <- penalty
  return(out)
}

timeSliceVAR_ENET <- function(data, p, opt) {

  t <- Sys.time()
  nr <- nrow(data)
  nc <- ncol(data)

  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  returnFit <- ifelse(!is.null(opt$returnFit), opt$returnFit, FALSE)
  methodCov <- ifelse(!is.null(opt$methodCov), opt$methodCov, "tiger")
  a  <- ifelse(!is.null(opt$alpha), opt$alpha, 1)
  l <- ifelse(!is.null(opt$leaveOut), opt$leaveOut, 10)
  ## TODO: Add the look ahead period > 1
  winLength <- nr - l
  horizon <- 1

  trDt <- transformData(data[1:winLength, ], p, opt)
  lam <- glmnet::glmnet(trDt$X, trDt$y, alpha = a)$lambda

  resTS <- matrix(0, ncol = l+1, nrow = length(lam))
  resTS[ , 1] <- lam

  for (i in 1:l) {
      d <- data[i:(winLength + i), ]
      fit <- varENET(d[1:(nrow(d)-1), ], p, lam, opt)
      resTS[, i+1] <- computeErrors(d, p, fit)
  }

  finalRes <- matrix(0, ncol = 3, nrow = length(lam))
  finalRes[,1] <- lam
  finalRes[,2] <- rowMeans(resTS[,2:(l+1)])
  for (k in 1:length(lam)) {
    finalRes[k,3] <- stats::sd(resTS[k,2:(l+1)])
  }

  ix <- which(finalRes[,2] == min(finalRes[,2]))[1]
  fit <- varENET(data, p, finalRes[ix, 1], opt)

  Avector <- stats::coef(fit, s = finalRes[ix, 1])
  A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)

  elapsed <- Sys.time() - t

  # If threshold = TRUE then set to zero all the entries that are smaller than
  # the threshold
  if (threshold == TRUE) {
    A <- applyThreshold(A, nr, nc, p)
  }
  
  # Get back the list of VAR matrices (of length p)
  A <- splitMatrix(A, p)

  # Now that we have the matrices compute the residuals
  res <- computeResiduals(trDt$series, A)

  # Create the output
  output = list()
  output$mu <- trDt$mu
  output$A <- A

  # Do you want the fit?
  if (returnFit == TRUE) {
    output$fit <- fit
  }

  output$lambda <- finalRes[ix, 1]
  output$mse <- finalRes[ix, 2]
  output$mseSD <- finalRes[ix, 3]
  output$time <- elapsed
  output$series <- trDt$series
  output$residuals <- res
  
  # Variance/Covariance estimation
  output$sigma <- estimateCovariance(res)
  
  output$penalty <- "ENET"
  output$method <- "timeSlice"
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"

  return(output)

  return(finalRes)
}

timeSliceVAR_SCAD <- function(data, p, opt, penalty) {

  t <- Sys.time()
  nr <- nrow(data)
  nc <- ncol(data)

  picasso <- ifelse(!is.null(opt$picasso), opt$picasso, FALSE)
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  returnFit <- ifelse(!is.null(opt$returnFit), opt$returnFit, FALSE)
  methodCov <- ifelse(!is.null(opt$methodCov), opt$methodCov, "tiger")  
  a  <- ifelse(!is.null(opt$alpha), opt$alpha, 1)
  ## TODO: Add the look ahead period > 1
  l <- ifelse(!is.null(opt$leaveOut), opt$leaveOut, 10)
  winLength <- nr - l
  horizon <- 1
  
  trDt <- transformData(data[1:winLength, ], p, opt)

  if (!picasso) {
    if (penalty == "SCAD") {
      lam <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y, family = "gaussian", penalty = "SCAD",
                            alpha = 1)$lambda
    } else if (penalty == "MCP") {
      lam <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y, family = "gaussian", penalty = "MCP",
                            alpha = 1)$lambda
    } else {

      lam <- sparsevar::scadReg(as(trDt$X, "dgCMatrix"), trDt$y, alpha = 1)$lambda
    }
  } else {
    lam <- picasso::picasso(trDt$X, trDt$y, method = "scad", nlambda = 100)$lambda
  }

  resTS <- matrix(0, ncol = l+1, nrow = length(lam))
  resTS[ , 1] <- lam

  for (i in 1:l) {

    d <- data[i:(winLength + i), ]
    if (!picasso) {
      if (penalty == "SCAD" | penalty == "SCAD2") {
        fit <- varSCAD(d[1:(nrow(d)-1), ], p, lam, opt, penalty)
        resTS[, i+1] <- computeErrors(d, p, fit, penalty = penalty)
      } else {
        fit <- varMCP(d[1:(nrow(d)-1), ], p, lam, opt)
        resTS[, i+1] <- computeErrors(d, p, fit, penalty = "MCP")
      }
    } else {
      trDt <- transformData(d, p, opt)
      fit <- picasso::picasso(trDt$X, trDt$y, method = "scad", lambda = lam)
      resTS[, i+1] <- computeErrorsPicasso(d, p, fit)
    }
  }

  finalRes <- matrix(0, ncol = 3, nrow = length(lam))
  finalRes[,1] <- lam
  finalRes[,2] <- rowMeans(resTS[,2:(l+1)])
  for (k in 1:length(lam)) {
    finalRes[k,3] <- stats::sd(resTS[k,2:(l+1)])
  }

  ix <- which(finalRes[,2] == min(finalRes[,2]))[1]

  if (!picasso) {
    if (penalty == "SCAD") {
      fit <- varSCAD(data, p, finalRes[ix,1], opt)
      Avector <- fit$beta[2:nrow(fit$beta), 1]
    } else if (penalty == "MCP") {
      fit <- varMCP(data, p, finalRes[ix,1], opt)
      Avector <- fit$beta[2:nrow(fit$beta), 1]
    } else {
      fit <- varSCAD(data, p, finalRes[ix,1], opt, penalty == "SCAD2")
      Avector <- fit$beta[1:nrow(fit$beta), 1]
    }
    A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
  } else {
    trDt <- transformData(data, p, opt)
    fit <- picasso::picasso(trDt$X, trDt$y, method = "scad", lambda = finalRes[ix,1])
    Avector <- fit$beta[, 1]
    A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
  }

  elapsed <- Sys.time() - t

  # If threshold = TRUE then set to zero all the entries that are smaller than
  # the threshold
  if (!is.null(opt$threshold)) {
    if (opt$threshold == TRUE) {
      tr <- 1 / sqrt(p*nc*log(nr))
      L <- abs(A) >= tr
      A <- A * L
    }
  }

  # Get back the list of VAR matrices (of length p)
  A <- splitMatrix(A, p)

  # Now that we have the matrices compute the residuals
  res <- computeResiduals(data, A)

  # Create the output
  output = list()
  output$mu <- trDt$mu
  output$A <- A

  # Do you want the fit?
  if (!is.null(opt$returnFit)) {
    if (opt$returnFit == TRUE) {
      output$fit <- fit
    }
  }

  output$lambda <- finalRes[ix, 1]
  output$mse <- finalRes[ix, 2]
  output$mseSD <- finalRes[ix, 3]
  output$time <- elapsed
  output$series <- trDt$series
  output$residuals <- res
  
  # Variance/Covariance estimation
  output$sigma <- estimateCovariance(res)
  
  output$penalty <- penalty
  output$method <- "timeSlice"
  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  return(output)

  return(finalRes)
}

computeErrors <- function(data, p, fit, penalty = "ENET") {

  nr <- nrow(data)
  nc <- ncol(data)
  l <- length(fit$lambda)

  err <- rep(0, ncol = 1, nrow = nr)

  for (i in 1:l) {

    if (penalty == "ENET") {
      Avector <- stats::coef(fit, s = fit$lambda[i])
      A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
    } else if (penalty == "SCAD" | penalty == "MCP") {
      Avector <- fit$beta[2:nrow(fit$beta), i]
      A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
    } else {
      Avector <- fit$beta[1:nrow(fit$beta), i]
      A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
    }

    A <- splitMatrix(A, p)

    n <- data[nr, ]

    f <- rep(0, nrow = nc, ncol = 1)
    tmpData <- data[((nr-1)- p + 1):(nr - 1), ]
    for (k in 1:p) {

      f <- f + A[[k]] %*% data[((nr-1) - (k-1)), ]

    }

    err[i] <- mean((f - n)^2)

  }

  return(err)
}

computeErrorsPicasso <- function(data, p, fit) {

  nr <- nrow(data)
  nc <- ncol(data)
  l <- length(fit$lambda)

  err <- rep(0, ncol = 1, nrow = nr)

  for (i in 1:l) {

    Avector <- fit$beta[, i]
    A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)

    A <- splitMatrix(A, p)

    n <- data[nr, ]

    f <- rep(0, nrow = nc, ncol = 1)
    tmpData <- data[((nr-1)- p + 1):(nr - 1), ]
    for (k in 1:p) {

      f <- f + A[[k]] %*% data[((nr-1) - (k-1)), ]

    }

    err[i] <- mean((f - n)^2)

  }

  return(err)
}
