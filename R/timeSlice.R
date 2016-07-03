#' @export
timeSliceVAR <- function(data, p = 1, penalty = "ENET", opt) {
  
  if (penalty == "ENET") {
    
    out <- timeSliceVAR_ENET(data, p, opt)
    
  } else if (penalty == "SCAD" | penalty == "MCP") {
    
    out <- timeSliceVAR_SCAD(data, p, opt, penalty)
  
  } else {
    
    stop("Unknown penalty. Possible values are \"ENET\", \"SCAD\" or \"MCP\".")
  
  }
  
  return(out)
}

timeSliceVAR_ENET <- function(data, p, opt) {
  
  t <- Sys.time()
  nr <- nrow(data)
  nc <- ncol(data)
  
  leaveOut <- ifelse(is.null(opt$leaveOut), 10, opt$leaveOut)
  winLength <- nr - leaveOut
  a  <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  
  horizon <- 1
  l <- 10
  
  trDt <- transformData(data[1:winLength, ], p, opt)
  lam <- glmnet::glmnet(trDt$X, trDt$y, alpha = a)$lambda
  
  resTS <- matrix(0, ncol = l+1, nrow = length(lam))
  resTS[ , 1] <- lam
  
  for (i in 1:l) {
    
      d <- data[i:(winLength + i), ]
      fit <- varENET2(d[1:(nrow(d)-1), ], p, lam, opt)
      resTS[, i+1] <- computeErrors(d, p, fit)
    
  }
  
  finalRes <- matrix(0, ncol = 2, nrow = length(lam))
  finalRes[,1] <- lam 
  finalRes[,2] <- rowMeans(resTS[,2:(l+1)])
  
  # TODO: check: one can have multiple mins. Choose the first one. Is this the 
  # better choice?
  ix <- which(finalRes[,2] == min(finalRes[,2]))[1]
  fit <- varENET2(data, p, finalRes[ix, 1], opt)
  
  Avector <- stats::coef(fit, s = finalRes[ix, 1])
  A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)

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
  output$time <- elapsed
  output$series <- trDt$series
  output$residuals <- res
  output$sigma <- cov(res)
  attr(output, "class") <- "var"
  attr(output, "type") <- "estimate"
  return(output)
  
  return(finalRes)
}

timeSliceVAR_SCAD <- function(data, p, opt, penalty) {
  
  t <- Sys.time()
  nr <- nrow(data)
  nc <- ncol(data)
  
  leaveOut <- ifelse(is.null(opt$leaveOut), 10, opt$leaveOut)
  winLength <- nr - leaveOut
  a  <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  
  horizon <- 1
  l <- 10
  
  trDt <- transformData(data[1:winLength, ], p, opt)
  
  if (penalty == "SCAD") {
    lam <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y, family = "gaussian", penalty = "SCAD",
                          alpha = 1)$lambda    
  } else {
    lam <- ncvreg::ncvreg(as.matrix(trDt$X), trDt$y, family = "gaussian", penalty = "MCP",
                          alpha = 1)$lambda
  }

  
  resTS <- matrix(0, ncol = l+1, nrow = length(lam))
  resTS[ , 1] <- lam
  
  for (i in 1:l) {
    
    d <- data[i:(winLength + i), ]
    if (penalty == "SCAD") {
      fit <- varSCAD2(d[1:(nrow(d)-1), ], p, lam, opt)
      resTS[, i+1] <- computeErrors(d, p, fit, penalty = "SCAD")
    } else {
      fit <- varMCP2(d[1:(nrow(d)-1), ], p, lam, opt)
      resTS[, i+1] <- computeErrors(d, p, fit, penalty = "MCP")
    }

  }
  
  finalRes <- matrix(0, ncol = 2, nrow = length(lam))
  finalRes[,1] <- lam 
  finalRes[,2] <- rowMeans(resTS[,2:(l+1)])
  
  # TODO: check: one can have multiple mins. Choose the first one. Is this the 
  # better choice?
  ix <- which(finalRes[,2] == min(finalRes[,2]))[1]
  
  if (penalty == "SCAD") {
    fit <- varSCAD2(data, p, finalRes[ix,1], opt)    
  } else {
    fit <- varMCP2(data, p, finalRes[ix,1], opt)
  }

  Avector <- fit$beta[2:nrow(fit$beta), 1]
  A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
  
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
  output$time <- elapsed
  output$series <- trDt$series
  output$residuals <- res
  output$sigma <- cov(res)
  attr(output, "class") <- "var"
  attr(output, "type") <- "estimate"
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