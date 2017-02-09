#' @title Impulse Response Function
#' 
#' @description A function to estimate the Impulse Response Function of a given VAR.
#' 
#' @usage impulseResponse(v, len = 20)
#' 
#' @param v the data in the for of a VAR
#' @param len length of the impulse response function
#' 
#' @return \code{irf} a 3d array containing the impulse response function.
#' 
#' @export
impulseResponse <- function(v, len = 20) {
  
  # Check if v is a VAR object
  if (!checkIsVar(v)) { 
    stop("Input v must be a VAR object")
  }
  
  # Numerical problems in the estimated variance covariance
  e <- eigen(v$sigma)$values
  if (!is.null(e[e<=0])) {
    P <- t(chol(v$sigma, pivot = TRUE))
  } else {
    P <- t(chol(v$sigma))
  }
  
  bigA <- companionVAR(v)
  
  out <- getIRF(v, bigA, len = len, P)
  out$cholP <- P  # Add Choleski factorization to the output
  
  return(out)
}

#' @title Check Impulse Zero
#' 
#' @description A function to find which entries of the impulse response function 
#' are zero.
#' 
#' @usage checkImpulseZero(irf)
#' 
#' @param irf irf output from impulseResponse function
#' 
#' @return a matrix containing the indices of the impulse response function that
#' are 0.
#' 
#' @export
checkImpulseZero <- function(irf) {
  
  nx <- dim(irf)[1]
  ny <- dim(irf)[2]
  nz <- dim(irf)[3]
  logicalIrf <- matrix(0, nx, ny)
  
  for (z in 1:nz) {
    logicalIrf <- logicalIrf + abs(irf[ , , z])
  }
  
  logicalIrf <- logicalIrf == 0
  return(which(logicalIrf == TRUE, arr.ind = TRUE))
}

#' @title Error bands for IRF
#' 
#' @description A function to estimate the confidence intervals for irf and oirf.
#' 
#' @usage errorBandsIRF(v, irf, alpha, M, resampling, ...)
#'
#' @param v a var object as from fitVAR or simulateVAR
#' @param irf irf output from impulseResponse function
#' @param alpha level of confidence (default \code{alpha = 0.01})
#' @param M number of bootstrapped series (default \code{M = 100})
#' @param resampling type of resampling: \code{"bootstrap"} or \code{"jackknife"}
#' @param ... some options for the estimation: \code{verbose = TRUE} or \code{FALSE}, 
#' \code{mode = "fast"} or \code{"slow"}, \code{threshold = TRUE} or \code{FALSE}.
#' 
#' @return a matrix containing the indices of the impulse response function that
#' are 0.
#' 
#' @export
errorBandsIRF <- function(v, irf, alpha = 0.01, M = 100, resampling = "bootstrap", ...) {

  opt <- list(...)
  verbose <- ifelse(!is.null(opt$verbose), opt$verbose, TRUE)
  mode <- ifelse(!is.null(opt$mode), opt$mode, "fast")
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  thresholdType <- ifelse(!is.null(opt$thresholdType), opt$thresholdType, "soft")

  if (resampling == "bootstrap"){
    lambda <- v$lambda 
    p <- length(v$A)
    nr <- ncol(v$series) 
    nc <- ncol(v$A[[1]]) 
    len <- dim(irf$irf)[3]
    
    irfs <- array(data = rep(0,len*nc^2*M), dim = c(nc,nc,len+1, M))
    oirfs <- array(data = rep(0,len*nc^2*M), dim = c(nc,nc,len+1, M))
    
    if (verbose == TRUE){
      cat("Step 1 of 2: bootstrapping series and re-estimating VAR...\n")
      pb <- utils::txtProgressBar(min = 0, max = M, style = 3)
    }
    
    for (k in 1:M) {
      # create Xs and Ys (temp variables)
      o <- bootstrappedVAR(v)
      
      if (mode == "fast") {
        if (v$penalty == "ENET"){
          # fit ENET to a specific value of lambda
          fit <- varENET(o, p, lambda, opt = list(method = v$method, penalty = v$penalty))
          Avector <- stats::coef(fit, s = lambda)
          A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
        } else if (v$penalty == "SCAD") {
          fit <- varSCAD(o, p, lambda, opt = list(method = v$method, penalty = v$penalty))    
          Avector <- fit$beta[2:nrow(fit$beta), 1]
          A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
        } else {
          fit <- varMCP(o, p, lambda, opt = list(method = v$method, penalty = v$penalty))    
          Avector <- fit$beta[2:nrow(fit$beta), 1]
          A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
        }
        
        if (threshold == TRUE) {
          applyThreshold(A, nr, nc, p, type = thresholdType)
        }
        
        M <- cbind(diag(x = 1, nrow = (nc*(p-1)), ncol = (nc*(p-1))), matrix(0, nrow = (nc*(p-1)), ncol = nc))
        bigA <- rbind(A,M)  
        
      } else {
        # fit ENET on a series of lambdas 
        if (threshold == TRUE) {
          fit <- fitVAR(o, p, penalty=v$penalty, method=v$method, threshold = TRUE)          
        } else {
          fit <- fitVAR(o, p, penalty=v$penalty, method=v$method)
        }
        bigA <- companionVAR(fit)
      }

      tmpRes <- getIRF(v, bigA, len, irf$cholP)
      irfs[,,,k] <- tmpRes$irf
      oirfs[,,,k] <- tmpRes$oirf
      
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, k)
      }
    }
    
    if (verbose == TRUE){
      close(pb)
      cat("Step 2 of 2: computing quantiles...\n")
      pb <- utils::txtProgressBar(min = 0, max = (nc*nc), style = 3)
    }
    
    irfUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
    irfLB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
    oirfUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
    oirfLB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
    irfQUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
    irfQLB <- irfQUB
    oirfQUB <- irfQUB
    oirfQLB <- irfQUB
    
    a <- alpha/2
    qLB <- stats::qnorm(a)
    qUB <- stats::qnorm((1-a))
    
    for (i in 1:nc) {
      for (j in 1:nc) {
        for (k in 1:len) {
          
          irfQUB[i,j,k]  <- stats::quantile(irfs[i,j,k,], probs = (1-a), na.rm = TRUE)
          oirfQUB[i,j,k]<- stats::quantile(oirfs[i,j,k,], probs = (1-a), na.rm = TRUE)
          irfQLB[i,j,k] <- stats::quantile(irfs[i,j,k,], probs = a, na.rm = TRUE)
          oirfQLB[i,j,k] <- stats::quantile(oirfs[i,j,k,], probs = a, na.rm = TRUE)
          
          irfUB[i,j,k] <- qUB*stats::sd(irfs[i,j,k,])
          oirfUB[i,j,k] <- qUB*stats::sd(oirfs[i,j,k,])
          irfLB[i,j,k] <- qLB*stats::sd(irfs[i,j,k,])
          oirfLB[i,j,k] <- qLB*stats::sd(oirfs[i,j,k,])
          
          
        }
        if (verbose == TRUE){
          utils::setTxtProgressBar(pb, (i-1)*nc + j)
        }
      }
    }
    
    if (verbose == TRUE){
      close(pb)
    }
    
    output <- list()
    
    output$irfUB <- irfUB
    output$oirfUB <- oirfUB
    output$irfLB <- irfLB
    output$oirfLB <- oirfLB
    
    output$irfQUB <- irfQUB
    output$oirfQUB <- oirfQUB
    output$irfQLB <- irfQLB
    output$oirfQLB <- oirfQLB
    
    attr(output, "class") <- "irfBands"
    attr(output, "resampling") <- "bootstrap"
    return(output)
    
  }  else if (resampling == "jackknife") {
    
    output <- jackknife(v, irf, alpha = alpha, ...)
    return(output)
    
  } else if (resampling == "bootstrapOLS") {
    output <- bootstrapOLS(v, irf, alpha = alpha, ...)
  } else {
    stop("Unknown resampling method. Possible values are \"bootstrap\" or \"jackknife\"")
  }
  

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
  thresholdType <- ifelse(!is.null(opt$thresholdType), opt$thresholdType, "soft")
  
  irfs <- array(data = rep(0,len*nc^2*nr), dim = c(nc,nc,len+1, nr))
  oirfs <- array(data = rep(0,len*nc^2*nr), dim = c(nc,nc,len+1, nr))
    
  if (verbose == TRUE){
    cat("Step 1 of 2: jack knifing series and re-estimating VAR...\n")
    pb <- utils::txtProgressBar(min = 0, max = nr, style = 3)
  }
  
  for (k in 1:nr) {
    # create Xs and Ys (temp variables)
    data <- v$series[-k, ]
    trDt <- transformData(data, p, opt = list(method = v$method, penalty = v$penalty))
    trDt$X <- trDt$X
    trDt$y <- trDt$y

    # data <- v$series
    # trDt <- transformData(data, p, opt = list(method = v$method, penalty = v$penalty))
    # trDt$X <- trDt$X[-k, ]
    # trDt$y <- trDt$y[-k]
    
    if (mode == "fast") {
      if (v$penalty == "ENET"){
        # fit ENET to a specific value of lambda
        fit <- varENET(data, p, lambda, opt = list(method = v$method, penalty = v$penalty))
        Avector <- stats::coef(fit, s = lambda)
        A <- matrix(Avector[2:length(Avector)], nrow = nc, ncol = nc*p, byrow = TRUE)
      } else if (v$penalty == "SCAD") {
        fit <- varSCAD(data, p, lambda, opt = list(method = v$method, penalty = v$penalty))    
        Avector <- fit$beta[2:nrow(fit$beta), 1]
        A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
      } else {
        fit <- varMCP(data, p, lambda, opt = list(method = v$method, penalty = v$penalty))    
        Avector <- fit$beta[2:nrow(fit$beta), 1]
        A <- matrix(Avector, nrow = nc, ncol = nc*p, byrow = TRUE)
      }
      
      if (threshold == TRUE) {
        applyThreshold(A, nr, nc, p, type = thresholdType)
      }
      
      M <- cbind(diag(x = 1, nrow = (nc*(p-1)), ncol = (nc*(p-1))), matrix(0, nrow = (nc*(p-1)), ncol = nc))
      bigA <- rbind(A,M)  
      
    } else {
      # fit ENET on a series of lambdas 
      if (threshold == TRUE) {
        fit <- fitVAR(data, p, penalty=v$penalty, method=v$method, threshold = TRUE)        
      } else {
        fit <- fitVAR(data, p, penalty=v$penalty, method=v$method)
      }
      bigA <- companionVAR(fit)
    }

    tmpRes <- getIRF(v, bigA, len, irf$cholP)
    irfs[,,,k] <- tmpRes$irf
    oirfs[,,,k] <- tmpRes$oirf
    
    if (verbose == TRUE){
      utils::setTxtProgressBar(pb, k)
    }
  }
  
  if (verbose == TRUE){
    close(pb)
    cat("Step 2 of 2: computing quantiles...\n")
    pb <- utils::txtProgressBar(min = 0, max = (nc*nc), style = 3)
  }
  
  irfUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  irfLB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  oirfUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  oirfLB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  
  a <- alpha/2
  qLB <- stats::qnorm(a)
  qUB <- stats::qnorm((1-a))
  
  for (i in 1:nc) {
    for (j in 1:nc) {
      for (k in 1:len) {

        irfUB[i,j,k] <- base::mean(irfs[i,j,k,]) + qUB*stats::sd(irfs[i,j,k,])
        oirfUB[i,j,k] <- base::mean(oirfs[i,j,k,]) + qUB*stats::sd(oirfs[i,j,k,])
        irfLB[i,j,k] <- base::mean(irfs[i,j,k,]) + qLB*stats::sd(irfs[i,j,k,])
        oirfLB[i,j,k] <- base::mean(oirfs[i,j,k,]) + qLB*stats::sd(oirfs[i,j,k,])
        
      }
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, (i-1)*nc + j)
      }
    }
  }
  
  if (verbose == TRUE){
    close(pb)
  }
  
  output <- list()
  
  output$irfUB <- irfUB
  output$oirfUB <- oirfUB
  output$irfLB <- irfLB
  output$oirfLB <- oirfLB
  
  attr(output, "class") <- "irfBands"
  attr(output, "resampling") <- "jackknife"
  return(output)
}


bootstrapOLS <- function(v, irf, mode = "fast", alpha, ...) {
  
  lambda <- v$lambda 
  p <- length(v$A)
  nc <- ncol(v$A[[1]]) 
  len <- dim(irf$irf)[3]
  nr <- nrow(v$series)
  
  opt <- list(...)
  verbose <- ifelse(!is.null(opt$verbose), opt$verbose, TRUE)
  mode <- ifelse(!is.null(opt$mode), opt$mode, "fast")
  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  thresholdType <- ifelse(!is.null(opt$thresholdType), opt$thresholdType, "soft")
  
  irfs <- array(data = rep(0,len*nc^2*nr), dim = c(nc,nc,len+1, nr))
  oirfs <- array(data = rep(0,len*nc^2*nr), dim = c(nc,nc,len+1, nr))
  
  if (verbose == TRUE){
    cat("Step 1 of 2: bootstrappingOLS series and re-estimating VAR...\n")
    pb <- utils::txtProgressBar(min = 0, max = nr, style = 3)
  }
  
  for (k in 1:nr) {
    # create Xs and Ys (temp variables)
    o <- bootstrappedVAR(v)
    
    N <- ncol(v$A[[1]])
    nobs <- nrow(v$series)
    
    bigA <- companionVAR(v)
    
    trDt <- transformData(o, p = p, opt = list(method = v$method, scale = FALSE, center = TRUE))
    
    nonZeroEntries <- as.matrix(bigA != 0)
    
    ## Create matrix R
    t <- as.vector(nonZeroEntries)
    n <- sum(t != 0)
    ix <- which(t != 0)
    j <- 1:n
    
    R <- matrix(0, ncol = n, nrow = length(t))
    for(zz in 1:n) {
      R[ix[zz],j[zz]] <- 1
    } 
    
    X <- as.matrix(trDt$X)
    y <- as.vector(t(o[-(1:p), ]))
    
    # Metodo A MANO
    s <- corpcor::invcov.shrink(v$residuals, verbose = FALSE)
    G <- t(o[-nobs, ])%*%o[-nobs, ] / nobs
    
    V <- solve(t(R)%*%(kronecker(G,s)%*%R))
    VV <- nonZeroEntries
    VV[nonZeroEntries] <- diag(V)
    G1 <- solve(t(R)%*%(kronecker(t(o[-nobs, ])%*%o[-nobs, ], s))%*%R) 
    G2 <- t(R)%*%(kronecker(t(o[-nobs,]), s))
    
    g <- G1 %*% G2#[ , (N+1):(length(y) + N)]
    ga <- g %*% y
    
    b1 <- vector(length = N*N)
    b1 <- R%*%ga
    A <- matrix(b1, ncol = N, byrow = F)

      
    M <- cbind(diag(x = 1, nrow = (nc*(p-1)), ncol = (nc*(p-1))), matrix(0, nrow = (nc*(p-1)), ncol = nc))
    bigA <- rbind(A,M)  
      
    tmpRes <- getIRF(v, bigA, len, irf$cholP)
    irfs[,,,k] <- tmpRes$irf
    oirfs[,,,k] <- tmpRes$oirf
    
    if (verbose == TRUE){
      utils::setTxtProgressBar(pb, k)
    }
  }
  
  if (verbose == TRUE){
    close(pb)
    cat("Step 2 of 2: computing quantiles...\n")
    pb <- utils::txtProgressBar(min = 0, max = (nc*nc), style = 3)
  }
  
  irfUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  irfLB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  oirfUB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  oirfLB <- array(data = rep(0,len*nc^2), dim = c(nc,nc,len))
  
  a <- alpha/2
  qLB <- stats::qnorm(a)
  qUB <- stats::qnorm((1-a))
  
  for (i in 1:nc) {
    for (j in 1:nc) {
      for (k in 1:len) {
        
        irfUB[i,j,k] <- base::mean(irfs[i,j,k,]) + qUB*stats::sd(irfs[i,j,k,])
        oirfUB[i,j,k] <- base::mean(oirfs[i,j,k,]) + qUB*stats::sd(oirfs[i,j,k,])
        irfLB[i,j,k] <- base::mean(irfs[i,j,k,]) + qLB*stats::sd(irfs[i,j,k,])
        oirfLB[i,j,k] <- base::mean(oirfs[i,j,k,]) + qLB*stats::sd(oirfs[i,j,k,])
        
      }
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, (i-1)*nc + j)
      }
    }
  }
  
  if (verbose == TRUE){
    close(pb)
  }
  
  output <- list()
  
  output$irfUB <- irfUB
  output$oirfUB <- oirfUB
  output$irfLB <- irfLB
  output$oirfLB <- oirfLB
  
  attr(output, "class") <- "irfBands"
  attr(output, "resampling") <- "jackknife"
  return(output)
}


getIRF <- function (v, bigA, len = 20, P) {
  
  nr <- nrow(v$A[[1]])
  
  irf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len+1))
  oirf <- array(data = rep(0,len*nr^2), dim = c(nr,nr,len+1))
  
  Atmp <- diag(nrow = nrow(bigA), ncol = ncol(bigA))
  
  irf[ , , 1] <- Atmp[1:nr, 1:nr] 
  oirf[ , , 1] <- Atmp[1:nr, 1:nr] %*% P
  
  for (k in 1:len) {
    Atmp <- Atmp %*% bigA
    irf[ , , (k+1)] <- as.matrix(Atmp[1:nr, 1:nr])
    oirf[ , , (k+1)] <- as.matrix(Atmp[1:nr, 1:nr] %*% P)
  }
  
  ## TODO: add cumulative response functions
  out <- list()
  out$irf <- irf
  out$oirf <- oirf
  attr(out, "class") <- "irf"
  
  return(out)
}
