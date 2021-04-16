twoStepOLS <- function(series, p = 1, penalty = "ENET", method = "cv", ...) {

  ## TODO: rewrite this function and add p>1 support

  ## First step: estimate VAR using LASSO
  fit <- fitVAR(data = series, p = p, penalty = penalty, method = method, ...)

  N <- ncol(fit$A[[1]])
  nobs <- nrow(fit$series)

  bigA <- companionVAR(fit)

  trDt <- transformData(fit$series, p = p, opt = list(method = method, scale = FALSE, center = TRUE))

  nonZeroEntries <- as.matrix(bigA != 0)

  ## Create matrix R
  t <- as.vector(nonZeroEntries)
  n <- sum(t != 0)
  ix <- which(t != 0)
  j <- 1:n

  R <- matrix(0, ncol = n, nrow = length(t))
  for (k in 1:n) {
    R[ix[k], j[k]] <- 1
  }

  X <- as.matrix(trDt$X)
  y <- as.vector(t(fit$series[-(1:p), ]))

  # Metodo A MANO
  s <- corpcor::invcov.shrink(fit$residuals, verbose = FALSE)
  G <- t(fit$series[-nobs, ]) %*% fit$series[-nobs, ] / nobs

  V <- solve(t(R) %*% (kronecker(G, s) %*% R))
  VV <- nonZeroEntries
  VV[nonZeroEntries] <- diag(V)
  G1 <- solve(t(R) %*% (kronecker(t(fit$series[-nobs, ]) %*% fit$series[-nobs, ], s)) %*% R)
  G2 <- t(R) %*% (kronecker(t(fit$series[-nobs, ]), s))

  g <- G1 %*% G2 # [ , (N+1):(length(y) + N)]
  ga <- g %*% y

  b1 <- vector(length = N * N)
  b1 <- R %*% ga
  A <- matrix(b1, ncol = N, byrow = F)

  varCov <- R %*% (solve(t(R) %*% (kronecker(G, s)) %*% R) / nobs) %*% t(R)
  varA <- matrix(diag(varCov), ncol = N, byrow = F)

  result <- list()
  attr(result, "class") <- "var"

  result$A <- splitMatrix(A, p)
  result$varA <- list(varA)

  uA <- result$A[[1]] + 2 * sqrt(result$varA[[1]])
  lA <- result$A[[1]] - 2 * sqrt(result$varA[[1]])
  L <- (uA < 0) | (lA > 0)
  result$cleanA <- result$A[[1]] * L
  result$residuals <- fit$residuals
  result$varCov <- varCov
  return(result)
}
