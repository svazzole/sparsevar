#' @export
scadReg <- function(X, y, family="gaussian", penalty="SCAD",
                   gamma=3.7, alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100,
                   lambda, eps=.001, max.iter=1000, convex=TRUE, dfmax=p+1, penalty.factor=rep(1, ncol(X)),
                   warn=TRUE, returnX=FALSE, ...) {
  # Coersion
  # if (class(X) != "matrix") {
  #   tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
  #   if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  # }

  # Error checking
  standardize <- FALSE
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  if (family=="binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data")
  if (family=="binomial" & !identical(sort(unique(y)), 0:1)) y <- as.numeric(y==max(y))

  ## Deprication support
  dots <- list(...)
  if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda
  
  ## Set up XX, yy, lambda
  if (standardize) {
    std <- standardize2(X)
    XX <- std[[1]]
    center <- as.numeric(std[[2]])
    scale <- as.numeric(std[[3]])
    nz <- which(scale > 1e-6)
    if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
    penalty.factor <- penalty.factor[nz]
  } else {
    XX <- X
  }

  p <- ncol(XX)

  if (family=="gaussian") {
    yy <- y - mean(y)
  } else {
    yy <- y
  }
  n <- length(yy)
  if (missing(lambda)) {
    #lambda <- setupLambda(if (standardize) XX else X, yy, family, alpha, lambda.min, nlambda, penalty.factor)
    lambda <- setupLambda2(as.matrix(XX), yy, family, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  ## Fit
  if (family=="gaussian" & standardize == FALSE) {

    #res <- cdfit_gaussianTEST(XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    res <- cdfit_gaussianTEST(XX, yy, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    a <- rep(mean(y),nlambda)
    #b <- matrix(res[[1]], p, nlambda)
    b <- t(res[[1]])
    b[is.nan(b)] <- 0
    loss <- res[[2]]
    iter <- res[[3]]

  } else if (family=="gaussian" & standardize==FALSE & 1==0) {
    beta <- cdfit_rawTEST(X, y, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)))
    # b <- matrix(res[[1]], p, nlambda)
    # loss <- res[[2]]
    # iter <- res[[3]]
  #else if (family=="binomial") {
  #   res <- .Call("cdfit_binomial", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
  #   a <- res[[1]]
  #   b <- matrix(res[[2]], p, nlambda)
  #   loss <- res[[3]]
  #   iter <- res[[4]]
  # } else if (family=="poisson") {
  #   res <- .Call("cdfit_poisson", XX, yy, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor, alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
  #   a <- res[[1]]
  #   b <- matrix(res[[2]], p, nlambda)
  #   loss <- res[[3]]
  #   iter <- res[[4]]
  # }
  }
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)

  if (family!="gaussian" | standardize==TRUE) a <- a[ind]
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  #if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for some values of lambda")

  ## Local convexity?
  #convex.min <- if (convex & standardize) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), family, penalty.factor, a=a) else NULL

  ## Unstandardize
  if (standardize) {

    beta <- b/scale[nz]
    val <- structure(list(beta = beta,
                          lambda = lambda,
                          center = center,
                          scale = scale))
    return(val)
    # beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda))
    # bb <- b/scale[nz]
    # beta[nz+1,] <- bb
    # beta[1,] <- a - crossprod(center[nz], bb)
    
  } else {
    beta <- if (family=="gaussian") b else rbind(a, b)
  }
  # 
  # ## Names
  # varnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  # if (family!="gaussian" | standardize==TRUE) varnames <- c("(Intercept)", varnames)
  # dimnames(beta) <- list(varnames, round(lambda,digits=4))
  # 
  # ## Output
  # val <- structure(list(beta = beta,
  #                       iter = iter,
  #                       lambda = lambda,
  #                       penalty = penalty,
  #                       family = family,
  #                       gamma = gamma,
  #                       alpha = alpha,
  #                       convex.min = convex.min,
  #                       loss = loss,
  #                       penalty.factor = penalty.factor,
  #                       n = n),
  #                  class = "ncvreg")
  # if (family=="poisson") val$y <- y
  # if (returnX) {
  #   val$X <- XX
  #   val$center <- center
  #   val$scale <- scale
  #   val$y <- yy
  # }
  # val
}

setupLambda2 <- function(X, y, family, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)
  
  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- glm(y~X[, -ind], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }
  if (family=="gaussian") {
    zmax <- maxprod(X, fit$residuals, ind, penalty.factor) / n
  } else {
    zmax <- maxprod(X, residuals(fit, "working") * fit$weights, ind, penalty.factor) / n
  }
  lambda.max <- zmax/alpha
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  } else {
    lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  }
  
  if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
  lambda
}

convexMin <- function(b, X, penalty, gamma, l2, family, penalty.factor, a, Delta=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)
  
  if (penalty=="MCP") {
    k <- 1/gamma
  } else if (penalty=="SCAD") {
    k <- 1/(gamma-1)
  } else if (penalty=="lasso") {
    return(NULL)
  }
  if (l==0) return(NULL)
  
  val <- NULL
  for (i in 1:l) {
    A1 <- if (i==1) rep(1,p) else b[,i]==0
    if (i==l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[,i+1]==0
      U <- A1&A2
      L2 <- l2[i+1]
    }
    if (sum(!U)==0) next
    Xu <- X[,!U]
    p.. <- k*(penalty.factor[!U]!=0) - L2*penalty.factor[!U]
    if (family=="gaussian") {
      if (any(A1!=A2)) {
        eigen.min <- min(eigen(crossprod(Xu)/n - diag(p.., length(p..), length(p..)))$values)
      }
    } else if (family=="binomial") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      pi. <- exp(eta)/(1+exp(eta))
      w <- as.numeric(pi.*(1-pi.))
      w[eta > log(.9999/.0001)] <- .0001
      w[eta < log(.0001/.9999)] <- .0001
      Xu <- sqrt(w) * cbind(1,Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*p..)))$values)
    } else if (family=="poisson") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      mu <- exp(eta)
      w <- as.numeric(mu)
      Xu <- sqrt(w) * cbind(1,Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*p..)))$values)
    } else if (family=="cox") {
      eta <- if (i==l) X%*%b[,i] else X%*%b[,i+1]
      haz <- drop(exp(eta))
      rsk <- rev(cumsum(rev(haz)))
      h <- haz*cumsum(Delta/rsk)
      xwxn <- crossprod(sqrt(h) * Xu)/n
      eigen.min <- min(eigen(xwxn-diag(diag(xwxn)*p.., nrow(xwxn), ncol(xwxn)))$values)
    }
    
    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}

