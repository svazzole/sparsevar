library(glmnet)
library(Matrix)

mcSimulations <- function(N, nobs = 250, nMC) {
  
  source("simulations.R")
  source("utils.R")
  
  rho <- 0.9
  sparsity <- 0.05
    
  results <- matrix(0, nMC, 8)
  
  pb <- txtProgressBar(min = 0, max = nMC, style = 3)
    
  for (i in 1:nMC){

      s <- simulateVAR(nobs = nobs, N = N, rho = rho, sparsity = sparsity)
      rets <- s$data$series
      genA <- s$A
      
      N <- ncol(rets)
      d <- nrow(rets)
      
      for (j in 1:N){
        rets[, j] <- scale(rets[, j])
      }
      
      tmpX <- rets[1:(d-1), ]
      tmpY <- rets[2:d, ] 
      
      I <- diag(N)
      #X <- Matrix(I %x% tmpX, nrow = d*N, ncol = N*N, sparse = TRUE)
      X <- Matrix(I %x% tmpX, sparse = TRUE)
      Y <- as.vector(tmpY)
      gc()
      
      t <- Sys.time()
      cvfit = cv.glmnet(X, Y, alpha = 1, nlambda = 100, type.measure="mse")
      elapsed <- Sys.time() - t
      
      Av <- coef(cvfit, s = "lambda.min")
      A <- matrix(Av[2:nrow(Av)], nrow = N, ncol = N, byrow = TRUE)
      Av2 <- coef(cvfit, s = "lambda.1se")
      A2 <- matrix(Av2[2:nrow(Av)], nrow = N, ncol = N, byrow = TRUE)
      
      L <- A
      L[L!=0] <- 1
      L[L==0] <- 0
      
      L2 <- A2
      L2[L2!=0] <- 1
      L2[L2==0] <- 0
      
      genL <- genA
      genL[genL!=0] <- 1
      genL[genL==0] <- 0
      
      results[i, 1] <- 1 - sum(abs(L-genL))/N^2
      results[i, 2] <- 1 - sum(abs(L2-genL))/N^2
      results[i, 3] <- sum(L)/N^2
      results[i, 4] <- sum(L2)/N^2
      results[i, 5] <- l2norm(L-genL)
      results[i, 6] <- l2norm(L2-genL)
      results[i, 7] <- frobNorm(L-genL)
      results[i, 8] <- frobNorm(L2-genL)
      
      setTxtProgressBar(pb, i)
    }
  
  close(pb)
  
  results <- as.data.frame(results)
  colnames(results) <- c("Accuracy_A", "Accuracy_A2", "Sparsity_A", "Sparsity_A2", "l2_A", "l2_A2", "frob_A", "frob_A2")
  
  return(results)
  
}
