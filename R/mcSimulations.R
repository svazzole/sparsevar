#' @title Monte Carlo simulations
#'
#' @description This function generates monte carlo simultaions of SVAR and its
#' estimation.
#' @param N dimension of the multivariate time series.
#' @param nobs number of observations to be generated.
#' @param nMC number of Monte Carlo simulations.
#' @param rho base value for the covariance.
#' @param sparsity density of non zero entries of the VAR matrices.
#' @param penalty penalty function to use for LS estimation. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}.
#' @param covariance type of covariance matrix to be used in the generation of the SVAR model.
#' @param options (TODO: complete)
#' 
#' @return A list containing ...
#' 
#' @author Simone Vazzoler
#'
#' @export
#' 
mcSimulations <- function(N, nobs = 250, nMC = 100, rho = 0.5, sparsity = 0.05, 
                          penalty = "ENET", covariance = "toeplitz", 
                          options = NULL, method = "normal") {

  results <- matrix(0, nMC, 5)
  
  pb <- txtProgressBar(min = 0, max = nMC, style = 3)
    
  for (i in 1:nMC){

      s <- simulateVAR(nobs = nobs, N = N, rho = rho, sparsity = sparsity, covariance = covariance, method = method)
      rets <- s$data$series
      genA <- s$A[[1]]
      
      res <- estimateVAR(rets, penalty = penalty, options = options)
      
      A <- res$A[[1]]
      
      L <- A
      L[L!=0] <- 1
      L[L==0] <- 0
      
      genL <- genA
      genL[genL!=0] <- 1
      genL[genL==0] <- 0
      
      results[i, 1] <- 1 - sum(abs(L-genL))/N^2   # accuracy    -(1 - sum(genL)/N^2)
      results[i, 2] <- abs(sum(L)/N^2 - sparsity) # sparsity
      results[i, 3] <- l2norm(A-genA) / l2norm(genA)
      results[i, 4] <- frobNorm(A-genA) / frobNorm(genA)
      results[i, 5] <- res$mse
      setTxtProgressBar(pb, i)
    }
  
  close(pb)
  
  results <- as.data.frame(results)
  colnames(results) <- c("accuracy", "sparDiff", "l2", "frob", "mse")
  
  return(results)
  
}
