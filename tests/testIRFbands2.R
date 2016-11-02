N <- 20
nMC <- 100

mse <- matrix(0, nrow = nMC, 2)

pb <- utils::txtProgressBar(min = 0, max = nMC, style = 3)

for(k in 1:nMC) {
  
  sim <- simulateVAR(N = 20, p = 1, sparsity = 0.1)
  est <- estimateVAR(sim$series, p = 1)
  
  irfSim <- impulseResponse(sim)
  irfEst <- impulseResponse(est)
  
  mseIRF <- matrix(0, nrow = N^2)
  mseOIRF <- matrix(0, nrow = N^2)
  
  for (i in 1:N) {
    for (j in 1:N) {
      mseIRF[(i-1)*N + j, ] <- mean((irfSim$irf[i,j,] - irfEst$irf[i,j,])^2)/mean((irfSim$irf[i,j,] + 1)^2)
      mseOIRF[(i-1)*N + j, ] <- mean((irfSim$oirf[i,j,] - irfEst$oirf[i,j,])^2)/mean((irfSim$oirf[i,j,]+1)^2)
    }
  }
  
  mse[k, 1] <- mean(mseIRF)
  mse[k, 2] <- mean(mseOIRF)
  
  utils::setTxtProgressBar(pb, k)
}

close(pb)

plot(mse[,1])
points(mse[,2], pch = 2, col = "red")

diff <- mse[,2] - mse[,1]
plot(diff)
mean(diff)

r <- mse[,2]/mse[,1]
mean(r)

