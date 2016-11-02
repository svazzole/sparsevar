N <- 30
p <- 2
sim <- simulateVAR(N = N, p = p, sparsity = 0.1)
est <- fitVAR(sim$series, p = p)
plotComparisonVAR(sim,est)
plotMatrix(est$sigma)

irf <- impulseResponse(est)
irfSim <- impulseResponse(sim)

bb <- errorBandsIRF(est, irf, M = 1000, alpha = 0.01)

i <- 2
j <- 2
plotIRF(irf, bb, i, j, type = "irf")
plotIRF(irf, bb, i, j, type = "oirf")
plotIRF(irf, bb, i, j, type = "irf", bands = "sd")
plotIRF(irf, bb, i, j, type = "oirf", bands = "sd")

s <- sample(1:N, 3)
plotIRFGrid(irf, bb, s, type = "irf")
plotIRFGrid(irf, bb, s, type = "oirf")

####
# Error on IRF
m <- min(irf$irf[i,j,], irfSim$irf[i,j,])
M <- max(irf$irf[i,j,], irfSim$irf[i,j,])
plot(irf$irf[i,j,], type = "l", col = "black", ylim = c(m,M))
lines(irfSim$irf[i,j,], col = "blue")

####
# Error on OIRF
m <- min(irf$oirf[i,j,], irfSim$oirf[i,j,])
M <- max(irf$oirf[i,j,], irfSim$oirf[i,j,])
plot(irf$oirf[i,j,], type = "l", col = "black", ylim = c(m,M))
lines(irfSim$oirf[i,j,], col = "blue")

####
# IRF with BANDS
m <- min(bb$irfLB[i,j,], irf$irf[i,j,])
M <- max(bb$irfUB[i,j,], irf$irf[i,j,])
plot(bb$irfUB[i,j,], type = "l", lty = 2, col = "red", ylim = c(m,M))
lines(irfSim$irf[i,j,], col = "blue")
lines(irf$irf[i,j,], type = "l")
lines(bb$irfLB[i,j,], lty = 2, col = "red")

####
# OIRF with BANDS
m <- min(bb$oirfLB[i,j,], irf$oirf[i,j,])
M <- max(bb$oirfUB[i,j,], irf$oirf[i,j,])
plot(bb$oirfUB[i,j,], type = "l", lty = 2, col = "red", ylim = c(m,M))
lines(irfSim$oirf[i,j,], col = "blue")
lines(irf$oirf[i,j,], type = "l")
lines(bb$oirfLB[i,j,], type = "l", lty = 2, col = "red")
