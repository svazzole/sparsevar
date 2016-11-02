devtools::load_all()
sim <- simulateVAR(N = 30, p = 2, n = 300, sparsity = 0.15, rho = 0.8)

fit1 <- fitVAR(sim$series, p = 1)
fit2 <- fitVAR(sim$series, p = 2)
fit3 <- fitVAR(sim$series, p = 3)
fit4 <- fitVAR(sim$series, p = 4)
fit5 <- fitVAR(sim$series, p = 5)
informCrit(list(fit1,fit2,fit3,fit4,fit5))

informCrit(list(fit1,fit2,fit3))

