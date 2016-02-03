
source("simulations.R")
source("mcSimulations.R")
source("estimateVAR.R")
source("utils.R")

direct <- "results/"

#########
# LASSO #
#########

N <- 10
rho <- 0.5
s <- 0.05
nobs <- 30

sim <- simulateVAR(N = N, nobs = nobs, rho = rho, sparsity = s)
genA <- sim$A
rets <- sim$data$series

optLambda1se <- list(parallel = TRUE, ncores = 2, lambda = "lambda.1se")
optPar <- list(parallel = TRUE, ncores = 2)

res <- estimateVAR(rets, options = optPar)

A <- res$A

L <- A
L[L!=0] <- 1
L[L==0] <- 0

genL <- genA
genL[genL!=0] <- 1
genL[genL==0] <- 0

cat("Relative Accuracy: ") 
1 - sum(abs(L-genL))/N^2 - (1 - sum(genL)/N^2)  # accuracy    -(1 - sum(genL)/N^2)
cat("Relative Sparsity: ")
abs(sum(L)/N^2 - sum(genL)/N^2)                 # sparsity
cat("Relative l2: ")
l2norm(A-genA) / l2norm(genA)
cat("Relative Froebenius: ")
frobNorm(A-genA) / frobNorm(genA)

par(mfrow = c(1,2))
image(genA[N:1, ])
image(A[N:1, ])

############

results <- mcSimulations(N = N, rho = rho, options = optLambda1se)
res <- mcSimulations(N = N, rho = rho, options = optPar)
colMeans(results)
colMeans(res)