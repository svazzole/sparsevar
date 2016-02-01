library(glmnet)
library(Matrix)
library(corrplot)

source("simulations.R")

N <- 50
s <- simulateVAR(nobs = 200, N = N, rho = 0.5, sparsity = 0.15)
rets <- s$data$series
genA <- s$A

#################################
N <- ncol(rets)
d <- nrow(rets)

for (i in 1:N){
  rets[, i] <- scale(rets[, i])
}

tmpX <- rets[1:(d-1), ]
tmpY <- rets[2:d, ] 

I <- diag(N)
X <- Matrix(I %x% tmpX, sparse = TRUE)
Y <- as.vector(tmpY)
gc()

t <- Sys.time()
cvfit = cv.glmnet(X, Y, alpha = 0.95, nlambda = 100, type.measure="mse")
elapsed <- Sys.time() - t

Av <- coef(cvfit, s = "lambda.min")
A <- matrix(Av[2:nrow(Av)], nrow = N, ncol = N, byrow = TRUE)
Av2 <- coef(cvfit, s = "lambda.1se")
A2 <- matrix(Av2[2:nrow(Av)], nrow = N, ncol = N, byrow = TRUE)
######################
L <- A
L[L!=0] <- 1
L[L==0] <- 0

L2 <- A2
L2[L2!=0] <- 1
L2[L2==0] <- 0

genL <- genA
genL[genL!=0] <- 1
genL[genL==0] <- 0

cat("Accuracy L:", 1 - sum(abs(L-genL))/N^2, sep = " ")
cat("Accuracy L2:", 1 - sum(abs(L2-genL))/N^2, sep = " ")
cat("Sparsities:\n", "L: ", sum(L)/N^2, "\n L2: ", sum(L2)/N^2, "\n genL: ", sum(genL)/N^2)

A <- atan(A) * (2/pi)
A2 <- atan(A2) * (2/pi)
genA <- atan(genA) * (2/pi)

par(mfrow = c(2,2))
image(A, col = rainbow(12))
#corrplot(A, method = "square")
image(genA, col = rainbow(12))
#corrplot(genA, method = "square")
plot(cvfit)
image(A2, col = rainbow(12))
#corrplot(A2, method = "square")
