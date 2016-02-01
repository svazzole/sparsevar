library(glmnet)
library(hdf5)
library(Matrix)

hdf5load("data/est/2004.h5", load = TRUE, verbosity = 0, tidy = FALSE)

p <- prices$block0_values
r <- rets$block0_values

nAssets <- 100
permNo <- p[1, 1:nAssets]

rets <- r[3:nrow(r), 1:nAssets]

N <- ncol(rets)
d <- nrow(rets)

rm(prices, r, p)

for (i in 1:N){
  rets[, i] <- scale(rets[, i])
}
  
tmpX <- rets[1:(d-1), ]
tmpY <- rets[2:d, ] 

I <- diag(N)
# X <- Matrix(I %x% tmpX, nrow = d*N, ncol = N*N, sparse = TRUE)
X <- I %x% tmpX
Y <- as.vector(tmpY)
gc()

t <- Sys.time()
cvfit = cv.glmnet(X, Y, alpha = 1, nlambda = 100, type.measure="mse")
elapsed <- Sys.time() - t

plot(cvfit)

Av <- coef(cvfit, s = "lambda.min")

if(sum(Av[2:N^2]!=0) == 0){
  Av <- coef(cvfit, s = "lambda.1se")
}

A <- t(matrix(Av[2:nrow(Av)], nrow = N, ncol = N, byrow = TRUE))
image(A)


#####################
# TEST

source("simulations.R")
source("estimateVAR.R")
sim <- simulateVAR(N = 200)
genA <- sim$A
data <- sim$data$series

resSCAD <- estimateVAR(data, penalty = "SCAD")
ASCAD <- resSCAD$A


resLASSO <- estimateVAR(data)
ALASSO <- resLASSO$A

resLASSO2 <- estimateVAR(data, options = list(lambda = "lambda.1se"))
ALASSO2 <- resLASSO2$A

par(mfrow = c(2,2))
image(ASCAD[30:1, ])
image(genA[200:1, ])
image(ALASSO[200:1, ])
image(ALASSO2[30:1, ])
