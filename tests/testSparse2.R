library(ncvreg)
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/sparsevar/src/scad.cpp")
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/ncvreg/src/maxprod.cpp")
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/ncvreg/src/standardize.cpp")
source(file = "/home/svazzole/workspace/r/sparsevar/R/scadReg.R")
n <- 200
p <- 23

x <- matrix(rnorm(n * p), n, p)
mult <- rep(1, ncol(x))
b <- rnorm(p)
y <- rnorm(n, x %*% b)
beta <- lm(y ~ x)$coef

scad <- scadReg(x, y, lambda = c(0), eps = .0001)
scad2 <- scadReg(x, y, lambda = c(0), eps = .0001)
scad3 <- scadReg(x, y, nlambda = 20, lambda.min = 0, eps = .0001)
ncv <- coef(ncvreg(x, y, nlambda = 20, penalty = "SCAD", eps = .0001))
beta
scad
scad2
ncv

check(scad, beta[2:nrow(beta)], tolerance = .01, check.attributes = FALSE)
check(mcp, beta, tolerance = .01, check.attributes = FALSE)

scad3$beta[, 21] - beta[2:24]

abs((scad3$beta[, 6] - scad3$beta[, 5]) / scad3$beta[, 5]) > .0001
abs((scad3_2$beta[, 6] - scad3_2$beta[, 5]) / scad3_2$beta[, 5]) > .0001
