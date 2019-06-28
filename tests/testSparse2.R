library(ncvreg)
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/sparsevar/src/scad.cpp")
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/ncvreg/src/maxprod.cpp")
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/ncvreg/src/standardize.cpp")
source(file = "/home/svazzole/workspace/r/sparsevar/R/scadReg.R")
n <- 200
p <- 23

X <- matrix(rnorm(n*p), n, p)
mult <- rep(1, ncol(X))
b <- rnorm(p)
y <- rnorm(n, X%*%b)
beta <- lm(y~X)$coef

# XX <- as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
scad <- scadReg(X,y,lambda=c(0),eps=.0001)
scad2 <- scadReg(X,y,lambda=c(0), eps=.0001)
scad3 <- scadReg(X,y, nlambda = 20, lambda.min = 0, eps = .0001)
ncv <- coef(ncvreg(X,y,nlambda = 20, penalty="SCAD",eps=.0001))
beta
scad
scad2
ncv

check(scad, beta[2:nrow(beta)], tolerance=.01, check.attributes=FALSE)
check(mcp, beta,tolerance=.01,check.attributes=FALSE)

scad3$beta[, 21] - beta[2:24]

abs((scad3$beta[,6]-scad3$beta[,5])/scad3$beta[,5])>.0001

abs((scad3_2$beta[,6]-scad3_2$beta[,5])/scad3_2$beta[,5])>.0001
