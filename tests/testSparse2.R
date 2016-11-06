library(ncvreg)
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/sparsevar/src/scad.cpp")
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/ncvreg/src/maxprod.cpp")
Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/ncvreg/src/standardize.cpp")
source(file = "/home/svazzole/workspace/r/sparsevar/R/scadReg.R")
n <- 200
p <- 20

X <- matrix(rnorm(n*p), n, p)
mult <- rep(1, ncol(X))
b <- rnorm(p)
y <- rnorm(n, X%*%b)
beta <- lm(y~X)$coef

XX <- as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
scad <- scadReg(X,y)
scad2 <- scadReg(XX,y,lambda=c(0), eps=.0001)
scad3 <- scadReg(X,y, nlambda = 10)
ncv <- coef(ncvreg(X,y,lambda=0,penalty="SCAD",eps=.0001))
beta
scad
scad2
ncv

check(scad, beta[2:nrow(beta)], tolerance=.01, check.attributes=FALSE)
check(mcp, beta,tolerance=.01,check.attributes=FALSE)

