## Test1 -----------------------------------------------------------------------
suppressMessages(library(Matrix))
i <- sample(1:100, 20)
j <- sample(1:100, 20)
x <- rnorm(20)
a <- sparseMatrix(dims = c(100,100), i, j, x = x)
print(a)
b <- rep(1, 100)
crossprod(a, b)

## Test2 -----------------------------------------------------------------------
sim <- simulate_var(N = 100)

tr_dt <- transformData(sim$series, 1, list())

x <- as(tr_dt$X, "dgCMatrix")
y <- tr_dt$y

Rcpp::sourceCpp(file = "/home/svazzole/workspace/r/sparsevar/src/scad.cpp")
z <- sparsevar::crossprod(x, y)
z <- crossprod(x, y)
