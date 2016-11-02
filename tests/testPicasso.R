################################################################
## Sparse linear regression
## Generate the design matrix and regression coefficient vector
n = 100
d = 400
X = matrix(rnorm(n*d), n, d)
beta = c(3,2,0,1.5,rep(0,d-4))
## Generate response using Gaussian noise, and fit sparse linear models
noise = rnorm(n)
Y = X%*%beta + noise
out.l1.cyclic = picasso(X, Y, nlambda=10)
out.l1.greedy = picasso(X, Y, nlambda=10, alg="greedy")
out.mcp.greedy = picasso(X, Y, nlambda=10, method="scad")
## Visualize the solution path
plot(out.l1.cyclic)
plot(out.l1.greedy)
plot(out.mcp.greedy)
