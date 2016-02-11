# Sparse VAR (SVAR)

Some R functions useful to estimate sparse VAR models.

The functions included are:
- `estimateVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+;
- `simulateVAR`: to generate a sparse VAR multivariate time series;
- `mcSimulations`: to generate Monte Carlo simulations of SVAR and the relative estimation;
- `createSparseMatrix`: used to create sparse matrices with a given density;
- `plotMatrix`: useful to plot sparse matrices;

## Installation

To install:
```
install.packages("devtools")
devtools::install_github("svazzole/svar")
```
Check [here](https://www.rstudio.com/products/rpackages/devtools/) to understand which are the dependencies of `devtools` for your OS.


## Estimation

Use `estimateVAR`. The arguments of the function are:
- `rets`: the multivariate time series (variables in columns, observations in rows);
- `p`: the order of the VAR model to be estimated;
- `penalty`: the penalty used in Least Squares (LS). Possible values are: `"ENET"`, `"SCAD"` or `"MCP"`;
- `options`: list of options. Some will depend on the penalty and some are global.

Global options:
- `parallel`: `TRUE` or `FALSE` (default). Parallel cross-validation (on the folds);
- `ncores`: if `parallel = TRUE` then you must specify the number of cores used for the parallelization (default = `1`).
- `nfolds`: 

`penalty = "ENET` options:
- `lambda`: `"lambda.min"` (default) or `"lambda.1se"`;
- `alpha`: a value in [0,1] (default `alpha = 1`). `alpha = 1` is LASSO regression, `alpha = 0` is Ridge LS;
- `type.measure`: `"mse"` (default) or `"mae"`;
- `nlambda`: number of lambdas used for cross validation.

`penalty = "SCAD"` or `"MCP"` options:
- `eps`: convergence tolerance

### Examples

```
results <- estimateVAR(rets)
```
will estimate VAR(1) process using LASSO regression on the dataset `rets`.


## Simulations

Use `simulateVAR`. The parameters for the function are:
- `N`: the dimension of the process;
- `nobs`: the number of observations of the process;
- `rho`: the variance/covariance "intensity";
- `sparsity`: the percentage of non zero elements in the matrix of the VAR;
- `method`: `"normal"` or `"bimodal"`.

```
sim <- simulateVAR(N = 100, nobs = 250, rho = 0.75, sparsity = 0.05, method = "normal")
```
