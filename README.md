# Sparse VAR (SVAR)

Some R functions useful to estimate sparse VAR models.

The functions included are:
- `simulateVAR`: to generate a sparse VAR multivariate time series
- `estimateVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+
- `mcSimulations`: to generate Monte Carlo simulations of SVAR and the relative estimation.

## Installation

To install:
```
install.packages("devtools")
devtools::install_github("svazzole/svar")
```
Check [here](https://www.rstudio.com/products/rpackages/devtools/) which are the dependencies of `devtools` for your OS.

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

## Estimation

Use `estimateVAR`. The parameters are:
- `rets`: the multivariate time series (variables in columns, observations in rows);
- `penalty`: the penalty used in Least Squares (LS). Possible values are: `"ENET"`, `"SCAD"` or `"MCP"`;
- `options`: list of options (...).
