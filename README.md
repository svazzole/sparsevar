# Sparse VAR (SVAR)

Some R functions useful to estimate sparse VAR models.

For the moment only VAR(1) is implemented.

The functions included are:
- `simulateVAR`: to generate a sparse VAR multivariate time series
- `estimateVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+
- `mcSimulations`: to generate Monte Carlo simulations of SVAR and the relative estimation.

## Simulations

Use `simulateVAR`. The parameters for the function are:
- `N`: the dimension of the process;
- `nobs`: the number of observations of the process;
- `rho`: the variance/covariance "intensity";
- `sparsity`: the percentage of non zero elements in the matrix of the VAR;
- `method`: `"normal"` or `"bimodal"`.

```
sim <- simulateVAR(N = 100,nobs = 250, rho = 0.75, sparsity = 0.05, method = "normal")
```

