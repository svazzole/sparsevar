## Sparse VAR (sparsevar)

[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![Version](https://img.shields.io/badge/version-1.0.0-oran.svg)](https://github.com/svazzole/sparsevar)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sparsevar)](https://cran.r-project.org/package=sparsevar)
[![Downloads](http://cranlogs.r-pkg.org/badges/sparsevar)](https://cran.r-project.org/package=sparsevar)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/sparsevar?color=brightgreen)](https://cran.r-project.org/package=sparsevar)

Some R functions useful to estimate sparse VAR / VECM models.

### Installation

To install the stable version from CRAN:

```r
install.package("sparsevar")
```

To install the developing version:

```r
install.packages("devtools")
devtools::install_github("svazzole/sparsevar", "master")
```

### Quick start

To load the `sparsevar` package simply type

```r
library(sparsevar)
```

Using the function included in the package, we simply generate a 20x20 VAR(2) process

```r
set.seed(1)
sim <- simulate_var(n = 20, p = 2)
```

This command will generate a model with two sparse matrices with 5% of non-zero entries and a Toeplitz variance-covariance matrix with rho = 0.5.
We can estimate the matrices of the process using for example

```r
fit <- fit_var(sim$series, p = 2, threshold = TRUE)
```

The results can be seen by plotting the two `var` objects

```r
plot_var(sim, fit)
```

the first row of the plot is made by the matrices of the simulated process and the second row is formed by their estimates.

The fit contains also the estimate of the variance/covariance matrix of the residuals

```r
plot_matrix(fit$sigma)
```

which can be compared with the covariance matrix of the errors of the generating process

```r
plot_matrix(sim$sigma)
```

### Usage

The functions included for model estimation are:

- `fit_var`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+;
- `fit_varx`: to estimate a sparse VAR-X model using ENET;
- `fit_vecm`: to estimate a sparse VECM (Vector Error Correction Model) using LS with penalty (again: ENET, SCAD or MC+);
- `impulse_response`: compute the impulse response function;
- `error_bands`: estimate the error bands for the IRF (using bootstrap).

For simulations:

- `simulate_var`: to generate a sparse VAR multivariate time series;
- `simulate_varx`: to generate a sparse VARX time series;
- `create_sparse_matrix`: used to create sparse matrices with a given density.

For plotting:

- `plot_matrix`: useful to plot matrices and sparse matrices;
- `plot_var`: plot all the matrices of the model or models in input;
- `plot_irf`: plot IRF function;
- `plot_grid_irf`: multiple plots of IRF.

### Papers using `sparsevar`

[[1](https://doi.org/10.1371/journal.pcbi.1005364)] Gibbons SM, Kearney SM, Smillie CS, Alm EJ (2017) Two dynamic regimes in the human gut microbiome. PLoS Comput Biol 13(2): e1005364.

[[2](https://doi.org/10.1016/j.insmatheco.2019.07.004)] Quentin Guibert, Olivier Lopez, Pierrick Piette. Forecasting mortality rate improvements with a high-dimensional VAR, Insurance: Mathematics and Economics, Volume 88, Pages 255-272, ISSN 0167-6687, 2019.

[[3](https://doi.org/10.3390/risks10110219)] Kularatne TD, Li J, Shi Y. Forecasting Mortality Rates with a Two-Step LASSO Based Vector Autoregressive Model. Risks. 2022; 10(11):219.

[[4](https://doi.org/10.1002/env.70041)] Schisa, V., and M. Farnè. 2025. “The Impact of Climatic Factors on Respiratory Pharmaceutical Demand: A Comparison of Forecasting Models for Greece.” Environmetrics36, no. 7: e70041.

[[5](https://doi.org/10.1007/978-3-031-85870-3_13)] Nason, G., Salnikov, D., Cortina-Borja, M. (2025). Modelling Clusters in Network Time Series with an Application to Presidential Elections in the USA. In: Trejos, J., Chadjipadelis, T., Grané, A., Villalobos, M. (eds) Data Science, Classification, and Artificial Intelligence for Modeling Decision Making. IFCS 2024. Studies in Classification, Data Analysis, and Knowledge Organization. Springer, Cham.

[[6](https://doi.org/10.1093/jrsssa/qnag048)] Tim J Boonen, Yuhuai Chen, Low-rank tensor autoregressive models for mortality modelling, Journal of the Royal Statistical Society Series A: Statistics in Society, 2026;, qnag048.

[[7](https://doi.org/10.1214/26-BA1589)]  Arkaprava Roy. Anindya Roy. Subhashis Ghosal. "Bayesian Learning of Relational Graphs in Semiparametric High-Dimensional Time Series." Bayesian Anal. Advance Publication 1 - 30, 2025.

[[8](https://doi.org/10.48550/arXiv.2410.22617)] Roy, A., Roy, A., & Ghosal, S. (2024). Relational Graph in Vector Autoregression: A Case Study on the Effect of the Great Recession on Connectivity of Economic Indicators. arXiv preprint arXiv:2410.22617.

[[9](https://doi.org/10.1038/s41598-021-03322-8)] Bitetto, A., Cerchiello, P. & Mertzanis, C. A data-driven approach to measuring epidemiological susceptibility risk around the world. Sci Rep 11, 24037 (2021).

### References

[[1](https://projecteuclid.org/euclid.aos/1434546214)] Basu, Sumanta; Michailidis, George. Regularized estimation in sparse high-dimensional time series models. Ann. Statist. 43 (2015), no. 4, 1535--1567. doi:10.1214/15-AOS1315.

[[2](https://books.google.it/books/?id=COUFCAAAQBAJ&redir_esc=y)] Lütkepohl, Helmut. New Introduction to Multiple Time Series Analysis. Springer Science & Business Media, 2005, ISBN 3540277528.
