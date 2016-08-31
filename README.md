## Sparse VAR (sparsevar) 
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) 
[![Version](https://img.shields.io/badge/version-0.0.8-oran.svg)](https://github.com/svazzole/sparsevar)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sparsevar)](https://cran.r-project.org/package=sparsevar)
[![Downloads](http://cranlogs.r-pkg.org/badges/sparsevar)](https://cran.r-project.org/package=sparsevar)

Some R functions useful to estimate sparse VAR / VECM models.

### Installation

To install the stable version from CRAN:
```{r}
install.package("sparsevar")
```

To install the developing version:
```{r}
install.packages("devtools")
devtools::install_github("svazzole/sparsevar")
```
Check [here](https://www.rstudio.com/products/rpackages/devtools/) to understand which are the dependencies of `devtools` for your OS.

### Quick start

To load the `sparsevar` package simply type
```{r}
library(sparsevar)
```

Using the function included in the package, we simply generate a 20x20 VAR(2) process
```{r}
set.seed(1)
sim <- simulateVAR(N = 20, p = 2)
```
This command will generate a model with two sparse matrices with 5% of non-zero entries and a Toeplitz variance-covariance matrix (with $\rho=0.5$).
We can estimate the matrices of the process using for example
```{r}
est <- fitVAR(sim$series, p = 2, threshold = TRUE)
```

The results can be seen by plotting the matrices
```{r}
plotComparisonVAR(sim, est)
```
the first row of the plot is made by the matrices of the simulated process and the second row is formed by their estimates.

One can also estimate the variance/covariance matrix of the residuals with 
```{r}
M <- cov(est$residuals)
plotMatrix(M)
```

and compare with the covariance matrix of the errors of the generating process
```{r}
plotMatrix(sim$S)
```

### Usage

The functions included for model estimation are:

- `fitVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+;
- `fitVECM`: to estimate a sparse VECM (Vector Error Correction Model) using LS with penalty (again: ENET, SCAD or MC+);
- `impulseResponse`: compute the impulse response function;
- `errorBands`: estimate the error bands for the IRF (using bootstrap);

For simulations:

- `simulateVAR`: to generate a sparse VAR multivariate time series;
- `createSparseMatrix`: used to create sparse matrices with a given density;

For plotting:

- `plotMatrix`: useful to plot sparse matrices;
- `plotVAR`: plot all the matrices of the model;
- `plotComparisonVAR`: plot the comparison between the matrices of the simulated model and the matrices of the estimate.
- `plotIRF`: plot IRF function;
- `plotGridIRF`: multiple plots of IRF.

### References
[[1](http://projecteuclid.org/euclid.aos/1434546214)] Basu, Sumanta; Michailidis, George. Regularized estimation in sparse high-dimensional time series models. Ann. Statist. 43 (2015), no. 4, 1535--1567. doi:10.1214/15-AOS1315. 
