## Sparse VAR (sparsevar)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![Version](https://img.shields.io/badge/version-0.1.0-oran.svg)](https://github.com/svazzole/sparsevar)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sparsevar)](https://cran.r-project.org/package=sparsevar)
[![Downloads](http://cranlogs.r-pkg.org/badges/sparsevar)](https://cran.r-project.org/package=sparsevar)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/sparsevar?color=brightgreen)](https://cran.r-project.org/package=sparsevar)
[![Build Status](https://travis-ci.org/svazzole/sparsevar.svg?branch=master)](https://travis-ci.org/svazzole/sparsevar)

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
Check [here](https://www.rstudio.com/products/rpackages/devtools/) to understand which are the dependencies of `devtools` for your OS.

### Quick start

To load the `sparsevar` package simply type
```r
library(sparsevar)
```

Using the function included in the package, we simply generate a 20x20 VAR(2) process
```r
set.seed(1)
sim <- simulateVAR(N = 20, p = 2)
```
This command will generate a model with two sparse matrices with 5% of non-zero entries and a Toeplitz variance-covariance matrix with rho = 0.5.
We can estimate the matrices of the process using for example
```r
fit <- fitVAR(sim$series, p = 2, threshold = TRUE)
```

The results can be seen by plotting the two `var` objects
```r
plotVAR(sim, fit)
```
the first row of the plot is made by the matrices of the simulated process and the second row is formed by their estimates.

The fit contains also the estimate of the variance/covariance matrix of the residuals
```r
plotMatrix(fit$sigma)
```

which can be compared with the covariance matrix of the errors of the generating process
```r
plotMatrix(sim$sigma)
```

### Usage

The functions included for model estimation are:

- `fitVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+;
- `fitVARX`: to estimate a sparse VAR-X model using ENET;
- `fitVECM`: to estimate a sparse VECM (Vector Error Correction Model) using LS with penalty (again: ENET, SCAD or MC+);
- `impulseResponse`: compute the impulse response function;
- `errorBands`: estimate the error bands for the IRF (using bootstrap).

For simulations:

- `simulateVAR`: to generate a sparse VAR multivariate time series;
- `simulateVARX`: to generate a sparse VARX time series;
- `createSparseMatrix`: used to create sparse matrices with a given density.

For plotting:

- `plotMatrix`: useful to plot matrices and sparse matrices;
- `plotVAR`: plot all the matrices of the model or models in input;
- `plotIRF`: plot IRF function;
- `plotGridIRF`: multiple plots of IRF.

### Papers using `sparsevar`
[[1](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005364)] Gibbons SM, Kearney SM, Smillie CS, Alm EJ (2017) Two dynamic regimes in the human gut microbiome. PLoS Comput Biol 13(2): e1005364.

[[2](https://doi.org/10.1016/j.insmatheco.2019.07.004)] Quentin Guibert, Olivier Lopez, Pierrick Piette, Forecasting mortality rate improvements with a high-dimensional VAR, Insurance: Mathematics and Economics, Volume 88, 2019, Pages 255-272, ISSN 0167-6687.

### References
[[1](http://projecteuclid.org/euclid.aos/1434546214)] Basu, Sumanta; Michailidis, George. Regularized estimation in sparse high-dimensional time series models. Ann. Statist. 43 (2015), no. 4, 1535--1567. doi:10.1214/15-AOS1315.

[[2](https://books.google.it/books/?id=COUFCAAAQBAJ&redir_esc=y)] Lütkepohl, Helmut. New Introduction to Multiple Time Series Analysis. Springer Science & Business Media, 2005, ISBN 3540277528.
