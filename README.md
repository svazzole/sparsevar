## Sparse VAR (sparsevar) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![Version](https://img.shields.io/badge/version-0.0.3-oran.svg)](https://github.com/svazzole/sparsevar)

Some R functions useful to estimate sparse VAR / VECM models.

### Installation

To install:
```
install.packages("devtools")
devtools::install_github("svazzole/sparsevar")
```
Check [here](https://www.rstudio.com/products/rpackages/devtools/) to understand which are the dependencies of `devtools` for your OS.

### Usage

The functions included are:
- `estimateVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+;
- `estimateVECM`: to estimate a sparse VECM (Vector Error Correction Model) using LS with penalty (again: ENET, SCAD or MC+);
- `simulateVAR`: to generate a sparse VAR multivariate time series;
- `mcSimulations`: to generate Monte Carlo simulations of SVAR and the relative estimation;
- `createSparseMatrix`: used to create sparse matrices with a given density;
- `plotMatrix`: useful to plot sparse matrices;

### References
[[1](http://projecteuclid.org/euclid.aos/1434546214)] Basu, Sumanta; Michailidis, George. Regularized estimation in sparse high-dimensional time series models. Ann. Statist. 43 (2015), no. 4, 1535--1567. doi:10.1214/15-AOS1315. 
