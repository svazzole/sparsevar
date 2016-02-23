## Sparse VAR (SVAR) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

Some R functions useful to estimate sparse VAR / VECM models.

The functions included are:
- `estimateVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+;
- `estimateVECM`: to estimate a sparse VECM (Vector Error Correction Model) using LS with penalty (again: ENET, SCAD or MC+);
- `simulateVAR`: to generate a sparse VAR multivariate time series;
- `mcSimulations`: to generate Monte Carlo simulations of SVAR and the relative estimation;
- `createSparseMatrix`: used to create sparse matrices with a given density;
- `plotMatrix`: useful to plot sparse matrices;

### Installation

To install:
```
install.packages("devtools")
devtools::install_github("svazzole/svar")
```
Check [here](https://www.rstudio.com/products/rpackages/devtools/) to understand which are the dependencies of `devtools` for your OS.


