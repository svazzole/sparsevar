# Sparse VAR (SVAR)

Some R functions useful to estimate sparse VAR models.

For the moment only VAR(1) is implemented.

The functions included are:
- `simulateVAR`: to generate a sparse VAR multivariate time series
- `estimateVAR`: to estimate a sparse VAR multivariate time series with ENET, SCAD or MC+
- `mcSimulations`: to generate Monte Carlo simulations of SVAR and the relative estimation.

### `simulateVAR`

Use:

