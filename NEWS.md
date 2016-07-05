# sparsevar 0.0.7

- Major code rewriting
- Remove dependecies from MTS and caret
- Added impulse response error bands (using bootstrap)
- Added plot functions for IRF
- New timeSlice estimation
- Removed repeated cross validation

# sparsevar 0.0.6

- Added impulse response function for VAR processes

# sparsevar 0.0.5

- Added timeSlice estimation
- Fixed normalization constant in creating sparse var matrix
- Fixed parallel backend in Windows 

# sparsevar 0.0.4

- Added as output the residuals of the estimation (for estimateVAR)
- Fixed parallel background in LASSO estimation
- Now repeatedCV returns MSE
