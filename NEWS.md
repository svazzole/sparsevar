# sparsevar 1.1.0 (in development)

## Breaking changes

- `fit_var()`, `fit_varx()`, `error_bands_irf()` and `jackknife()` now reject
  unrecognised names passed through `...`, suggesting the nearest known option.
  Previously any misspelled option was silently ignored, so
  `fit_var(d, p = 1, thersholdd_type = "hard")` ran and quietly returned the
  default answer. Code that passes an option this package does not recognise
  will now fail loudly; that code was already not doing what it looked like.
- The legacy `thresholdType` spelling is no longer accepted anywhere. It has
  had no effect in `fit_var()` since version 1.0.0 and was never documented in
  any release; it is now an error suggesting `threshold_type`.
- `horizon` is rejected unless it is 1. It has been documented since 0.0.7 but
  never read — the time-slice loop only ever looked one step ahead.
- `threshold_type` now defaults to `"theory"` instead of `"soft"`. This changes
  results for anyone setting `threshold = TRUE` without naming a type: `"soft"`
  is close to a no-op, evaluating to 0.0018 at the largest setting in the
  package benchmark. A note naming the new default is shown once per session
  when the default is used; it will be removed in 1.2.0. Set
  `threshold_type = "soft"` for the previous behaviour.

## New features

- `threshold_type = "theory"`, a per-coefficient cutoff
  `tau_ij = (sigma_i / s_j) * C * sqrt(log(n^2 p) / (T - p))`, where `sigma_i`
  is the residual standard deviation of equation `i` from the un-thresholded
  fit and `s_j` the standard deviation of series `j`. Unlike the two legacy
  cutoffs it is derived from the estimation theory and carries the units of
  the coefficients, so it is not distorted when series are measured on
  different scales. It is now the default.
- `threshold_const`, the constant `C` above, default 1. Calibrated by
  simulation; the optimum is flat over roughly [0.75, 1.5]. Applies only to
  `threshold_type = "theory"`, and warns if given with another type.
- `threshold` accepts a positive number as an explicit cutoff, in addition to
  `TRUE` and `FALSE`.
- `?fit_var` now documents all three cutoffs with their formulas, states that
  the rule applied is always hard thresholding, and states that thresholding
  controls over-selection rather than recovering the true support.

## Bug fixes

- `threshold` given a number was silently swallowed: `threshold = 0.05`
  returned bit-for-bit the `"soft"` answer rather than thresholding at 0.05,
  because any non-zero number tested as `TRUE` and the value was then dropped.
  Numeric cutoffs are now honoured, and `threshold = 0`, a negative number, a
  vector, or a non-numeric value are errors instead of being silently coerced.
- All three cutoffs now use the number of design rows, `T - p`, where they
  previously used `T`. This shifts the `"soft"` and `"hard"` cutoffs slightly —
  about 0.2% at `p = 1` and 1.0% at `p = 5`.

- `fit_var(..., scale = TRUE)` returned forecasts on the standardized scale
  while adding back the raw-scale column means, so they were wrong by a factor
  of each column's standard deviation. `transform_data()` now records the
  column standard deviations as `sd`, fitted objects carry them, and
  `compute_forecasts()` inverts the standardization before restoring the means.
- `test_granger()` now stops when given error bands computed with a resampling
  scheme other than `"bootstrap"`. Jackknife perturbations are much smaller
  than sampling variability, so the studentized statistic was inflated and the
  test over-rejected.
- `error_bands_irf(resampling = "bootstrapOLS")` labelled its output
  `"jackknife"`. It now reports `"bootstrapOLS"`.
- `method = "timeSlice"` selected its penalty parameter by scoring candidate
  coefficients, which are fitted on centered data, against the raw uncentered
  held-out observation. The validation error was therefore dominated by the
  series' means and the selected lambda depended on where the data happened to
  sit: adding a constant to each column of a simulated series moved lambda from
  0.00055 to 0.16, and the fitted coefficients by 0.54 in absolute value. The
  held-out observation is now put on the same scale as the fit. Results are
  unchanged for series that are already centered, and change for any series
  with a non-zero mean — which is most real data.
- SCAD and MCP fits with `method = "timeSlice"` computed their residuals from
  the raw, uncentered series rather than the centered series the coefficients
  were fitted on, inflating the residual standard deviation by about 60% on a
  pure location shift. They also returned `T` rows of residuals against a
  `$series` of `T - leaveOut` rows. Both now match the ENET time-slice path,
  which was already correct: residuals are computed on `$series` and have the
  same number of rows. `$sigma`, which is estimated from these residuals, was
  affected too.
- SCAD and MCP fits with `method = "timeSlice"` refitted at a single lambda,
  which `ncvreg` warns risks converging to an inferior local optimum. The final
  refit now runs along the full decreasing lambda path and selects the chosen
  value from it.
- Class checks use `inherits()` instead of comparing `attr(x, "class")` to a
  string, which errored on objects carrying more than one class.
- Removed a dead `center` override in `fit_vecm()` that never reached
  `fit_var()`. Behaviour is unchanged: the data was, and still is, centered.

# sparsevar 1.0.0

- Linted code: snake case for functions and variables
- Removed optional picasso dependency; SCAD estimation uses ncvreg only
- Minor bug fixing

# sparsevar 0.1.0

- Fix bug in plotIRF
- Linted code
- Fixed knitr/markdown/rmarkdown problems as in https://github.com/yihui/knitr/issues/1864

# sparsevar 0.0.11

- Added CV on a predefined list of lambdas (thanks to PierrickPiette) 

# sparsevar 0.0.10

- Added plotVECM function
- Removed plotComparisonVAR (substituted by plotVAR)
- Added the option to generate VARs with given matrices
- Fixed AIC
- Fixed problems with error bands options
- Added tests
- Added the function computeForecasts

# sparsevar 0.0.9

- Fast SCAD estimation for timeSlice SCAD
- Added function to compute VAR forecasts
- Added information criteria (AIC, Schwartz and Hannan-Quinn)
- Fixed mean estimation for timeSlice

# sparsevar 0.0.7

- Major code rewriting
- Remove dependencies from MTS and caret
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
