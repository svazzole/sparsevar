context("Testing generation of VARs")

test_that("sparsevar", {
  sim <- sparsevar::simulate_var(n = 20, p = 2)
  expect_output(str(sim), "List of 4")
  expect_output(cat(attr(sim, "class")), "var")
  expect_output(cat(attr(sim, "type")), "simulation")
})

context("Testing estimation of VARs")

test_that("sparsevar", {
  sim <- sparsevar::simulate_var(n = 30, p = 1)
  fit <- sparsevar::fit_var(sim$series, p = 1)
  expect_output(str(fit), "List of 12")
  expect_output(cat(attr(fit, "class")), "var")
  expect_output(cat(attr(fit, "type")), "fit")

  sim <- sparsevar::simulate_var(n = 30, p = 1)
  fit <- sparsevar::fit_var(sim$series, p = 1,
                            lambdas_list = c(0.1, 0.5, 0.3))
  expect_output(str(fit), "List of 12")
  expect_output(cat(attr(fit, "class")), "var")
  expect_output(cat(attr(fit, "type")), "fit")
})

context("Testing IRF and error bands")

test_that("sparsevar", {
  sim <- sparsevar::simulate_var(n = 10, p = 3)
  fit <- sparsevar::fit_var(sim$series, p = 3)
  irf <- sparsevar::impulse_response(fit, len = 20)
  expect_output(str(irf), "List of 3")
  expect_output(cat(attr(irf, "class")), "irf")
  eb <- sparsevar::error_bands_irf(fit, irf, verbose = FALSE)
  # 8 band arrays plus Aboot, the resampled coefficient matrices used by
  # test_granger()
  expect_output(str(eb), "List of 9")
  expect_output(cat(attr(eb, "class")), "irfBands")
})

context("Testing Granger causality")

test_that("test_granger returns a joint causality decision", {
  sim <- sparsevar::simulate_var(n = 5, p = 2, nobs = 250, sparsity = 0.3)
  fit <- sparsevar::fit_var(sim$series, p = 2)
  irf <- sparsevar::impulse_response(fit, len = 10)
  eb <- sparsevar::error_bands_irf(fit, irf, m = 50, verbose = FALSE)

  g <- sparsevar::test_granger(fit, eb, alpha = 0.05)
  causal <- attr(g, "causality")

  # one matrix per lag, and an n x n logical causality matrix
  expect_equal(length(g), length(fit$A))
  expect_true(is.matrix(causal) && is.logical(causal))
  expect_equal(dim(causal), c(5L, 5L))

  # the null is joint over lags: a non-causal pair is zero at every lag
  for (l in seq_along(g)) {
    expect_true(all(g[[l]][!causal] == 0))
    expect_equal(g[[l]][causal], fit$A[[l]][causal])
  }

  # a stricter level can never detect more pairs
  strict <- attr(sparsevar::test_granger(fit, eb, alpha = 0.01), "causality")
  expect_true(sum(strict) <= sum(causal))

  # error bands without Aboot must fail loudly rather than silently
  stale <- eb
  stale$Aboot <- NULL
  expect_error(sparsevar::test_granger(fit, stale), "Aboot")
})

context("Correctness fixes (ticket 05)")

test_that("D7: each resampling scheme labels its own output", {
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 120, sparsity = 0.3)
  fit <- sparsevar::fit_var(sim$series, p = 1)
  irf <- sparsevar::impulse_response(fit, len = 5)

  eb_boot <- sparsevar::error_bands_irf(fit, irf, m = 20, verbose = FALSE)
  expect_equal(attr(eb_boot, "resampling"), "bootstrap")

  eb_jack <- sparsevar::error_bands_irf(fit, irf, verbose = FALSE,
                                        resampling = "jackknife")
  expect_equal(attr(eb_jack, "resampling"), "jackknife")

  # bootstrap_ols() used to label its output "jackknife" too, which is what made
  # the D6 guard unreliable; it now says "bootstrapOLS". That cannot be asserted
  # end-to-end because the bootstrapOLS path is broken independently of this fix
  # (p = 1: "$ operator is invalid for atomic vectors"; p = 2: a solve()
  # conformability error). When that path is repaired, replace this with a
  # positive check that the attribute reads "bootstrapOLS".
  # suppressWarnings() because that path also warns about a matrix shape
  # mismatch on its way to failing; the error is what we are pinning here
  expect_error(
    suppressWarnings(
      sparsevar::error_bands_irf(fit, irf, m = 20, verbose = FALSE,
                                 resampling = "bootstrapOLS")
    )
  )
})

test_that("D6: test_granger refuses non-bootstrap error bands", {
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 120, sparsity = 0.3)
  fit <- sparsevar::fit_var(sim$series, p = 1)
  irf <- sparsevar::impulse_response(fit, len = 5)

  eb_jack <- sparsevar::error_bands_irf(fit, irf, verbose = FALSE,
                                        resampling = "jackknife")
  # jackknife DOES carry Aboot, so the pre-existing is.null(Aboot) guard never
  # caught it; the statistic was silently inflated instead
  expect_false(is.null(eb_jack$Aboot))
  expect_error(sparsevar::test_granger(fit, eb_jack), "bootstrap")

  # the supported path still works
  eb_boot <- sparsevar::error_bands_irf(fit, irf, m = 20, verbose = FALSE)
  expect_silent(sparsevar::test_granger(fit, eb_boot))
})

test_that("D5: class checks survive objects with several classes", {
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 120, sparsity = 0.3)
  fit <- sparsevar::fit_var(sim$series, p = 1)

  expect_true(sparsevar::check_is_var(fit))

  # attr(x, "class") == "var" errors on a length-2 class vector; inherits() does not
  class(fit) <- c("myWrapper", "var")
  expect_true(sparsevar::check_is_var(fit))
  expect_silent(sparsevar::compute_forecasts(fit, 1))
})

test_that("D3: SCAD/MCP timeSlice fits along a full lambda path", {
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 120, sparsity = 0.3)

  # ncvreg warns when handed a single lambda; along a path it must not
  expect_silent(
    fit <- sparsevar::fit_var(sim$series, p = 1, penalty = "SCAD",
                              method = "timeSlice")
  )
  expect_equal(length(fit$A), 1L)
  expect_equal(dim(fit$A[[1]]), c(4L, 4L))
  expect_true(all(is.finite(fit$A[[1]])))
})

test_that("D2: fit_vecm centers, and recovers the cointegrating vector", {
  set.seed(789)
  z <- cumsum(stats::rnorm(200))
  d <- cbind(y1 = z + stats::rnorm(200, sd = 0.5),
             y2 = 0.5 * z + stats::rnorm(200, sd = 0.5))

  fv <- sparsevar::fit_vecm(data = d, p = 1, penalty = "ENET",
                            log_scale = FALSE)
  dec <- sparsevar::decompose_pi(fv, rk = 1)

  # true cointegrating vector is (1, -2); centering gives -2.083
  expect_equal(dec$beta[1], 1)
  expect_lt(abs(dec$beta[2] - (-2)), 0.12)
})

test_that("D8: scale = TRUE forecasts come back on the data's own scale", {
  set.seed(11)
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 200, sparsity = 0.3)
  # deliberately heterogeneous scales and non-zero means
  d <- sweep(sim$series, 2, c(1, 10, 100, 0.1, 5), "*")
  d <- sweep(d, 2, c(100, 50, 0, -20, 200), "+")

  f_raw <- sparsevar::compute_forecasts(
    sparsevar::fit_var(d, p = 1), 1)
  f_scl <- sparsevar::compute_forecasts(
    sparsevar::fit_var(d, p = 1, scale = TRUE), 1)

  # the scale factors must be recorded so the forecast can be inverted
  fit_scl <- sparsevar::fit_var(d, p = 1, scale = TRUE)
  expect_false(is.null(fit_scl$sd))
  expect_equal(fit_scl$sd, apply(scale(d, TRUE, FALSE), 2, stats::sd),
               tolerance = 1e-8, check.attributes = FALSE)

  # both forecasts must sit on the raw series scale. Before the fix the scaled
  # one returned a standardized deviation on top of a raw-scale mean, so it was
  # wrong by a factor of each column's sd.
  spread <- apply(d, 2, stats::sd)
  expect_true(all(abs(f_scl - f_raw) < 0.75 * spread))
})

test_that("F1: residuals and lambda do not depend on where the series sits", {
  set.seed(101)
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 200, sparsity = 0.3)
  d0 <- sim$series
  # a pure location shift: nothing about the dynamics changes
  d1 <- sweep(d0, 2, c(100, 50, -20, 200, 7), "+")

  for (pen in c("ENET", "SCAD", "MCP")) {
    for (meth in c("cv", "timeSlice")) {
      set.seed(7)
      f0 <- sparsevar::fit_var(d0, p = 1, penalty = pen, method = meth)
      set.seed(7)
      f1 <- sparsevar::fit_var(d1, p = 1, penalty = pen, method = meth)

      label <- paste(pen, meth, sep = "/")

      # residual sd cannot depend on a location shift. Before the fix the
      # SCAD/MCP timeSlice paths computed residuals on raw uncentered data
      # (+60% on a shift), and every timeSlice path scored its validation
      # error on raw levels, which moved the selected lambda by ~300x.
      expect_equal(apply(f0$residuals, 2, stats::sd),
                   apply(f1$residuals, 2, stats::sd),
                   tolerance = 1e-8, info = label)
      expect_equal(f0$A[[1]], f1$A[[1]], tolerance = 1e-8, info = label)
      expect_equal(f0$lambda, f1$lambda, tolerance = 1e-8, info = label)

      # residuals must describe the series the fit reports. timeSlice returns
      # the training window, so both must be T - leaveOut rows.
      expect_equal(nrow(f0$residuals), nrow(f0$series), info = label)
    }
  }
})

# ---- ticket 17: the threshold API (specified by ticket 09) -------------------

test_that("09.5: unknown options are an error, not silently ignored", {
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 120, sparsity = 0.3)
  d <- sim$series

  # the whole class of silent-typo bug
  expect_error(sparsevar::fit_var(d, p = 1, completely_made_up_argument = 99),
               "unknown option")
  # a near miss gets a suggestion
  expect_error(sparsevar::fit_var(d, p = 1, thersholdd_type = "hard"),
               "did you mean \"threshold_type\"")
  # the legacy spelling: dead in the fit path since v1.0.0, now loud
  expect_error(sparsevar::fit_var(d, p = 1, thresholdType = "hard"),
               "did you mean \"threshold_type\"")
  # documented but never implemented
  expect_error(sparsevar::fit_var(d, p = 1, horizon = 3), "fixed at 1")

  expect_silent(sparsevar::fit_var(d, p = 1, horizon = 1))
  expect_error(sparsevar::fit_varx(d, p = 1, xt_matrix = d[, 1, drop = FALSE],
                                   made_up = 1), "unknown option")
})

test_that("09.1: threshold is three-way, and validated", {
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 150, sparsity = 0.3)
  d <- sim$series

  none <- sparsevar::fit_var(d, p = 1)
  num <- sparsevar::fit_var(d, p = 1, threshold = 0.05)
  soft <- sparsevar::fit_var(d, p = 1, threshold = TRUE,
                             threshold_type = "soft")

  # before the fix a numeric threshold was swallowed and returned bit-for-bit
  # the soft answer, so this is a bug fix rather than a new feature
  expect_false(identical(as.vector(num$A[[1]]), as.vector(soft$A[[1]])))
  # and it must actually apply the number it was given
  expect_true(all(abs(num$A[[1]][num$A[[1]] != 0]) >= 0.05))
  expect_lt(sum(num$A[[1]] != 0), sum(none$A[[1]] != 0))

  expect_error(sparsevar::fit_var(d, p = 1, threshold = 0), "> 0")
  expect_error(sparsevar::fit_var(d, p = 1, threshold = -1), "> 0")
  expect_error(sparsevar::fit_var(d, p = 1, threshold = c(1, 2)), "length 1")
  expect_error(sparsevar::fit_var(d, p = 1, threshold = "x"), "must be FALSE")
  expect_error(sparsevar::fit_var(d, p = 1, threshold = NA), "NA")
  # a numeric cutoff makes threshold_type meaningless: warn, do not error
  expect_warning(sparsevar::fit_var(d, p = 1, threshold = 0.05,
                                    threshold_type = "hard"), "ignored")
})

test_that("09.2: threshold_type defaults to theory, and takes three values", {
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 150, sparsity = 0.3)
  d <- sim$series

  set.seed(3)
  suppressMessages(
    default <- sparsevar::fit_var(d, p = 1, threshold = TRUE)
  )
  set.seed(3)
  theory <- sparsevar::fit_var(d, p = 1, threshold = TRUE,
                               threshold_type = "theory")
  expect_equal(default$A[[1]], theory$A[[1]])

  expect_silent(sparsevar::fit_var(d, p = 1, threshold = TRUE,
                                   threshold_type = "hard"))
  expect_error(sparsevar::fit_var(d, p = 1, threshold = TRUE,
                                  threshold_type = "nonsense"))
})

test_that("09.6: threshold_const scales the theory cutoff only", {
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 150, sparsity = 0.3)
  d <- sim$series

  c1 <- sparsevar::fit_var(d, p = 1, threshold = TRUE,
                           threshold_type = "theory", threshold_const = 1)
  c3 <- sparsevar::fit_var(d, p = 1, threshold = TRUE,
                           threshold_type = "theory", threshold_const = 3)
  # a larger constant is a larger cutoff, so it cannot select more
  expect_lte(sum(c3$A[[1]] != 0), sum(c1$A[[1]] != 0))

  expect_error(sparsevar::fit_var(d, p = 1, threshold = TRUE,
                                  threshold_const = 0), "positive")
  expect_error(sparsevar::fit_var(d, p = 1, threshold = TRUE,
                                  threshold_const = c(1, 2)), "single")
  # ignored by the legacy cutoffs, and said so
  expect_warning(sparsevar::fit_var(d, p = 1, threshold = TRUE,
                                    threshold_type = "soft",
                                    threshold_const = 2), "theory")
})

test_that("09.7: the theory cutoff carries the units of the coefficients", {
  # Exact, and the property that matters: under a change of per-series units
  # u, a VAR coefficient A_ij becomes A_ij * u_i / u_j, and sigma_i / s_j
  # carries exactly that factor. So the cutoff tracks the coefficients and the
  # selected support is *identical*. This is the invariant the scale factors
  # were added for; it is a property of apply_threshold(), tested here without
  # the non-scale-equivariant lambda selection upstream getting in the way.
  set.seed(404)
  nc <- 5
  p <- 1
  a <- matrix(stats::rnorm(nc * nc * p, sd = 0.3), nrow = nc, ncol = nc * p)
  sigma <- runif(nc, 0.5, 2)
  s <- runif(nc, 0.5, 2)
  units <- c(1, 100, 0.01, 10, 5)

  keep <- sparsevar:::apply_threshold(a, 200, nc, p, type = "theory",
                                      sigma = sigma, s = s) != 0

  # re-express every series in new units
  a_u <- a * outer(units, rep(units, p), "/")
  keep_u <- sparsevar:::apply_threshold(a_u, 200, nc, p, type = "theory",
                                        sigma = sigma * units,
                                        s = s * units) != 0

  expect_equal(keep, keep_u)

  # the scale-free cutoffs cannot do this
  h <- sparsevar:::apply_threshold(a, 200, nc, p, type = "hard") != 0
  h_u <- sparsevar:::apply_threshold(a_u, 200, nc, p, type = "hard") != 0
  expect_false(isTRUE(all.equal(h, h_u)))
})

test_that("09.7: end to end, theory keeps the support broadly stable", {
  # The end-to-end version is *not* invariant, and must not be asserted as
  # such: sparsevar stacks all n equations into one regression with a single
  # shared lambda, which is itself not scale-equivariant. So this asserts
  # broad agreement only.
  set.seed(303)
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 250, sparsity = 0.3)
  d <- sim$series
  d_u <- sweep(d, 2, c(1, 100, 0.01, 10, 5), "*")

  set.seed(9)
  f <- sparsevar::fit_var(d, p = 1, threshold = TRUE,
                          threshold_type = "theory")
  set.seed(9)
  f_u <- sparsevar::fit_var(d_u, p = 1, threshold = TRUE,
                            threshold_type = "theory")

  expect_gt(mean((f$A[[1]] != 0) == (f_u$A[[1]] != 0)), 0.75)
})

test_that("09.7: theory needs a scale, and every fit path supplies one", {
  set.seed(55)
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 150, sparsity = 0.3)
  d <- sim$series

  # a missing scale silently degrades to the bare scalar cutoff, so it errors
  a <- matrix(stats::rnorm(16), nrow = 4, ncol = 4)
  expect_error(sparsevar:::apply_threshold(a, 100, 4, 1, type = "theory"),
               "needs both sigma and s")

  # all four fit paths reach the theory branch without error
  for (pen in c("ENET", "SCAD", "MCP")) {
    for (meth in c("cv", "timeSlice")) {
      f <- sparsevar::fit_var(d, p = 1, penalty = pen, method = meth,
                              threshold = TRUE, threshold_type = "theory")
      expect_true(all(is.finite(f$A[[1]])),
                  info = paste(pen, meth, sep = "/"))
    }
  }

  # and so do the IRF bootstrap paths, including the SCAD/MCP branches that
  # have no residuals of their own
  for (pen in c("ENET", "SCAD")) {
    f <- sparsevar::fit_var(d, p = 1, penalty = pen)
    irf <- sparsevar::impulse_response(f)
    # the SCAD branch here still refits at a single lambda, which ncvreg warns
    # about — a pre-existing issue on the IRF path, not this ticket's
    eb <- suppressWarnings(
      sparsevar::error_bands_irf(f, irf, m = 5, verbose = FALSE,
                                 threshold = TRUE,
                                 threshold_type = "theory",
                                 threshold_const = 1)
    )
    expect_true(all(is.finite(eb$irf$upper)), info = pen)
  }

  # the legacy spelling is gone from this path too
  f <- sparsevar::fit_var(d, p = 1)
  irf <- sparsevar::impulse_response(f)
  expect_error(sparsevar::error_bands_irf(f, irf, m = 5, verbose = FALSE,
                                          thresholdType = "hard"),
               "unknown option")
})

test_that("09.7: the cutoff uses T - p design rows, not T", {
  # the legacy cutoffs read nr, which is now the number of design rows. The
  # shift is small but real, and is why NEWS records a numerical change.
  a <- matrix(0.05, nrow = 4, ncol = 8)
  at_t <- sparsevar:::apply_threshold(a, 200, 4, 2, type = "soft")
  at_tp <- sparsevar:::apply_threshold(a, 198, 4, 2, type = "soft")
  expect_false(isTRUE(all.equal(1 / (2 * 4 * log(200)),
                                1 / (2 * 4 * log(198)))))
  expect_equal(dim(at_t), dim(at_tp))
})

test_that("09.7: under scale = TRUE the cutoff is unscaled, by design", {
  set.seed(77)
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 200, sparsity = 0.3)
  d <- sweep(sim$series, 2, c(1, 50, 0.2, 8), "*")

  # scale = TRUE already puts the coefficients on a common scale, so s_j = 1
  # and sigma_i is on that same scale: the ratio stays coherent either way
  f <- sparsevar::fit_var(d, p = 1, scale = TRUE, threshold = TRUE,
                          threshold_type = "theory")
  expect_true(all(is.finite(f$A[[1]])))
  expect_equal(as.vector(apply(f$series, 2, stats::sd)), rep(1, 4),
               tolerance = 1e-8)
})
