test_that("21.1: print.var returns a compact summary, not the raw list", {
  set.seed(101)
  sim <- sparsevar::simulate_var(n = 5, p = 2, nobs = 120, sparsity = 0.3)
  f <- sparsevar::fit_var(sim$series, p = 2, penalty = "ENET", method = "cv")

  out <- capture.output(print(f))

  # the defect this method exists to fix: the default print emitted 211 lines
  # of list internals for a 4-series VAR(1). Anything in this range is compact.
  expect_lt(length(out), 12)
  expect_match(out[1], "^Sparse VAR\\(2\\)$")
  expect_true(any(grepl("Series: 5", out, fixed = TRUE)))
  expect_true(any(grepl("Observations: 120", out, fixed = TRUE)))
  expect_true(any(grepl("ENET", out, fixed = TRUE)))
})

test_that("21.1: print methods return their argument invisibly", {
  set.seed(102)
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 100, sparsity = 0.3)
  f <- sparsevar::fit_var(sim$series, p = 1)

  # capture.output keeps the method's own output out of the test log
  capture.output(res <- withVisible(print(f)))

  expect_false(res$visible)
  expect_identical(res$value, f)
})

test_that("21.1: the reported non-zero count matches the coefficients", {
  set.seed(103)
  sim <- sparsevar::simulate_var(n = 5, p = 2, nobs = 120, sparsity = 0.3)
  f <- sparsevar::fit_var(sim$series, p = 2)

  nz <- sum(vapply(f$A, function(a) sum(a != 0), numeric(1)))
  total <- sum(vapply(f$A, length, integer(1)))

  out <- capture.output(print(f))
  line <- out[grepl("Non-zero coefficients", out)]

  expect_length(line, 1)
  expect_true(grepl(paste0(nz, " of ", total), line, fixed = TRUE))
})

test_that("21.1: print.var works on a simulated object, which has no penalty", {
  set.seed(104)
  sim <- sparsevar::simulate_var(n = 4, p = 1, nobs = 100, sparsity = 0.3)

  out <- capture.output(print(sim))

  # simulated objects carry no penalty/lambda, so that line must be omitted
  # rather than printed as NA
  expect_match(out[1], "^Sparse VAR\\(1\\)$")
  expect_false(any(grepl("Penalty", out, fixed = TRUE)))
  expect_false(any(grepl("NA", out, fixed = TRUE)))
})

test_that("21.3: summary.var reports density per lag, and they sum correctly", {
  set.seed(105)
  sim <- sparsevar::simulate_var(n = 5, p = 2, nobs = 120, sparsity = 0.3)
  f <- sparsevar::fit_var(sim$series, p = 2)
  s <- summary(f)

  expect_s3_class(s, "summary.var")
  expect_equal(nrow(s$density), 2L)
  expect_equal(s$density$non_zero,
               vapply(f$A, function(a) sum(a != 0), integer(1)))
  expect_equal(s$density$density, s$density$non_zero / s$density$entries)
})

test_that("21.3: summary residual_sd is per-equation, not the shrunk sigma", {
  set.seed(106)
  sim <- sparsevar::simulate_var(n = 5, p = 2, nobs = 200, sparsity = 0.3)
  f <- sparsevar::fit_var(sim$series, p = 2)
  s <- summary(f)

  expect_length(s$residual_sd, 5L)

  # corpcor::cov.shrink pulls the variances toward a common target, so the
  # diagonal of $sigma is (near) constant across equations. Reporting that as
  # a per-equation figure would be a constant dressed as an estimate, so
  # summary must NOT be using it.
  shrunk <- sqrt(diag(as.matrix(f$sigma)))
  expect_equal(stats::sd(shrunk), 0, tolerance = 1e-6)
  expect_gt(stats::sd(s$residual_sd), 1e-4)

  # and it must drop the first p rows, which hold series values not residuals
  expected <- apply(f$residuals[3:nrow(f$residuals), , drop = FALSE], 2,
                    stats::sd)
  expect_equal(unname(s$residual_sd), unname(expected))
})

test_that("21.2: print.vecm reports the rank of Pi", {
  set.seed(107)
  d <- matrix(cumsum(stats::rnorm(600)), ncol = 3)
  f <- suppressWarnings(sparsevar::fit_vecm(d, p = 1))

  out <- capture.output(print(f))

  expect_lt(length(out), 12)
  expect_match(out[1], "^Sparse VECM\\(1\\)$")
  expect_true(any(grepl("Rank of Pi", out, fixed = TRUE)))
  expect_true(any(grepl("Gamma", out, fixed = TRUE)))
})

test_that("21.2: print.varx separates endogenous from exogenous", {
  set.seed(108)
  sx <- sparsevar::simulate_varx(n = 4, k = 2, p = 1, m = 1, nobs = 120,
                                 sparsity_a1 = 0.3)
  f <- sparsevar::fit_varx(sx$series, p = 1, xt_matrix = sx$Xt, m = 1)

  # fit_varx() tags its result "var", not "varx" -- see the structural
  # delegation in print.var(). This test pins that a fitted VARX still prints
  # as a VARX, so the workaround cannot be removed without noticing.
  expect_identical(class(f), "var")

  out <- capture.output(print(f))

  expect_lt(length(out), 12)
  expect_match(out[1], "^Sparse VARX\\(1\\)$")
  # 4 endogenous and 2 exogenous, taken from $newA and $B rather than from
  # $series, which stacks both blocks into 6 columns
  expect_true(any(grepl("Endogenous series: 4", out, fixed = TRUE)))
  expect_true(any(grepl("Exogenous: 2", out, fixed = TRUE)))
})

test_that("21: the methods are registered for S3 dispatch", {
  expect_gte(length(methods(class = "var")), 2L)
  expect_true(any(grepl("print", format(methods(class = "var")))))
  expect_true(any(grepl("summary", format(methods(class = "var")))))
})
