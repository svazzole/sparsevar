context("Test fit_var produces correct objects")

set.seed(123)

test_that("fitVAR SCAD via cv produces var object", {
  sim <- sparsevar::simulate_var(n = 6, p = 1, nobs = 40, sparsity = 0.1)
  fit <- sparsevar::fit_var(sim$series, p = 1, penalty = "SCAD", nfolds = 2)

  expect_s3_class(fit, "var")
  expect_equal(fit$penalty, "SCAD")
  expect_equal(fit$method, "cv")
  expect_length(fit$A, 1)
  expect_equal(dim(fit$A[[1]]), c(6, 6))
  expect_true(is.numeric(fit$lambda))
})

test_that("fitVAR MCP via cv produces var object", {
  sim <- sparsevar::simulate_var(n = 5, p = 1, nobs = 35, sparsity = 0.1)
  fit <- sparsevar::fit_var(sim$series, p = 1, penalty = "MCP", nfolds = 2)

  expect_s3_class(fit, "var")
  expect_equal(fit$penalty, "MCP")
  expect_equal(fit$method, "cv")
  expect_length(fit$A, 1)
  expect_equal(dim(fit$A[[1]]), c(5, 5))
  expect_true(is.numeric(fit$lambda))
})
