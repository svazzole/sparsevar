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
  fit <- sparsevar::fitVAR(sim$series, p = 1)
  expect_output(str(fit), "List of 11")
  expect_output(cat(attr(fit, "class")), "var")
  expect_output(cat(attr(fit, "type")), "fit")
  
  sim <- sparsevar::simulate_var(n = 30, p = 1)
  fit <- sparsevar::fit_var(sim$series, p = 1, 
                           lambdas_list = c(0.1, 0.5, 0.3))
  expect_output(str(fit), "List of 11")
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
  expect_output(str(eb), "List of 8")
  expect_output(cat(attr(eb, "class")), "irfBands")
})
