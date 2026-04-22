# Test Hessian accuracy after selective zeroing optimization
#
# This test verifies that the Hessian computation is correct after
# implementing selective zeroing of Hessian arrays (only zeroing
# free parameter entries instead of full nparam×nparam arrays).

test_that("Hessian is correct with selective zeroing optimization", {
  skip_on_cran()

  set.seed(42)

  # Create simple data
  n <- 200
  f <- rnorm(n)
  x1 <- rnorm(n)

  # Simulate outcomes
  T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
  T2 <- 1.5 + 0.8*f + rnorm(n, 0, 0.6)
  Y <- 1.0 + 0.5*x1 + 0.7*f + rnorm(n, 0, 0.4)

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, Y = Y, eval = 1)

  # Create model with fixed loading for identification
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0,  # Fixed loading for identification
    evaluation_indicator = "eval"
  )

  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Free loading
    evaluation_indicator = "eval"
  )

  mc_Y <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_,  # Free loading
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_Y), factor = fm)

  # Get initial parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Get the full parameter vector including names
  full_params <- init_result$init_params

  # All parameters are free (factor variance can be estimated because T1 has fixed loading)
  param_fixed <- rep(FALSE, length(full_params))

  # Use the actual parameter values from initialization as test point
  test_params <- full_params

  cat("\nNumber of parameters:", length(test_params), "\n")
  cat("Parameter names:", paste(names(test_params), collapse = ", "), "\n")

  # Check Hessian accuracy
  hess_check <- check_hessian_accuracy(ms, dat, test_params,
                                       param_fixed = param_fixed,
                                       tol = 1e-3, verbose = TRUE, n_quad = 8)

  expect_true(hess_check$pass,
              info = sprintf("Hessian check failed with max error: %.2e", hess_check$max_error))

  # Also check gradient for completeness
  grad_check <- check_gradient_accuracy(ms, dat, test_params,
                                        param_fixed = param_fixed,
                                        tol = 1e-3, verbose = TRUE, n_quad = 8)

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed with max error: %.2e", grad_check$max_error))
})
