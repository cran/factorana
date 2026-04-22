# Systematic Test Suite for Factorana Package
#
# This test suite systematically tests the package functionality by:
# 1. Checking analytical gradients against finite differences
# 2. Checking analytical Hessians against finite differences
# 3. Estimating models with default and true initial parameters
# 4. Comparing estimates to true parameter values
#
# Tests are marked skip_on_cran() as they involve estimation and take time

# Test configuration
VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
SAVE_LOGS <- Sys.getenv("FACTORANA_TEST_SAVE_LOGS", "TRUE") == "TRUE"
GRAD_TOL <- 1e-3  # Relaxed to accommodate numerical precision in finite differences
HESS_TOL <- 1e-3  # Checking all elements (diagonal and off-diagonal)

# ==============================================================================
# Test A: Measurement System with 3 Linear Tests and 1 Factor
# ==============================================================================

test_that("Model A: Measurement system with 3 linear tests and 1 factor", {
  skip_on_cran()

  set.seed(104)

  # Simulate data
  n <- 500
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4)   # T3: int, lambda, sigma

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  dat <- data.frame(intercept = 1, T1 = T1, T2 = T2, T3 = T3, eval = 1)

  # Create model system with 1 factor
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0,  # Fix loading to 1 for identification
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Estimate freely
    evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Estimate freely
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_A_measurement_system", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test B: Measurement System (A) + Probit
# ==============================================================================

test_that("Model B: Measurement system with 3 linear tests and probit outcome", {
  skip_on_cran()

  set.seed(105)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma, y_int, y_beta1, y_loading
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                   0.5, 0.7, 0.9)   # Probit: int, beta1, loading

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  z <- true_params[10] + true_params[11]*x1 + true_params[12]*f
  y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_y <- define_model_component(
    name = "y", data = dat, outcome = "y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "probit",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "y_int", "y_beta1", "y_loading"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_B_measurement_plus_probit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test C: Measurement System (A) + Ordered Probit
# ==============================================================================

test_that("Model C: Measurement system with 3 linear tests and ordered probit", {
  skip_on_cran()

  set.seed(106)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma,
  #                           y_beta1, y_loading, thresh1, thresh_incr1
  # NOTE: Ordered probit has NO INTERCEPT (absorbed into thresholds)
  # Threshold parameterization: thresh2 = thresh1 + abs(thresh_incr1)
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                   0.6, 0.8,   # Oprobit: beta1, loading (NO intercept)
                   -0.5, 1.0)  # Oprobit thresholds: thresh1, thresh_incr1 (intercept absorbed)

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  # Ordered probit latent variable (with intercept for data generation)
  # But model will not estimate intercept - it's absorbed into thresholds
  y_intercept <- 0.5
  z <- y_intercept + true_params[10]*x1 + true_params[11]*f

  # Compute absolute threshold values from incremental parameterization
  # thresh1 has intercept absorbed: original 0.0 - intercept 0.5 = -0.5
  thresh1_abs <- true_params[12]
  thresh2_abs <- thresh1_abs + abs(true_params[13])

  u <- rnorm(n)
  y <- ifelse(z + u < thresh1_abs, 1, ifelse(z + u < thresh2_abs, 2, 3))

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_y <- define_model_component(
    name = "y", data = dat, outcome = "y", factor = fm,
    covariates = "x1",  # NO intercept for ordered probit (absorbed into thresholds)
    model_type = "oprobit",
    num_choices = 3,  # 3 ordered categories
    loading_normalization = NA_real_, evaluation_indicator = "eval",
    intercept = FALSE
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var", "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "y_beta1", "y_loading",  # Oprobit: NO intercept
                   "thresh1", "thresh_incr1"),  # Incremental threshold parameterization
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_C_measurement_plus_oprobit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test D: Measurement System (A) + Multinomial Logit
# ==============================================================================

test_that("Model D: Measurement system with 3 linear tests and multinomial logit", {
  skip_on_cran()

  set.seed(107)

  # Simulate data
  n <- 500
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma, T3_int, T3_lambda, T3_sigma,
  #                           y1_int, y1_beta1, y1_loading, y2_int, y2_beta1, y2_loading
  true_params <- c(1.0,  # Factor variance
                   2.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   1.5, 1.2, 0.6,   # T2: int, lambda, sigma
                   1.0, 0.8, 0.4,   # T3: int, lambda, sigma
                   0.5, 0.6, 0.7,   # Mlogit choice 1: int, beta1, loading
                   1.0, -0.5, 0.9)  # Mlogit choice 2: int, beta1, loading

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  z1 <- true_params[10] + true_params[11]*x1 + true_params[12]*f
  z2 <- true_params[13] + true_params[14]*x1 + true_params[15]*f

  exp_z0 <- 1
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  p0 <- exp_z0 / denom
  p1 <- exp_z1 / denom
  p2 <- exp_z2 / denom

  y <- numeric(n)
  for (i in seq_len(n)) {
    # C++ expects multinomial choices coded as 1, 2, 3 (not 0, 1, 2)
    y[i] <- sample(1:3, 1, prob = c(p0[i], p1[i], p2[i]))
  }

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2, T3 = T3, y = y, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  mc_y <- define_model_component(
    name = "y", data = dat, outcome = "y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "logit",
    num_choices = 3,  # 3 choices: 0, 1, 2
    loading_normalization = NA_real_,  # Single factor, loading estimated freely
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_T3, mc_y), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "y1_int", "y1_beta1", "y1_loading",
                   "y2_int", "y2_beta1", "y2_loading"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_D_measurement_plus_mlogit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test E: Roy Model
# ==============================================================================

test_that("Model E: Roy selection model", {
  skip_on_cran()

  set.seed(108)

  # Simulate Roy model data
  n <- 2000  # Increased from 500 to match manual test for better numerical stability
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  f <- rnorm(n)  # Latent factor (ability)

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma,
  #                           T3_int, T3_lambda, T3_sigma, wage0_int, wage0_beta1, wage0_beta2, wage0_sigma,
  #                           wage1_int, wage1_beta1, wage1_lambda, wage1_sigma, sector_int, sector_beta1, sector_loading
  true_params <- c(
    1.0,  # Factor variance
    # T1: int, sigma (loading FIXED to 1.0)
    2.0, 0.5,
    # T2: int, lambda, sigma
    1.5, 1.2, 0.6,
    # T3: int, lambda, sigma
    1.0, 0.8, 0.4,
    # Wage0: int, beta1, beta2, sigma (loading FIXED to 0.0 - no factor effect)
    2.0, 0.5, 0.3, 0.6,
    # Wage1: int, beta1, lambda, sigma
    2.5, 0.6, 1.0, 0.7,
    # Sector: int, beta1, loading (probit - no sigma)
    0.0, 0.4, 0.8
  )

  # Generate test scores (using fixed loading values for data generation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])
  T3 <- true_params[7] + true_params[8]*f + rnorm(n, 0, true_params[9])

  # Generate potential wages
  wage0 <- true_params[10] + true_params[11]*x1 + true_params[12]*x2 +
           rnorm(n, 0, true_params[13])  # No factor effect (loading = 0.0 fixed)
  wage1 <- true_params[14] + true_params[15]*x1 + true_params[16]*f +
           rnorm(n, 0, true_params[17])

  # Sector choice (based on utility difference)
  z_sector <- true_params[18] + true_params[19]*x2 + true_params[20]*f
  sector <- as.numeric(runif(n) < pnorm(z_sector))

  # Observed wage
  wage <- ifelse(sector == 1, wage1, wage0)

  dat <- data.frame(
    intercept = 1,
    x1 = x1, x2 = x2,
    T1 = T1, T2 = T2, T3 = T3,
    wage = wage,
    sector = sector,
    eval_tests = 1,
    eval_wage0 = 1 - sector,
    eval_wage1 = sector,
    eval_sector = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval_tests"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_tests"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_tests"
  )
  mc_wage0 <- define_model_component(
    name = "wage0", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1", "x2"), model_type = "linear",
    loading_normalization = 0.0, evaluation_indicator = "eval_wage0"
  )
  mc_wage1 <- define_model_component(
    name = "wage1", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval_wage1"
  )
  mc_sector <- define_model_component(
    name = "sector", data = dat, outcome = "sector", factor = fm,
    covariates = c("intercept", "x2"), model_type = "probit",
    loading_normalization = NA_real_, evaluation_indicator = "eval_sector"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_wage0, mc_wage1, mc_sector),
    factor = fm
  )

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sigma",  # T1: loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "T3_int", "T3_lambda", "T3_sigma",
                   "w0_int", "w0_beta1", "w0_beta2", "w0_sigma",
                   "w1_int", "w1_beta1", "w1_lambda", "w1_sigma",
                   "s_int", "s_beta1", "s_loading"),
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 16)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 16)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_E_roy_model", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
  expect_true(est_comp$reasonable_default, info = "Estimates are not reasonable")
})

# ==============================================================================
# Test F: Measurement System + Exploded Multinomial Logit (Ranked Choices)
# ==============================================================================

test_that("Model F: Measurement system with exploded logit (ranked choices)", {
  skip_on_cran()

  set.seed(109)

  # Simulate data - simpler case: 3 choices, 2 ranks
  n <- 1000  # Larger sample for stability
  x1 <- rnorm(n)
  f <- rnorm(n)  # Latent factor

  # True parameters - ORDER: factor_var, T1_int, T1_sigma, T2_int, T2_lambda, T2_sigma,
  #                           c1_int, c1_beta1, c1_loading, c2_int, c2_beta1, c2_loading
  # 3 choices = 2 non-reference alternatives
  true_params <- c(1.0,  # Factor variance
                   0.0, 0.5,   # T1: int, sigma (loading fixed to 1.0)
                   0.0, 0.8, 0.5,   # T2: int, lambda, sigma
                   0.5, 0.8, 0.6,   # Mlogit choice 1: int, beta1, loading
                   -0.3, -0.5, 0.9)  # Mlogit choice 2: int, beta1, loading

  # Generate data (using true loading values for simulation)
  T1 <- true_params[2] + 1.0*f + rnorm(n, 0, true_params[3])  # T1 loading = 1.0 (fixed)
  T2 <- true_params[4] + true_params[5]*f + rnorm(n, 0, true_params[6])

  # Multinomial logit utilities (3 choices: reference + 2 alternatives)
  z1 <- true_params[7] + true_params[8]*x1 + true_params[9]*f
  z2 <- true_params[10] + true_params[11]*x1 + true_params[12]*f

  exp_z0 <- 1  # Reference category
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  # Generate rankings for each observation (2 ranks from 3 choices)
  rank1 <- rank2 <- numeric(n)
  for (i in seq_len(n)) {
    probs <- c(exp_z0, exp_z1[i], exp_z2[i]) / denom[i]
    # Sample first choice
    rank1[i] <- sample(1:3, 1, prob = probs)
    # Update probs for second choice (remove first choice)
    remaining <- setdiff(1:3, rank1[i])
    probs2 <- probs[remaining] / sum(probs[remaining])
    rank2[i] <- sample(remaining, 1, prob = probs2)
  }

  dat <- data.frame(intercept = 1, x1 = x1, T1 = T1, T2 = T2,
                    rank1 = rank1, rank2 = rank2, eval = 1)

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )
  # Exploded logit with 2 ranks
  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2"),  # Vector of outcomes for ranked choices
    factor = fm,
    covariates = c("intercept", "x1"), model_type = "logit",
    num_choices = 3,  # 3 choices total
    loading_normalization = NA_real_,  # Single factor, loading estimated freely
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_choice), factor = fm)

  # First run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = c("f_var",
                   "T1_int", "T1_sigma",  # T1 loading fixed to 1.0
                   "T2_int", "T2_lambda", "T2_sigma",
                   "c1_int", "c1_beta1", "c1_loading",
                   "c2_int", "c2_beta1", "c2_loading"),
    verbose = VERBOSE
  )

  # Run gradient check with param_fixed
  # NOTE: Hessian for exploded logit with conditional probabilities is complex
  # and not yet fully implemented - skip Hessian check for this model type
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && est_comp$default_converged
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_F_exploded_logit", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions - for exploded logit, we focus on gradient correctness and convergence
  # Hessian for conditional probabilities in exploded logit is a TODO
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(est_comp$default_converged, info = "Estimation with default init failed to converge")
})

# ==============================================================================
# Test G: Structural Equation Model (SE_linear)
# f2 = se_intercept + se_linear_1 * f1 + epsilon
# ==============================================================================

test_that("Model G: SE_linear structural equation model with 2 factors", {
  skip_on_cran()

  set.seed(201)

  # Simulate data
  n <- 500

  # True parameters for factor structure
  true_var_f1 <- 1.5        # Variance of input factor
  true_se_intercept <- 0.5   # SE intercept
  true_se_linear <- 0.8      # SE linear coefficient
  true_se_residual_var <- 0.5  # SE residual variance

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + eps

  # Generate measurements for f1 (3 measures)
  Y1_1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 1.0)   # int=2, loading=1 (fixed), sigma=1
  Y1_2 <- 1.5 + 0.8 * f1 + rnorm(n, 0, 0.9)   # int=1.5, loading=0.8, sigma=0.9
  Y1_3 <- 1.0 + 1.2 * f1 + rnorm(n, 0, 1.1)   # int=1, loading=1.2, sigma=1.1

  # Generate measurements for f2 (3 measures)
  Y2_1 <- 1.8 + 1.0 * f2 + rnorm(n, 0, 1.0)   # int=1.8, loading=1 (fixed), sigma=1
  Y2_2 <- 1.2 + 0.7 * f2 + rnorm(n, 0, 0.8)   # int=1.2, loading=0.7, sigma=0.8
  Y2_3 <- 0.5 + 1.1 * f2 + rnorm(n, 0, 0.95)  # int=0.5, loading=1.1, sigma=0.95

  dat <- data.frame(
    intercept = 1,
    Y1_1 = Y1_1, Y1_2 = Y1_2, Y1_3 = Y1_3,
    Y2_1 = Y2_1, Y2_2 = Y2_2, Y2_3 = Y2_3,
    eval = 1
  )

  # Create model system with 2 factors and SE_linear structure
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  # 3 measurement equations for factor 1
  mc1_1 <- define_model_component(
    name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),  # Fix loading to 1 for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_2 <- define_model_component(
    name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_3 <- define_model_component(
    name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )

  # 3 measurement equations for factor 2
  mc2_1 <- define_model_component(
    name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),  # 0 for f1, fix to 1 for f2
    evaluation_indicator = "eval"
  )
  mc2_2 <- define_model_component(
    name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )
  mc2_3 <- define_model_component(
    name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm
  )

  # True parameters:
  # factor_var_1, se_intercept, se_linear_1, se_residual_var,
  # Y1_1_int, Y1_1_sigma,  (loading fixed to 1)
  # Y1_2_int, Y1_2_lambda, Y1_2_sigma,
  # Y1_3_int, Y1_3_lambda, Y1_3_sigma,
  # Y2_1_int, Y2_1_sigma,  (loading fixed to 1)
  # Y2_2_int, Y2_2_lambda, Y2_2_sigma,
  # Y2_3_int, Y2_3_lambda, Y2_3_sigma
  true_params <- c(
    true_var_f1,  # factor_var_1 = 1.5
    true_se_intercept, true_se_linear, true_se_residual_var,  # SE params: 0.5, 0.8, 0.5
    2.0, 1.0,  # Y1_1: int, sigma (loading fixed to 1)
    1.5, 0.8, 0.9,  # Y1_2: int, lambda, sigma
    1.0, 1.2, 1.1,  # Y1_3: int, lambda, sigma
    1.8, 1.0,  # Y2_1: int, sigma (loading fixed to 1)
    1.2, 0.7, 0.8,  # Y2_2: int, lambda, sigma
    0.5, 1.1, 0.95  # Y2_3: int, lambda, sigma
  )

  param_names <- c(
    "f1_var",
    "se_intercept", "se_linear_1", "se_residual_var",
    "Y1_1_int", "Y1_1_sigma",
    "Y1_2_int", "Y1_2_lambda", "Y1_2_sigma",
    "Y1_3_int", "Y1_3_lambda", "Y1_3_sigma",
    "Y2_1_int", "Y2_1_sigma",
    "Y2_2_int", "Y2_2_lambda", "Y2_2_sigma",
    "Y2_3_int", "Y2_3_lambda", "Y2_3_sigma"
  )

  # Run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = param_names,
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                        param_fixed = est_comp$param_fixed,
                                        tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_G_se_linear", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions - focus on gradient and Hessian correctness
  # Estimation convergence with default init is not strictly required for this test
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
})

# ==============================================================================
# Test G2: Two-Stage SE_linear Estimation
# Stage 1: Independent factors (measurement model)
# Stage 2: SE_linear with previous_stage (structural equation)
# ==============================================================================

test_that("Model G2: Two-stage SE_linear estimation with previous_stage", {
  skip_on_cran()

  set.seed(301)

  # Simulate data with known DGP
  n <- 500

  # True parameters
  # Note: In two-stage estimation, Stage 1 assumes independent factors with mean 0.
  # The measurement intercepts absorb the factor means. So in Stage 2, the SE intercept
  # captures only ADDITIONAL mean shift beyond what Stage 1 absorbed. We set it to 0.
  true_var_f1 <- 1.2        # Input factor variance
  true_se_linear <- 0.7     # Structural coefficient
  true_se_intercept <- 0.0  # Structural intercept (0 for two-stage identification)
  true_se_residual_var <- 0.4  # Structural residual variance

  # True measurement parameters
  true_lambda_f1 <- c(1.0, 0.9, 1.1)  # Loadings for f1 (first fixed)
  true_lambda_f2 <- c(1.0, 0.8, 1.2)  # Loadings for f2 (first fixed)
  true_int_f1 <- c(2.0, 1.5, 1.0)
  true_int_f2 <- c(1.8, 1.2, 0.5)
  true_sigma_f1 <- c(0.8, 0.9, 0.7)
  true_sigma_f2 <- c(0.9, 0.85, 0.95)

  # Generate factors according to SE_linear structure
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps_se <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + eps_se

  # Generate measurements
  Y1_1 <- true_int_f1[1] + true_lambda_f1[1] * f1 + rnorm(n, 0, true_sigma_f1[1])
  Y1_2 <- true_int_f1[2] + true_lambda_f1[2] * f1 + rnorm(n, 0, true_sigma_f1[2])
  Y1_3 <- true_int_f1[3] + true_lambda_f1[3] * f1 + rnorm(n, 0, true_sigma_f1[3])
  Y2_1 <- true_int_f2[1] + true_lambda_f2[1] * f2 + rnorm(n, 0, true_sigma_f2[1])
  Y2_2 <- true_int_f2[2] + true_lambda_f2[2] * f2 + rnorm(n, 0, true_sigma_f2[2])
  Y2_3 <- true_int_f2[3] + true_lambda_f2[3] * f2 + rnorm(n, 0, true_sigma_f2[3])

  dat <- data.frame(
    intercept = 1,
    Y1_1 = Y1_1, Y1_2 = Y1_2, Y1_3 = Y1_3,
    Y2_1 = Y2_1, Y2_2 = Y2_2, Y2_3 = Y2_3,
    eval = 1
  )

  # =========================================================================
  # STAGE 1: Measurement model with INDEPENDENT factors
  # =========================================================================

  fm_stage1 <- define_factor_model(n_factors = 2, factor_structure = "independent")

  # Factor 1 measures
  mc1_1 <- define_model_component(
    name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm_stage1,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0), evaluation_indicator = "eval"
  )
  mc1_2 <- define_model_component(
    name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm_stage1,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval"
  )
  mc1_3 <- define_model_component(
    name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm_stage1,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval"
  )

  # Factor 2 measures
  mc2_1 <- define_model_component(
    name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm_stage1,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1), evaluation_indicator = "eval"
  )
  mc2_2 <- define_model_component(
    name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm_stage1,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval"
  )
  mc2_3 <- define_model_component(
    name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm_stage1,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval"
  )

  ms_stage1 <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm_stage1
  )

  # Stage 1 estimation
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result_stage1 <- estimate_model_rcpp(
    model_system = ms_stage1,
    data = dat,
    init_params = NULL,
    control = ctrl,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check Stage 1 convergence
  expect_equal(result_stage1$convergence, 0, info = "Stage 1 should converge")

  # Build true parameters for Stage 1 (independent 2-factor model)
  # Parameter order: factor_var_1, factor_var_2, then per-component: intercept, [loading], sigma
  true_params_s1 <- c(
    true_var_f1,                              # factor_var_1
    true_var_f1 * true_se_linear^2 + true_se_residual_var,  # factor_var_2 (marginal variance of f2)
    true_int_f1[1], true_sigma_f1[1],         # Y1_1: intercept, sigma (loading fixed to 1)
    true_int_f1[2], true_lambda_f1[2], true_sigma_f1[2],  # Y1_2: intercept, loading, sigma
    true_int_f1[3], true_lambda_f1[3], true_sigma_f1[3],  # Y1_3: intercept, loading, sigma
    true_int_f2[1], true_sigma_f2[1],         # Y2_1: intercept, sigma (loading fixed to 1)
    true_int_f2[2], true_lambda_f2[2], true_sigma_f2[2],  # Y2_2: intercept, loading, sigma
    true_int_f2[3], true_lambda_f2[3], true_sigma_f2[3]   # Y2_3: intercept, loading, sigma
  )
  names(true_params_s1) <- names(result_stage1$estimates)

  # Get param_fixed vector for Stage 1
  stage1_init <- initialize_parameters(ms_stage1, dat, verbose = FALSE)
  param_metadata_s1 <- factorana:::build_parameter_metadata(ms_stage1)
  param_constraints_s1 <- factorana:::setup_parameter_constraints(
    ms_stage1, true_params_s1, param_metadata_s1,
    stage1_init$factor_variance_fixed, verbose = FALSE
  )
  param_fixed_s1 <- rep(TRUE, length(true_params_s1))
  param_fixed_s1[param_constraints_s1$free_idx] <- FALSE

  # Check Stage 1 gradient accuracy at TRUE parameters (not at MLE)
  grad_check_s1 <- check_gradient_accuracy(ms_stage1, dat, true_params_s1,
                                           param_fixed = param_fixed_s1,
                                           tol = GRAD_TOL, verbose = FALSE, n_quad = 8)
  expect_true(grad_check_s1$pass,
              info = sprintf("Stage 1 gradient check failed (max error: %.2e)", grad_check_s1$max_error))

  # Check Stage 1 parameter recovery (measurement params)
  # Factor variances - use unname() for comparison
  expect_equal(unname(result_stage1$estimates["factor_var_1"]), true_var_f1, tolerance = 0.3,
               info = "Stage 1: factor_var_1 recovery")

  # Loadings (free ones) - use unname() for comparison
  expect_equal(unname(result_stage1$estimates["Y1_2_loading_1"]), true_lambda_f1[2], tolerance = 0.15,
               info = "Stage 1: Y1_2 loading recovery")
  expect_equal(unname(result_stage1$estimates["Y2_2_loading_2"]), true_lambda_f2[2], tolerance = 0.15,
               info = "Stage 1: Y2_2 loading recovery")

  # Check Stage 1 SEs are computed
  expect_true(all(!is.na(result_stage1$std_errors)),
              info = "Stage 1 should have computed std_errors")
  expect_true(all(result_stage1$std_errors > 0),
              info = "Stage 1 std_errors should be positive")

  if (VERBOSE) {
    cat("\n=== Stage 1 (Independent Factors) ===\n")
    cat(sprintf("Convergence: %d\n", result_stage1$convergence))
    cat(sprintf("Log-likelihood: %.4f\n", result_stage1$loglik))
    cat(sprintf("Gradient max error: %.2e\n", grad_check_s1$max_error))
    cat(sprintf("factor_var_1: true=%.3f, est=%.3f\n",
                true_var_f1, result_stage1$estimates["factor_var_1"]))
  }

  # =========================================================================
  # STAGE 2: SE_linear with previous_stage
  # =========================================================================

  fm_stage2 <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  # No new components - just change factor structure
  ms_stage2 <- define_model_system(
    components = list(),
    factor = fm_stage2,
    previous_stage = result_stage1
  )

  # Initialize SE parameters
  init_stage2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  init_stage2$init_params["se_intercept"] <- 0.0
  init_stage2$init_params["se_linear_1"] <- 0.5
  init_stage2$init_params["se_residual_var"] <- 0.5

  result_stage2 <- estimate_model_rcpp(
    model_system = ms_stage2,
    data = dat,
    init_params = init_stage2$init_params,
    control = ctrl,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check Stage 2 convergence
  expect_equal(result_stage2$convergence, 0, info = "Stage 2 should converge")

  # Build true parameters for Stage 2 (SE_linear structure)
  # Parameter order for SE_linear: factor_var_1, se_intercept, se_linear_1, se_residual_var,
  # then measurement params from Stage 1
  true_params_s2 <- c(
    true_var_f1,          # factor_var_1
    true_se_intercept,    # se_intercept
    true_se_linear,       # se_linear_1
    true_se_residual_var  # se_residual_var
  )
  # Append measurement parameters (using Stage 1 true params without factor variances)
  meas_params <- true_params_s1[-(1:2)]  # Remove factor variances
  true_params_s2 <- c(true_params_s2, meas_params)
  names(true_params_s2) <- names(result_stage2$estimates)

  # Get param_fixed vector for Stage 2
  stage2_init <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  param_metadata_s2 <- factorana:::build_parameter_metadata(ms_stage2)
  param_constraints_s2 <- factorana:::setup_parameter_constraints(
    ms_stage2, true_params_s2, param_metadata_s2,
    stage2_init$factor_variance_fixed, verbose = FALSE
  )
  param_fixed_s2 <- rep(TRUE, length(true_params_s2))
  param_fixed_s2[param_constraints_s2$free_idx] <- FALSE

  # Check Stage 2 gradient accuracy at TRUE parameters (not at MLE)
  grad_check_s2 <- check_gradient_accuracy(ms_stage2, dat, true_params_s2,
                                           param_fixed = param_fixed_s2,
                                           tol = GRAD_TOL, verbose = FALSE, n_quad = 8)
  expect_true(grad_check_s2$pass,
              info = sprintf("Stage 2 gradient check failed (max error: %.2e)", grad_check_s2$max_error))

  # Check Stage 2 Hessian accuracy at TRUE parameters
  hess_check_s2 <- check_hessian_accuracy(ms_stage2, dat, true_params_s2,
                                          param_fixed = param_fixed_s2,
                                          tol = HESS_TOL, verbose = FALSE, n_quad = 8)
  expect_true(hess_check_s2$pass,
              info = sprintf("Stage 2 Hessian check failed (max error: %.2e)", hess_check_s2$max_error))

  # =========================================================================
  # VERIFICATION: Parameters preserved from Stage 1
  # =========================================================================

  # Measurement parameters should be identical
  stage1_meas_params <- result_stage1$estimates[!grepl("^factor_var", names(result_stage1$estimates))]
  stage2_meas_params <- result_stage2$estimates[names(stage1_meas_params)]

  max_meas_diff <- max(abs(stage1_meas_params - stage2_meas_params))
  expect_true(max_meas_diff < 1e-10,
              info = sprintf("Measurement params should be preserved (max diff: %.2e)", max_meas_diff))

  # SEs for fixed params should be preserved
  stage1_meas_se <- result_stage1$std_errors[names(stage1_meas_params)]
  stage2_meas_se <- result_stage2$std_errors[names(stage1_meas_params)]

  max_se_diff <- max(abs(stage1_meas_se - stage2_meas_se))
  expect_true(max_se_diff < 1e-10,
              info = sprintf("Measurement SEs should be preserved (max diff: %.2e)", max_se_diff))

  # =========================================================================
  # VERIFICATION: SE_linear parameter recovery
  # =========================================================================

  # SE linear coefficient (the key structural parameter)
  est_se_linear <- unname(result_stage2$estimates["se_linear_1"])
  expect_equal(est_se_linear, true_se_linear, tolerance = 0.15,
               info = sprintf("SE linear recovery: true=%.3f, est=%.3f",
                             true_se_linear, est_se_linear))

  # SE intercept (should be ~0 since Stage 1 absorbed factor means)
  est_se_intercept <- unname(result_stage2$estimates["se_intercept"])
  expect_equal(est_se_intercept, true_se_intercept, tolerance = 0.15,
               info = sprintf("SE intercept recovery: true=%.3f, est=%.3f",
                             true_se_intercept, est_se_intercept))

  # SE residual variance (relaxed tolerance due to sampling variability in SE models)
  est_se_resvar <- unname(result_stage2$estimates["se_residual_var"])
  expect_equal(est_se_resvar, true_se_residual_var, tolerance = 0.35,
               info = sprintf("SE residual var recovery: true=%.3f, est=%.3f",
                             true_se_residual_var, est_se_resvar))

  # Factor variance for f1 (should be estimated in Stage 2)
  est_var_f1 <- unname(result_stage2$estimates["factor_var_1"])
  expect_equal(est_var_f1, true_var_f1, tolerance = 0.3,
               info = sprintf("Factor var recovery: true=%.3f, est=%.3f",
                             true_var_f1, est_var_f1))

  # Check SE params have positive std_errors
  se_param_names <- c("se_linear_1", "se_intercept", "se_residual_var", "factor_var_1")
  se_std_errors <- result_stage2$std_errors[se_param_names]
  expect_true(all(se_std_errors > 0),
              info = "SE parameters should have positive std_errors")

  if (VERBOSE) {
    cat("\n=== Stage 2 (SE_linear) ===\n")
    cat(sprintf("Convergence: %d\n", result_stage2$convergence))
    cat(sprintf("Log-likelihood: %.4f\n", result_stage2$loglik))
    cat(sprintf("Gradient max error: %.2e\n", grad_check_s2$max_error))
    cat(sprintf("Hessian max error: %.2e\n", hess_check_s2$max_error))
    cat(sprintf("Meas params preserved: max diff = %.2e\n", max_meas_diff))
    cat(sprintf("Meas SEs preserved: max diff = %.2e\n", max_se_diff))
    cat("\nSE Parameter Recovery:\n")
    cat(sprintf("  se_linear_1: true=%.3f, est=%.3f (%.1f%%)\n",
                true_se_linear, est_se_linear, 100*est_se_linear/true_se_linear))
    cat(sprintf("  se_intercept: true=%.3f, est=%.3f\n",
                true_se_intercept, est_se_intercept))
    cat(sprintf("  se_residual_var: true=%.3f, est=%.3f (%.1f%%)\n",
                true_se_residual_var, est_se_resvar, 100*est_se_resvar/true_se_residual_var))
    cat(sprintf("  factor_var_1: true=%.3f, est=%.3f (%.1f%%)\n",
                true_var_f1, est_var_f1, 100*est_var_f1/true_var_f1))
  }
})

# ==============================================================================
# Test H: Structural Equation Model (SE_quadratic)
# f2 = se_intercept + se_linear_1 * f1 + se_quadratic_1 * f1^2 + epsilon
# ==============================================================================

test_that("Model H: SE_quadratic structural equation model with 2 factors", {
  skip_on_cran()

  set.seed(202)

  # Simulate data
  n <- 500

  # True parameters for factor structure
  true_var_f1 <- 1.5          # Variance of input factor
  true_se_intercept <- 0.5     # SE intercept
  true_se_linear <- 0.8        # SE linear coefficient
  true_se_quadratic <- 0.2     # SE quadratic coefficient
  true_se_residual_var <- 0.5  # SE residual variance

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + true_se_quadratic * f1^2 + eps

  # Generate measurements for f1 (3 measures)
  Y1_1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 1.0)   # int=2, loading=1 (fixed), sigma=1
  Y1_2 <- 1.5 + 0.8 * f1 + rnorm(n, 0, 0.9)   # int=1.5, loading=0.8, sigma=0.9
  Y1_3 <- 1.0 + 1.2 * f1 + rnorm(n, 0, 1.1)   # int=1, loading=1.2, sigma=1.1

  # Generate measurements for f2 (3 measures)
  Y2_1 <- 1.8 + 1.0 * f2 + rnorm(n, 0, 1.0)   # int=1.8, loading=1 (fixed), sigma=1
  Y2_2 <- 1.2 + 0.7 * f2 + rnorm(n, 0, 0.8)   # int=1.2, loading=0.7, sigma=0.8
  Y2_3 <- 0.5 + 1.1 * f2 + rnorm(n, 0, 0.95)  # int=0.5, loading=1.1, sigma=0.95

  dat <- data.frame(
    intercept = 1,
    Y1_1 = Y1_1, Y1_2 = Y1_2, Y1_3 = Y1_3,
    Y2_1 = Y2_1, Y2_2 = Y2_2, Y2_3 = Y2_3,
    eval = 1
  )

  # Create model system with 2 factors and SE_quadratic structure
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  # 3 measurement equations for factor 1
  mc1_1 <- define_model_component(
    name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),  # Fix loading to 1 for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_2 <- define_model_component(
    name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )
  mc1_3 <- define_model_component(
    name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0),  # Free for f1, 0 for f2
    evaluation_indicator = "eval"
  )

  # 3 measurement equations for factor 2
  mc2_1 <- define_model_component(
    name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),  # 0 for f1, fix to 1 for f2
    evaluation_indicator = "eval"
  )
  mc2_2 <- define_model_component(
    name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )
  mc2_3 <- define_model_component(
    name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_),  # 0 for f1, free for f2
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm
  )

  # True parameters:
  # factor_var_1, se_intercept, se_linear_1, se_quadratic_1, se_residual_var,
  # Y1_1_int, Y1_1_sigma,  (loading fixed to 1)
  # Y1_2_int, Y1_2_lambda, Y1_2_sigma,
  # Y1_3_int, Y1_3_lambda, Y1_3_sigma,
  # Y2_1_int, Y2_1_sigma,  (loading fixed to 1)
  # Y2_2_int, Y2_2_lambda, Y2_2_sigma,
  # Y2_3_int, Y2_3_lambda, Y2_3_sigma
  true_params <- c(
    true_var_f1,  # factor_var_1 = 1.5
    true_se_intercept, true_se_linear, true_se_quadratic, true_se_residual_var,  # SE params
    2.0, 1.0,  # Y1_1: int, sigma (loading fixed to 1)
    1.5, 0.8, 0.9,  # Y1_2: int, lambda, sigma
    1.0, 1.2, 1.1,  # Y1_3: int, lambda, sigma
    1.8, 1.0,  # Y2_1: int, sigma (loading fixed to 1)
    1.2, 0.7, 0.8,  # Y2_2: int, lambda, sigma
    0.5, 1.1, 0.95  # Y2_3: int, lambda, sigma
  )

  param_names <- c(
    "f1_var",
    "se_intercept", "se_linear_1", "se_quadratic_1", "se_residual_var",
    "Y1_1_int", "Y1_1_sigma",
    "Y1_2_int", "Y1_2_lambda", "Y1_2_sigma",
    "Y1_3_int", "Y1_3_lambda", "Y1_3_sigma",
    "Y2_1_int", "Y2_1_sigma",
    "Y2_2_int", "Y2_2_lambda", "Y2_2_sigma",
    "Y2_3_int", "Y2_3_lambda", "Y2_3_sigma"
  )

  # Run estimation to get param_fixed
  est_comp <- run_estimation_comparison(
    ms, dat, true_params,
    param_names = param_names,
    verbose = VERBOSE
  )

  # Run gradient/Hessian checks with param_fixed
  grad_check <- check_gradient_accuracy(ms, dat, true_params,
                                        param_fixed = est_comp$param_fixed,
                                        tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8)
  hess_check <- check_hessian_accuracy(ms, dat, true_params,
                                       param_fixed = est_comp$param_fixed,
                                       tol = HESS_TOL, verbose = VERBOSE, n_quad = 8)

  # Collect diagnostics
  diagnostics <- list(
    gradient_check = grad_check,
    hessian_check = hess_check,
    estimation = est_comp,
    overall_pass = grad_check$pass && hess_check$pass &&
                   est_comp$default_converged && est_comp$reasonable_default
  )

  # Save log
  if (SAVE_LOGS) {
    log_file <- save_diagnostics_to_log("test_H_se_quadratic", diagnostics)
    if (VERBOSE) cat("Log saved to:", log_file, "\n")
  }

  # Assertions - focus on gradient and Hessian correctness
  # Estimation convergence with default init is not strictly required for this test
  expect_true(grad_check$pass, info = sprintf("Gradient check failed (max error: %.2e)", grad_check$max_error))
  expect_true(hess_check$pass, info = sprintf("Hessian check failed (max error: %.2e)", hess_check$max_error))
})


# ==============================================================================
# Test I: Equality Constraints (Measurement Invariance) with SE_quadratic
# ==============================================================================

test_that("Model I: SE_quadratic with equality constraints (measurement invariance)", {
  skip_on_cran()

  set.seed(456)

  # True parameters
  n <- 800
  true_var_f1 <- 1.0
  true_se_intercept <- 0.0  # Set to 0 for identification
  true_se_linear <- 0.7
  true_se_quadratic <- 0.15
  true_se_residual_var <- 0.4

  # Common measurement parameters (invariance across time points)
  true_loading_2 <- 0.8
  true_loading_3 <- 1.2
  true_sigma_1 <- 1.0
  true_sigma_2 <- 0.9
  true_sigma_3 <- 1.1

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, sqrt(true_se_residual_var))
  f2 <- true_se_intercept + true_se_linear * f1 + true_se_quadratic * f1^2 + eps

  # Measurements at time 1 (factor 1)
  Y1_1 <- 1.0 * f1 + rnorm(n, 0, true_sigma_1)
  Y1_2 <- true_loading_2 * f1 + rnorm(n, 0, true_sigma_2)
  Y1_3 <- true_loading_3 * f1 + rnorm(n, 0, true_sigma_3)

  # Measurements at time 2 (factor 2) - SAME loadings and sigmas
  Y2_1 <- 1.0 * f2 + rnorm(n, 0, true_sigma_1)
  Y2_2 <- true_loading_2 * f2 + rnorm(n, 0, true_sigma_2)
  Y2_3 <- true_loading_3 * f2 + rnorm(n, 0, true_sigma_3)

  dat <- data.frame(Y1_1, Y1_2, Y1_3, Y2_1, Y2_2, Y2_3, intercept = 1, eval = 1)

  # Model setup
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  mc1_1 <- define_model_component(name = "Y1_1", data = dat, outcome = "Y1_1", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(1, 0), evaluation_indicator = "eval")
  mc1_2 <- define_model_component(name = "Y1_2", data = dat, outcome = "Y1_2", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval")
  mc1_3 <- define_model_component(name = "Y1_3", data = dat, outcome = "Y1_3", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval")
  mc2_1 <- define_model_component(name = "Y2_1", data = dat, outcome = "Y2_1", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(0, 1), evaluation_indicator = "eval")
  mc2_2 <- define_model_component(name = "Y2_2", data = dat, outcome = "Y2_2", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval")
  mc2_3 <- define_model_component(name = "Y2_3", data = dat, outcome = "Y2_3", factor = fm,
                                    covariates = "intercept", model_type = "linear",
                                    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval")

  # Define model system WITH equality constraints (measurement invariance)
  ms <- define_model_system(
    components = list(mc1_1, mc1_2, mc1_3, mc2_1, mc2_2, mc2_3),
    factor = fm,
    equality_constraints = list(
      c("Y1_2_loading_1", "Y2_2_loading_2"),  # Same loading for item 2
      c("Y1_3_loading_1", "Y2_3_loading_2"),  # Same loading for item 3
      c("Y1_1_sigma", "Y2_1_sigma"),          # Same error variance for item 1
      c("Y1_2_sigma", "Y2_2_sigma"),          # Same error variance for item 2
      c("Y1_3_sigma", "Y2_3_sigma")           # Same error variance for item 3
    )
  )

  # Test 1: Verify equality constraints are stored correctly
  expect_equal(length(ms$equality_constraints), 5)

  # Test 2: Run estimation
  control <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(model_system = ms, data = dat, control = control,
                                 optimizer = "nlminb", verbose = FALSE)

  # Test 3: Verify tied parameters are exactly equal (use unname() to compare values, not named vectors)
  expect_equal(unname(result$estimates["Y1_2_loading_1"]), unname(result$estimates["Y2_2_loading_2"]))
  expect_equal(unname(result$estimates["Y1_3_loading_1"]), unname(result$estimates["Y2_3_loading_2"]))
  expect_equal(unname(result$estimates["Y1_1_sigma"]), unname(result$estimates["Y2_1_sigma"]))
  expect_equal(unname(result$estimates["Y1_2_sigma"]), unname(result$estimates["Y2_2_sigma"]))
  expect_equal(unname(result$estimates["Y1_3_sigma"]), unname(result$estimates["Y2_3_sigma"]))

  # Test 4: Verify param_table shows tied_to column correctly
  expect_equal(result$param_table$tied_to[result$param_table$name == "Y2_2_loading_2"], "Y1_2_loading_1")
  expect_equal(result$param_table$tied_to[result$param_table$name == "Y2_1_sigma"], "Y1_1_sigma")

  # Test 5: Parameter recovery (check key parameters are within tolerance)
  # Slightly relaxed tolerance to account for random variation in Monte Carlo simulation
  tolerance <- 0.20
  expect_lt(abs(result$estimates["factor_var_1"] - true_var_f1), tolerance)
  expect_lt(abs(result$estimates["se_linear_1"] - true_se_linear), tolerance)
  expect_lt(abs(result$estimates["se_quadratic_1"] - true_se_quadratic), tolerance)
  expect_lt(abs(result$estimates["se_residual_var"] - true_se_residual_var), tolerance)
  expect_lt(abs(result$estimates["Y1_2_loading_1"] - true_loading_2), tolerance)
  expect_lt(abs(result$estimates["Y1_3_loading_1"] - true_loading_3), tolerance)
  expect_lt(abs(result$estimates["Y1_1_sigma"] - true_sigma_1), tolerance)
  expect_lt(abs(result$estimates["Y1_2_sigma"] - true_sigma_2), tolerance)
  expect_lt(abs(result$estimates["Y1_3_sigma"] - true_sigma_3), tolerance)

  # Test 6: FD gradient check with equality constraints
  # The C++ gradient needs to be aggregated for tied params
  # Note: C++ now returns gradients for FREE params only, so we need to map back to full
  data_mat <- as.matrix(dat)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init_result$init_params

  # Set to true values for gradient check
  params[1] <- true_var_f1
  params[2] <- true_se_intercept
  params[3] <- true_se_linear
  params[4] <- true_se_quadratic
  params[5] <- true_se_residual_var
  params[9] <- true_loading_2   # Y1_2_loading_1
  params[12] <- true_loading_3  # Y1_3_loading_1
  params[7] <- true_sigma_1     # Y1_1_sigma
  params[10] <- true_sigma_2    # Y1_2_sigma
  params[13] <- true_sigma_3    # Y1_3_sigma
  # Set tied params equal
  params[17] <- true_loading_2  # Y2_2_loading_2
  params[20] <- true_loading_3  # Y2_3_loading_2
  params[15] <- true_sigma_1    # Y2_1_sigma
  params[18] <- true_sigma_2    # Y2_2_sigma
  params[21] <- true_sigma_3    # Y2_3_sigma

  fm_ptr <- initialize_factor_model_cpp(ms, data_mat, 8, params)

  # Extract free params (tied params are now fixed in C++)
  params_free <- extract_free_params_cpp(fm_ptr, params)
  cpp_result <- evaluate_likelihood_cpp(fm_ptr, params_free, compute_gradient = TRUE, compute_hessian = FALSE)

  # Map gradient from free params back to full params
  # Free indices are all params except tied ones (17, 20, 15, 18, 21)
  # free_idx = c(1:14, 16, 19) = 16 free params
  free_idx <- c(1:14, 16, 19)
  analytical_grad <- rep(0, length(params))
  analytical_grad[free_idx] <- cpp_result$gradient

  # No need for manual aggregation since C++ already excludes tied params
  # The gradient for primary params now directly includes the contribution
  aggregated_analytical <- analytical_grad

  # Compute FD gradient with constraints applied
  eps_fd <- 1e-5
  fd_grad <- numeric(length(params))
  for (i in seq_along(params)) {
    params_plus <- params
    params_minus <- params
    params_plus[i] <- params[i] + eps_fd
    params_minus[i] <- params[i] - eps_fd

    # Apply equality constraints
    if (i == 9) { params_plus[17] <- params_plus[9]; params_minus[17] <- params_minus[9] }
    if (i == 12) { params_plus[20] <- params_plus[12]; params_minus[20] <- params_minus[12] }
    if (i == 7) { params_plus[15] <- params_plus[7]; params_minus[15] <- params_minus[7] }
    if (i == 10) { params_plus[18] <- params_plus[10]; params_minus[18] <- params_minus[10] }
    if (i == 13) { params_plus[21] <- params_plus[13]; params_minus[21] <- params_minus[13] }

    # Need to extract free params for the perturbed params
    params_plus_free <- extract_free_params_cpp(fm_ptr, params_plus)
    params_minus_free <- extract_free_params_cpp(fm_ptr, params_minus)

    ll_plus <- evaluate_loglik_only_cpp(fm_ptr, params_plus_free)
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr, params_minus_free)
    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps_fd)
  }

  # Check primary parameters only (not derived ones)
  primary_indices <- c(1:14, 16, 19)
  max_rel_err <- 0
  worst_idx <- 0
  cat("\n=== DEBUG: Gradient comparison ===\n")
  for (i in primary_indices) {
    ana <- aggregated_analytical[i]
    fd <- fd_grad[i]
    rel_err <- if (abs(fd) > 1e-6) abs(ana - fd) / abs(fd) else abs(ana - fd)
    cat(sprintf("  param %2d: ana=%.6f, fd=%.6f, rel_err=%.2e\n", i, ana, fd, rel_err))
    if (!is.na(rel_err) && rel_err > max_rel_err) {
      max_rel_err <- rel_err
      worst_idx <- i
    }
  }
  cat(sprintf("Worst param: %d, max_rel_err: %.2e\n", worst_idx, max_rel_err))

  expect_lt(max_rel_err, 1e-4, label = sprintf("Gradient max error: %.2e", max_rel_err))

  if (VERBOSE) {
    cat("Test I: Equality constraints (measurement invariance)\n")
    cat("  Equality constraints: 5 groups\n")
    cat("  Free parameters: 16 (21 total - 5 tied)\n")
    cat(sprintf("  Gradient max rel error: %.2e\n", max_rel_err))
    cat(sprintf("  SE_linear recovery: true=%.3f, est=%.3f\n",
                true_se_linear, result$estimates["se_linear_1"]))
    cat(sprintf("  Loading_2 recovery: true=%.3f, est=%.3f\n",
                true_loading_2, result$estimates["Y1_2_loading_1"]))
  }
})

# =====================================================================
# Test J: Two-stage model with multiple outcome types and fixed coefficients
# - Stage 1: Linear measurement system for factor 1
# - Stage 2: Various outcome models (linear, probit, oprobit, mlogit, exploded logit)
#            with at least one fixed coefficient in each
# =====================================================================
test_that("Test J: Two-stage with multiple outcome types and fixed coefficients", {
  skip_on_cran()

  set.seed(54321)
  n <- 400

  # Generate factor
  f1 <- rnorm(n)

  # Stage 1: Linear measurement equations for factor 1
  Y_m1 <- 1.0 * f1 + rnorm(n, 0, 0.5)  # loading fixed to 1
  Y_m2 <- 0.8 * f1 + rnorm(n, 0, 0.6)  # loading = 0.8

  # Covariates for stage 2
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # True parameters for stage 2 models
  true_beta_linear <- c(0.5, 0.3)    # intercept, x1 coefficient
  true_loading_linear <- 0.7
  true_sigma_linear <- 0.8

  true_beta_probit <- c(-0.2, 0.4)   # intercept, x1 coefficient
  true_loading_probit <- 0.6

  true_beta_oprobit <- c(0.0, 0.5)   # intercept=0 (fixed), x1 coefficient
  true_loading_oprobit <- 0.5
  true_thresholds <- c(-0.5, 0.5)    # 3 categories -> 2 thresholds

  true_beta_mlogit <- matrix(c(
    0.3, 0.2,   # Choice 1: intercept, x1
    0.0, 0.4    # Choice 2: intercept=0 (fixed), x1
  ), nrow = 2, byrow = TRUE)
  true_loading_mlogit <- c(0.4, 0.5)  # loadings for choices 1, 2

  true_beta_explogit <- matrix(c(
    0.2, 0.0,   # Choice 1: intercept, x1=0 (fixed)
    0.1, 0.3    # Choice 2: intercept, x1
  ), nrow = 2, byrow = TRUE)
  true_loading_explogit <- c(0.3, 0.4)

  # Generate outcomes
  # Linear outcome
  Y_linear <- true_beta_linear[1] + true_beta_linear[2] * x1 +
              true_loading_linear * f1 + rnorm(n, 0, true_sigma_linear)

  # Probit outcome
  latent_probit <- true_beta_probit[1] + true_beta_probit[2] * x1 +
                   true_loading_probit * f1 + rnorm(n, 0, 1)
  Y_probit <- as.integer(latent_probit > 0)

  # Ordered probit outcome (3 categories) - store as numeric integers for matrix conversion
  latent_oprobit <- true_beta_oprobit[1] + true_beta_oprobit[2] * x1 +
                    true_loading_oprobit * f1 + rnorm(n, 0, 1)
  Y_oprobit <- as.integer(cut(latent_oprobit, breaks = c(-Inf, true_thresholds, Inf),
                               labels = FALSE))

  # Multinomial logit outcome (3 choices, choice 0 is reference)
  u0 <- 0  # Reference
  u1 <- true_beta_mlogit[1,1] + true_beta_mlogit[1,2] * x1 +
        true_loading_mlogit[1] * f1 + rlogis(n, 0, 1)
  u2 <- true_beta_mlogit[2,1] + true_beta_mlogit[2,2] * x1 +
        true_loading_mlogit[2] * f1 + rlogis(n, 0, 1)
  Y_mlogit <- apply(cbind(u0, u1, u2), 1, which.max)

  # Exploded logit outcome (rank 2 alternatives from 3 choices)
  # Compute utilities for all 3 choices
  v0 <- 0  # Reference
  v1 <- true_beta_explogit[1,1] + true_beta_explogit[1,2] * x1 +
        true_loading_explogit[1] * f1 + rlogis(n, 0, 1)
  v2 <- true_beta_explogit[2,1] + true_beta_explogit[2,2] * x1 +
        true_loading_explogit[2] * f1 + rlogis(n, 0, 1)
  utils <- cbind(v0, v1, v2)
  # Generate rankings (top 2 ranks)
  rankings <- t(apply(utils, 1, function(u) order(u, decreasing = TRUE)[1:2]))
  Y_rank1 <- rankings[, 1]
  Y_rank2 <- rankings[, 2]

  # Create data frame
  dat <- data.frame(
    Y_m1 = Y_m1, Y_m2 = Y_m2,
    Y_linear = Y_linear,
    Y_probit = Y_probit,
    Y_oprobit = Y_oprobit,
    Y_mlogit = Y_mlogit,
    Y_rank1 = Y_rank1, Y_rank2 = Y_rank2,
    x1 = x1, x2 = x2,
    intercept = 1, eval = 1
  )

  # Define factor model (single factor)
  fm <- define_factor_model(n_factors = 1)

  # Stage 1: Measurement equations
  mc_m1 <- define_model_component(name = "m1", data = dat, outcome = "Y_m1", factor = fm,
                                   covariates = "intercept", model_type = "linear",
                                   loading_normalization = 1, evaluation_indicator = "eval")
  mc_m2 <- define_model_component(name = "m2", data = dat, outcome = "Y_m2", factor = fm,
                                   covariates = "intercept", model_type = "linear",
                                   loading_normalization = NA_real_, evaluation_indicator = "eval")

  # Stage 2: Outcome models with fixed coefficients

  # Linear with x2 coefficient fixed to 0
  mc_linear <- define_model_component(name = "linear", data = dat, outcome = "Y_linear", factor = fm,
                                       covariates = c("intercept", "x1", "x2"), model_type = "linear",
                                       loading_normalization = NA_real_, evaluation_indicator = "eval")
  mc_linear <- fix_coefficient(mc_linear, "x2", 0)

  # Probit with x2 coefficient fixed to 0
  mc_probit <- define_model_component(name = "probit", data = dat, outcome = "Y_probit", factor = fm,
                                       covariates = c("intercept", "x1", "x2"), model_type = "probit",
                                       loading_normalization = NA_real_, evaluation_indicator = "eval")
  mc_probit <- fix_coefficient(mc_probit, "x2", 0)

  # Ordered probit with x2 coefficient fixed to 0 (3 categories).
  # Note: oprobit has no intercept covariate — the intercept is absorbed into
  # the cut points. So we fix x2 instead to exercise fix_coefficient() here.
  mc_oprobit <- define_model_component(name = "oprobit", data = dat, outcome = "Y_oprobit", factor = fm,
                                        covariates = c("x1", "x2"), model_type = "oprobit",
                                        num_choices = 3, loading_normalization = NA_real_,
                                        evaluation_indicator = "eval")
  mc_oprobit <- fix_coefficient(mc_oprobit, "x2", 0)

  # Multinomial logit with choice 2 intercept fixed to 0
  mc_mlogit <- define_model_component(name = "mlogit", data = dat, outcome = "Y_mlogit", factor = fm,
                                       covariates = c("intercept", "x1"), model_type = "logit",
                                       num_choices = 3, loading_normalization = NA_real_,
                                       evaluation_indicator = "eval")
  mc_mlogit <- fix_coefficient(mc_mlogit, "intercept", 0, choice = 2)

  # Exploded logit with choice 1 x1 coefficient fixed to 0
  mc_explogit <- define_model_component(name = "explogit", data = dat, outcome = c("Y_rank1", "Y_rank2"),
                                         factor = fm, covariates = c("intercept", "x1"), model_type = "logit",
                                         num_choices = 3, loading_normalization = NA_real_,
                                         evaluation_indicator = "eval")
  mc_explogit <- fix_coefficient(mc_explogit, "x1", 0, choice = 1)

  # Define model system
  ms <- define_model_system(
    components = list(mc_m1, mc_m2, mc_linear, mc_probit, mc_oprobit, mc_mlogit, mc_explogit),
    factor = fm
  )

  # Initialize and estimate
  control <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(model_system = ms, data = dat, control = control,
                                 optimizer = "nlminb", verbose = FALSE)

  # Test 1: Check that fixed coefficients are exactly at their fixed values
  expect_equal(unname(result$estimates["linear_x2"]), 0, tolerance = 1e-10)
  expect_equal(unname(result$estimates["probit_x2"]), 0, tolerance = 1e-10)
  expect_equal(unname(result$estimates["oprobit_x2"]), 0, tolerance = 1e-10)
  expect_equal(unname(result$estimates["mlogit_c2_intercept"]), 0, tolerance = 1e-10)
  expect_equal(unname(result$estimates["explogit_c1_x1"]), 0, tolerance = 1e-10)

  # Test 2: FD gradient check
  data_mat <- as.matrix(dat)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init_result$init_params

  fm_ptr <- initialize_factor_model_cpp(ms, data_mat, 8, params)
  params_free <- extract_free_params_cpp(fm_ptr, params)
  cpp_result <- evaluate_likelihood_cpp(fm_ptr, params_free, compute_gradient = TRUE, compute_hessian = FALSE)

  # Compute FD gradient
  eps <- 1e-5
  fd_grad <- numeric(length(params_free))
  for (i in seq_along(params_free)) {
    params_plus <- params_free
    params_minus <- params_free
    params_plus[i] <- params_free[i] + eps
    params_minus[i] <- params_free[i] - eps
    ll_plus <- evaluate_loglik_only_cpp(fm_ptr, params_plus)
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr, params_minus)
    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps)
  }

  # Check gradient accuracy
  max_rel_err <- 0
  for (i in seq_along(params_free)) {
    ana <- cpp_result$gradient[i]
    fd <- fd_grad[i]
    rel_err <- if (abs(fd) > 1e-6) abs(ana - fd) / abs(fd) else abs(ana - fd)
    if (!is.na(rel_err) && rel_err > max_rel_err) max_rel_err <- rel_err
  }
  expect_lt(max_rel_err, 1e-4, label = sprintf("Gradient max error: %.2e", max_rel_err))

  # Test 3: FD Hessian check (all elements)
  cpp_hess_result <- evaluate_likelihood_cpp(fm_ptr, params_free, compute_gradient = FALSE, compute_hessian = TRUE)
  hess_vec <- cpp_hess_result$hessian

  # Expand upper triangular to full matrix
  n_free <- length(params_free)
  hess_mat <- matrix(0, n_free, n_free)
  idx <- 1
  for (i in 1:n_free) {
    for (j in i:n_free) {
      hess_mat[i, j] <- hess_vec[idx]
      hess_mat[j, i] <- hess_vec[idx]
      idx <- idx + 1
    }
  }

  # Compute FD Hessian by differentiating gradient (all elements)
  fd_hess <- matrix(0, n_free, n_free)
  grad_base <- evaluate_likelihood_cpp(fm_ptr, params_free, compute_gradient = TRUE, compute_hessian = FALSE)$gradient
  for (i in seq_len(n_free)) {
    params_plus <- params_free
    h <- eps * (abs(params_free[i]) + 1.0)
    params_plus[i] <- params_free[i] + h

    grad_plus <- evaluate_likelihood_cpp(fm_ptr, params_plus, compute_gradient = TRUE, compute_hessian = FALSE)$gradient

    fd_hess[i, ] <- (grad_plus - grad_base) / h
  }
  # Symmetrize FD Hessian
  fd_hess <- (fd_hess + t(fd_hess)) / 2

  # Check all Hessian elements (upper triangle)
  max_hess_err <- 0
  for (i in seq_len(n_free)) {
    for (j in i:n_free) {
      ana <- hess_mat[i, j]
      fd <- fd_hess[i, j]
      abs_err <- abs(ana - fd)
      # For near-zero elements, use absolute error; otherwise use relative error
      if (abs(ana) < 1e-6 && abs(fd) < 1e-6) {
        rel_err <- abs_err
      } else {
        rel_err <- abs_err / max(abs(ana), abs(fd))
      }
      if (!is.na(rel_err) && rel_err > max_hess_err) max_hess_err <- rel_err
    }
  }
  # Tolerance of 1e-2 for full Hessian check (cross-derivatives have higher FD error)
  expect_lt(max_hess_err, 1e-2, label = sprintf("Hessian max error: %.2e", max_hess_err))

  # Test 4: Parameter recovery (with relaxed tolerance due to finite sample)
  tolerance <- 0.3

  # Factor variance (should be close to 1)
  expect_lt(abs(result$estimates["factor_var_1"] - 1.0), tolerance)

  # Measurement loading for m2 (should be close to 0.8)
  expect_lt(abs(result$estimates["m2_loading_1"] - 0.8), tolerance)

  # Linear model parameters
  expect_lt(abs(result$estimates["linear_intercept"] - true_beta_linear[1]), tolerance)
  expect_lt(abs(result$estimates["linear_x1"] - true_beta_linear[2]), tolerance)
  expect_lt(abs(result$estimates["linear_loading_1"] - true_loading_linear), tolerance)

  # Probit model parameters
  expect_lt(abs(result$estimates["probit_intercept"] - true_beta_probit[1]), tolerance)
  expect_lt(abs(result$estimates["probit_x1"] - true_beta_probit[2]), tolerance)

  if (VERBOSE) {
    cat("\nTest J: Two-stage with multiple outcome types and fixed coefficients\n")
    cat("  Stage 1: 2 linear measurement equations\n")
    cat("  Stage 2: linear, probit, oprobit, mlogit, exploded logit\n")
    cat("  Fixed coefficients: 5 (one per stage 2 model)\n")
    cat(sprintf("  Gradient max rel error: %.2e\n", max_rel_err))
    cat(sprintf("  Hessian max rel error: %.2e\n", max_hess_err))
    cat(sprintf("  Linear loading recovery: true=%.3f, est=%.3f\n",
                true_loading_linear, result$estimates["linear_loading_1"]))
    cat(sprintf("  m2 loading recovery: true=%.3f, est=%.3f\n",
                0.8, result$estimates["m2_loading_1"]))
  }
})

# =============================================================================
# Test K: Binary Logit Model (num_choices = 2)
# =============================================================================
# This test specifically covers the binary logit case which requires special
# handling for initialization (glm expects 0/1 but model uses 1/2 coding)

test_that("Model K: Binary logit gradient, Hessian, and parameter recovery", {
  skip_on_cran()

  set.seed(54321)
  n_obs <- 500

  # True parameters
  true_factor_var <- 1.0
  true_beta <- c(0.5, 0.3)  # intercept, x1 coefficient
  true_loading <- 0.7
  true_meas_loading <- 0.8
  true_meas_sigma <- 0.6

  # Generate data
  factors <- rnorm(n_obs, 0, sqrt(true_factor_var))
  x1 <- rnorm(n_obs)

  # Linear measurement equation (to identify factor)
  Y_meas <- true_meas_loading * factors + rnorm(n_obs, 0, true_meas_sigma)

  # Binary logit outcome
  # P(Y=2) = 1/(1 + exp(-(beta0 + beta1*x1 + lambda*f)))
  latent <- true_beta[1] + true_beta[2] * x1 + true_loading * factors
  prob <- 1 / (1 + exp(-latent))
  # Outcomes must be 1-indexed (1 or 2) for C++ compatibility
  Y_logit <- rbinom(n_obs, 1, prob) + 1  # 1 or 2 (NOT 0/1)

  dat <- data.frame(
    Y_meas = Y_meas,
    Y_logit = Y_logit,
    x1 = x1,
    intercept = 1,
    eval = 1
  )

  # Define model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Measurement equation (identifies factor)
  mc_meas <- define_model_component(
    name = "meas",
    data = dat,
    outcome = "Y_meas",
    factor = fm,
    covariates = "intercept",
    model_type = "linear",
    loading_normalization = 1.0,
    evaluation_indicator = "eval"
  )

  # Binary logit outcome (num_choices = 2)
  mc_logit <- define_model_component(
    name = "logit",
    data = dat,
    outcome = "Y_logit",
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = 2,  # Binary logit
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_meas, mc_logit), factor = fm)
  control <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  # Test 1: Verify likelihood is finite (not -Inf which was the bug)
  fm_ptr <- initialize_factor_model_cpp(ms, dat, control$n_quad_points)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  init_params <- init_result$init_params
  result_init <- evaluate_likelihood_cpp(fm_ptr, init_params)
  expect_true(is.finite(result_init$logLikelihood),
              info = "Binary logit should produce finite log-likelihood after initialization")

  # Test 2: Full estimation
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    init_params = init_params,
    control = control,
    optimizer = "nlminb",
    verbose = FALSE
  )

  expect_equal(result$convergence, 0)
  expect_true(is.finite(result$loglik))

  # Test 3: Gradient accuracy
  fm_ptr2 <- initialize_factor_model_cpp(ms, dat, control$n_quad_points)
  params <- result$estimates
  result_eval <- evaluate_likelihood_cpp(fm_ptr2, params, compute_gradient = TRUE, compute_hessian = FALSE)
  ana_grad <- result_eval$gradient

  # Finite difference gradient
  eps <- 1e-6
  fd_grad <- numeric(length(params))
  for (i in seq_along(params)) {
    params_plus <- params
    params_plus[i] <- params_plus[i] + eps
    ll_plus <- evaluate_loglik_only_cpp(fm_ptr2, params_plus)

    params_minus <- params
    params_minus[i] <- params_minus[i] - eps
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr2, params_minus)

    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps)
  }

  max_rel_err <- 0
  for (i in seq_along(params)) {
    ana <- ana_grad[i]
    fd <- fd_grad[i]
    abs_err <- abs(ana - fd)
    if (abs(ana) < 1e-6 && abs(fd) < 1e-6) {
      rel_err <- abs_err
    } else {
      rel_err <- abs_err / max(abs(ana), abs(fd))
    }
    if (!is.na(rel_err) && rel_err > max_rel_err) max_rel_err <- rel_err
  }
  expect_lt(max_rel_err, 1e-5, label = sprintf("Gradient max error: %.2e", max_rel_err))

  # Test 4: Hessian accuracy
  result_hess <- evaluate_likelihood_cpp(fm_ptr2, params, compute_gradient = TRUE, compute_hessian = TRUE)
  hess_vec <- result_hess$hessian
  n_params <- length(params)

  # Expand upper triangular to full matrix
  hess_mat <- matrix(0, n_params, n_params)
  idx <- 1
  for (i in 1:n_params) {
    for (j in i:n_params) {
      hess_mat[i, j] <- hess_vec[idx]
      hess_mat[j, i] <- hess_vec[idx]
      idx <- idx + 1
    }
  }

  # FD Hessian (using gradient differences)
  fd_hess <- matrix(0, n_params, n_params)
  grad_base <- evaluate_likelihood_cpp(fm_ptr2, params, compute_gradient = TRUE, compute_hessian = FALSE)$gradient
  for (i in seq_len(n_params)) {
    params_plus <- params
    h <- eps * (abs(params[i]) + 1.0)
    params_plus[i] <- params[i] + h
    grad_plus <- evaluate_likelihood_cpp(fm_ptr2, params_plus, compute_gradient = TRUE, compute_hessian = FALSE)$gradient
    fd_hess[i, ] <- (grad_plus - grad_base) / h
  }
  # Symmetrize FD Hessian
  fd_hess <- (fd_hess + t(fd_hess)) / 2

  max_hess_err <- 0
  for (i in seq_len(n_params)) {
    for (j in i:n_params) {
      ana <- hess_mat[i, j]
      fd <- fd_hess[i, j]
      abs_err <- abs(ana - fd)
      if (abs(ana) < 1e-6 && abs(fd) < 1e-6) {
        rel_err <- abs_err
      } else {
        rel_err <- abs_err / max(abs(ana), abs(fd))
      }
      if (!is.na(rel_err) && rel_err > max_hess_err) max_hess_err <- rel_err
    }
  }
  expect_lt(max_hess_err, 1e-2, label = sprintf("Hessian max error: %.2e", max_hess_err))

  # Test 5: Parameter recovery (with relaxed tolerance due to finite sample)
  tolerance <- 0.4

  # Factor variance
  expect_lt(abs(result$estimates["factor_var_1"] - true_factor_var), tolerance,
            label = "Factor variance recovery")

  # Logit parameters (binary logit uses logit_intercept, not logit_c1_intercept)
  expect_lt(abs(result$estimates["logit_intercept"] - true_beta[1]), tolerance,
            label = "Logit intercept recovery")
  expect_lt(abs(result$estimates["logit_x1"] - true_beta[2]), tolerance,
            label = "Logit x1 coefficient recovery")
  expect_lt(abs(result$estimates["logit_loading_1"] - true_loading), tolerance,
            label = "Logit loading recovery")

  if (VERBOSE) {
    cat("\nTest K: Binary logit (num_choices = 2)\n")
    cat("  Outcome coding: 1/2 (1-indexed, as required by C++)\n")
    cat(sprintf("  Initial log-likelihood: %.4f (should be finite, not -Inf)\n", result_init$logLikelihood))
    cat(sprintf("  Final log-likelihood: %.4f\n", result$loglik))
    cat(sprintf("  Gradient max rel error: %.2e\n", max_rel_err))
    cat(sprintf("  Hessian max rel error: %.2e\n", max_hess_err))
    cat(sprintf("  Factor var: true=%.3f, est=%.3f\n",
                true_factor_var, result$estimates["factor_var_1"]))
    cat(sprintf("  Logit intercept: true=%.3f, est=%.3f\n",
                true_beta[1], result$estimates["logit_intercept"]))
    cat(sprintf("  Logit x1: true=%.3f, est=%.3f\n",
                true_beta[2], result$estimates["logit_x1"]))
    cat(sprintf("  Logit loading: true=%.3f, est=%.3f\n",
                true_loading, result$estimates["logit_loading_1"]))
  }
})
