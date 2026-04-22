# Tests for Quadratic Factor Terms (factor_spec = "quadratic")
#
# This test suite validates:
# 1. Parameter initialization with quadratic loadings
# 2. Analytical gradients vs finite differences for quadratic terms
# 3. Analytical Hessians vs finite differences for quadratic terms
# 4. Model estimation with quadratic factor terms

# Test configuration
VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
GRAD_TOL <- 1e-3
HESS_TOL <- 1e-3

# ==============================================================================
# Test 1: Basic parameter structure with factor_spec = "quadratic"
# ==============================================================================

test_that("factor_spec='quadratic' creates correct parameter structure", {
  skip_on_cran()

  set.seed(201)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Create model with quadratic factor spec
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "test", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = NA_real_,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  # Check component structure
  expect_equal(mc$factor_spec, "quadratic")
  expect_equal(mc$n_quadratic_loadings, 1)  # 1 factor -> 1 quadratic loading
  expect_equal(mc$n_interaction_loadings, 0)  # No interactions for single factor

  # Create model system and initialize parameters
  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check that quadratic loading parameter exists
  expect_true("test_loading_quad_1" %in% init$param_names)

  # Verify parameter count: factor_var + intercept + x + linear_loading + quad_loading + sigma
  # But linear loading is free (NA_real_), so: 1 + 2 + 1 + 1 + 1 = 6
  # Wait - loading_normalization=NA means loading is estimated, but for single factor,
  # factor variance is fixed to 1.0, so:
  # factor_var (fixed but present) + intercept + x + loading + quad_loading + sigma
  n_expected <- length(init$param_names)
  expect_gte(n_expected, 5)  # At minimum: intercept, x, loading, quad_loading, sigma
})

# ==============================================================================
# Test 2: factor_spec='quadratic' for two factors
# ==============================================================================

test_that("factor_spec='quadratic' with multiple factors", {
  skip_on_cran()

  set.seed(202)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Create 2-factor model
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "test", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_),  # Fix first, estimate second
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  # Check component structure
  expect_equal(mc$factor_spec, "quadratic")
  expect_equal(mc$n_quadratic_loadings, 2)  # 2 factors -> 2 quadratic loadings
  expect_equal(mc$n_interaction_loadings, 0)  # No interactions (factor_spec != "interactions" or "full")

  # Create model system and initialize
  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check quadratic loading parameters exist
  expect_true("test_loading_quad_1" %in% init$param_names)
  expect_true("test_loading_quad_2" %in% init$param_names)
})

# ==============================================================================
# Test 3: Gradient accuracy for linear model with quadratic factor terms
# ==============================================================================

test_that("Gradient accuracy for linear model with quadratic factor terms", {
  skip_on_cran()

  set.seed(203)

  # Generate data with quadratic factor effect
  n <- 300
  f <- rnorm(n)  # Latent factor

  # True model: Y = intercept + x*beta + lambda*f + lambda_quad*f^2 + error
  true_intercept <- 2.0
  true_beta <- 0.5
  true_lambda <- 1.0  # Fixed for identification
  true_lambda_quad <- 0.3
  true_sigma <- 0.5

  x <- rnorm(n)
  Y <- true_intercept + true_beta * x + true_lambda * f + true_lambda_quad * f^2 + rnorm(n, 0, true_sigma)

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = 1.0,  # Fix linear loading for identification
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)

  # Initialize
  init <- initialize_parameters(ms, dat)

  # Test parameters close to true values
  # Parameter order: factor_var, intercept, x, quad_loading, sigma
  # (linear loading is fixed to 1.0, so not in free params)
  test_params <- init$init_params  # Use initialized values

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, test_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for quadratic linear model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 4: Hessian accuracy for linear model with quadratic factor terms
# ==============================================================================

test_that("Hessian accuracy for linear model with quadratic factor terms", {
  skip_on_cran()

  set.seed(204)

  # Generate data
  n <- 300
  f <- rnorm(n)

  true_intercept <- 2.0
  true_beta <- 0.5
  true_lambda <- 1.0
  true_lambda_quad <- 0.3
  true_sigma <- 0.5

  x <- rnorm(n)
  Y <- true_intercept + true_beta * x + true_lambda * f + true_lambda_quad * f^2 + rnorm(n, 0, true_sigma)

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run Hessian check
  hess_check <- check_hessian_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = HESS_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(hess_check$pass,
              info = sprintf("Hessian check failed for quadratic linear model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 5: Gradient accuracy for probit model with quadratic factor terms
# ==============================================================================

test_that("Gradient accuracy for probit model with quadratic factor terms", {
  skip_on_cran()

  set.seed(205)

  # Generate data
  n <- 300
  f <- rnorm(n)

  true_intercept <- 0.5
  true_beta <- 0.7
  true_lambda <- 1.0  # Fixed
  true_lambda_quad <- 0.4

  x <- rnorm(n)
  z <- true_intercept + true_beta * x + true_lambda * f + true_lambda_quad * f^2
  Y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for quadratic probit model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 6: Gradient accuracy for ordered probit with quadratic factor terms
# ==============================================================================

test_that("Gradient accuracy for ordered probit with quadratic factor terms", {
  skip_on_cran()

  set.seed(206)

  # Generate data
  n <- 400
  f <- rnorm(n)

  true_beta <- 0.6
  true_lambda <- 1.0  # Fixed
  true_lambda_quad <- 0.3
  thresh1 <- -0.5
  thresh2 <- 0.5

  x <- rnorm(n)
  z <- true_beta * x + true_lambda * f + true_lambda_quad * f^2
  u <- rnorm(n)
  Y <- ifelse(z + u < thresh1, 1, ifelse(z + u < thresh2, 2, 3))

  dat <- data.frame(x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = "x",  # No intercept for oprobit
    model_type = "oprobit",
    num_choices = 3,
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval",
    intercept = FALSE
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for quadratic oprobit model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 7: Gradient accuracy for logit model with quadratic factor terms
# ==============================================================================

test_that("Gradient accuracy for logit model with quadratic factor terms", {
  skip_on_cran()

  set.seed(207)

  # Generate data for multinomial logit (3 choices)
  n <- 300
  f <- rnorm(n)
  x <- rnorm(n)

  # True parameters for choices 2 and 3 (choice 1 is reference)
  true_int1 <- 0.3
  true_int2 <- -0.2
  true_beta1 <- 0.5
  true_beta2 <- 0.3
  true_lambda1 <- 1.0  # Fixed
  true_lambda2 <- 0.8
  true_lambda_quad1 <- 0.2
  true_lambda_quad2 <- 0.15

  # Compute utilities for each choice
  v1 <- rep(0, n)  # Reference category
  v2 <- true_int1 + true_beta1 * x + true_lambda1 * f + true_lambda_quad1 * f^2
  v3 <- true_int2 + true_beta2 * x + true_lambda2 * f + true_lambda_quad2 * f^2

  # Multinomial logit probabilities
  exp_v <- cbind(exp(v1), exp(v2), exp(v3))
  probs <- exp_v / rowSums(exp_v)

  # Sample outcomes
  Y <- integer(n)
  for (i in 1:n) {
    Y[i] <- sample(1:3, 1, prob = probs[i, ])
  }

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "logit",
    num_choices = 3,
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for quadratic logit model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 8: Hessian accuracy for probit model with quadratic factor terms
# ==============================================================================

test_that("Hessian accuracy for probit model with quadratic factor terms", {
  skip_on_cran()

  set.seed(211)

  # Generate data
  n <- 300
  f <- rnorm(n)

  true_intercept <- 0.5
  true_beta <- 0.7
  true_lambda <- 1.0  # Fixed
  true_lambda_quad <- 0.4

  x <- rnorm(n)
  z <- true_intercept + true_beta * x + true_lambda * f + true_lambda_quad * f^2
  Y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run Hessian check
  hess_check <- check_hessian_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = HESS_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(hess_check$pass,
              info = sprintf("Hessian check failed for quadratic probit model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 9: Hessian accuracy for ordered probit with quadratic factor terms
# ==============================================================================

test_that("Hessian accuracy for ordered probit with quadratic factor terms", {
  skip_on_cran()

  set.seed(212)

  # Generate data
  n <- 400
  f <- rnorm(n)

  true_beta <- 0.6
  true_lambda <- 1.0  # Fixed
  true_lambda_quad <- 0.3
  thresh1 <- -0.5
  thresh2 <- 0.5

  x <- rnorm(n)
  z <- true_beta * x + true_lambda * f + true_lambda_quad * f^2
  u <- rnorm(n)
  Y <- ifelse(z + u < thresh1, 1, ifelse(z + u < thresh2, 2, 3))

  dat <- data.frame(x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = "x",  # No intercept for oprobit
    model_type = "oprobit",
    num_choices = 3,
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval",
    intercept = FALSE
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run Hessian check
  hess_check <- check_hessian_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = HESS_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(hess_check$pass,
              info = sprintf("Hessian check failed for quadratic oprobit model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 10: Hessian accuracy for logit model with quadratic factor terms
# ==============================================================================

test_that("Hessian accuracy for logit model with quadratic factor terms", {
  skip_on_cran()

  set.seed(213)

  # Generate data for multinomial logit (3 choices)
  n <- 300
  f <- rnorm(n)
  x <- rnorm(n)

  # True parameters for choices 2 and 3 (choice 1 is reference)
  true_int1 <- 0.3
  true_int2 <- -0.2
  true_beta1 <- 0.5
  true_beta2 <- 0.3
  true_lambda1 <- 1.0  # Fixed
  true_lambda2 <- 0.8
  true_lambda_quad1 <- 0.2
  true_lambda_quad2 <- 0.15

  # Compute utilities for each choice
  v1 <- rep(0, n)  # Reference category
  v2 <- true_int1 + true_beta1 * x + true_lambda1 * f + true_lambda_quad1 * f^2
  v3 <- true_int2 + true_beta2 * x + true_lambda2 * f + true_lambda_quad2 * f^2

  # Multinomial logit probabilities
  exp_v <- cbind(exp(v1), exp(v2), exp(v3))
  probs <- exp_v / rowSums(exp_v)

  # Sample outcomes
  Y <- integer(n)
  for (i in 1:n) {
    Y[i] <- sample(1:3, 1, prob = probs[i, ])
  }

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "logit",
    num_choices = 3,
    loading_normalization = 1.0,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Run Hessian check
  hess_check <- check_hessian_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = HESS_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(hess_check$pass,
              info = sprintf("Hessian check failed for quadratic logit model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 11: Model estimation with quadratic factor terms (linear)
# ==============================================================================

test_that("Model estimation with quadratic factor terms recovers parameters", {
  skip_on_cran()

  set.seed(208)

  # Generate data with known parameters
  n <- 1000
  f <- rnorm(n)

  # True parameters
  true_factor_var <- 1.0  # Fixed for identification
  true_intercept <- 2.0
  true_beta <- 0.5
  true_lambda <- 1.0  # Fixed for identification
  true_lambda_quad <- 0.3
  true_sigma <- 0.5

  x <- rnorm(n)
  Y <- true_intercept + true_beta * x + true_lambda * f + true_lambda_quad * f^2 + rnorm(n, 0, true_sigma)

  # Add measurement equations to identify the factor
  T1 <- 1.0 + 1.0 * f + rnorm(n, 0, 0.3)  # Loading fixed to 1.0
  T2 <- 1.5 + 1.2 * f + rnorm(n, 0, 0.4)

  dat <- data.frame(
    intercept = 1, x = x, Y = Y, T1 = T1, T2 = T2, eval = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Measurement equations (standard linear)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0,  # Fix for identification
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  # Outcome with quadratic factor effect
  mc_Y <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = NA_real_,  # Estimate this one
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_Y), factor = fm)

  # Estimate model
  control <- define_estimation_control(num_cores = 1)
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = control,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0, info = "Model with quadratic factors failed to converge")

  # Check that quadratic loading was estimated
  expect_true("Y_loading_quad_1" %in% names(result$estimates),
              info = "Quadratic loading parameter not found in estimates")

  # Check reasonable estimates
  est_quad <- result$estimates["Y_loading_quad_1"]
  expect_true(abs(est_quad - true_lambda_quad) < 0.3,
              info = sprintf("Quadratic loading estimate (%.3f) far from true value (%.3f)",
                            est_quad, true_lambda_quad))
})

# ==============================================================================
# Test 12: Multiple model types with quadratic factors in same system
# ==============================================================================

test_that("Multiple model types with quadratic factors work together", {
  skip_on_cran()

  set.seed(209)

  n <- 500
  f <- rnorm(n)
  x <- rnorm(n)

  # Linear outcome with quadratic factor effect
  Y_lin <- 2.0 + 0.5 * x + 1.0 * f + 0.3 * f^2 + rnorm(n, 0, 0.5)

  # Probit outcome with quadratic factor effect
  z_prob <- 0.5 + 0.7 * x + 0.8 * f + 0.2 * f^2
  Y_prob <- as.numeric(runif(n) < pnorm(z_prob))

  # Standard measurement equation (no quadratic)
  T1 <- 1.0 + 1.0 * f + rnorm(n, 0, 0.3)

  dat <- data.frame(
    intercept = 1, x = x,
    Y_lin = Y_lin, Y_prob = Y_prob, T1 = T1,
    eval = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Standard measurement (for identification)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0,
    evaluation_indicator = "eval"
  )

  # Linear with quadratic
  mc_lin <- define_model_component(
    name = "Y_lin", data = dat, outcome = "Y_lin", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = NA_real_,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  # Probit with quadratic
  mc_prob <- define_model_component(
    name = "Y_prob", data = dat, outcome = "Y_prob", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = NA_real_,
    factor_spec = "quadratic",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_lin, mc_prob), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check both quadratic parameters exist
  expect_true("Y_lin_loading_quad_1" %in% init$param_names)
  expect_true("Y_prob_loading_quad_1" %in% init$param_names)

  # Run gradient check to verify derivatives work across model types
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for mixed model types with quadratic factors (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 13: factor_spec="linear" (default) has no quadratic loadings
# ==============================================================================

test_that("factor_spec='linear' creates no quadratic loadings", {
  skip_on_cran()

  set.seed(210)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Default factor_spec should be "linear"
  mc <- define_model_component(
    name = "test", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  expect_equal(mc$factor_spec, "linear")
  expect_equal(mc$n_quadratic_loadings, 0)
  expect_equal(mc$n_interaction_loadings, 0)

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # No quadratic parameter should exist
  expect_false(any(grepl("loading_quad", init$param_names)))
})
