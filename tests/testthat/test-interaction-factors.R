# Tests for Interaction Factor Terms (factor_spec = "interactions" or "full")
#
# This test suite validates:
# 1. Parameter initialization with interaction loadings
# 2. Analytical gradients vs finite differences for interaction terms
# 3. Analytical Hessians vs finite differences for interaction terms
# 4. Model estimation with interaction factor terms
#
# Note: Interaction terms require at least 2 factors (f_j * f_k for j < k)

# Test configuration
VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
GRAD_TOL <- 1e-3
HESS_TOL <- 1e-3

# ==============================================================================
# Test 1: Basic parameter structure with factor_spec = "interactions"
# ==============================================================================

test_that("factor_spec='interactions' creates correct parameter structure", {
  skip_on_cran()

  set.seed(301)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Create 2-factor model (minimum for interactions)
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "test", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_),  # Fix first, estimate second
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  # Check component structure
  expect_equal(mc$factor_spec, "interactions")
  expect_equal(mc$n_quadratic_loadings, 0)  # No quadratics for "interactions"
  expect_equal(mc$n_interaction_loadings, 1)  # 2 factors -> 1 interaction (f1*f2)

  # Create model system and initialize parameters
  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check that interaction loading parameter exists
  expect_true("test_loading_inter_1_2" %in% init$param_names)

  # No quadratic parameters should exist
  expect_false(any(grepl("loading_quad", init$param_names)))
})

# ==============================================================================
# Test 2: factor_spec='full' for two factors (both quadratic and interactions)
# ==============================================================================

test_that("factor_spec='full' creates quadratic AND interaction parameters", {
  skip_on_cran()

  set.seed(302)
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
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval"
  )

  # Check component structure
  expect_equal(mc$factor_spec, "full")
  expect_equal(mc$n_quadratic_loadings, 2)  # 2 factors -> 2 quadratic loadings
  expect_equal(mc$n_interaction_loadings, 1)  # 2 factors -> 1 interaction

  # Create model system and initialize
  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check both quadratic and interaction parameters exist
  expect_true("test_loading_quad_1" %in% init$param_names)
  expect_true("test_loading_quad_2" %in% init$param_names)
  expect_true("test_loading_inter_1_2" %in% init$param_names)
})

# ==============================================================================
# Test 3: factor_spec='interactions' with 3 factors
# ==============================================================================

test_that("factor_spec='interactions' with 3 factors creates correct interactions", {
  skip_on_cran()

  set.seed(303)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Create 3-factor model
  fm <- define_factor_model(n_factors = 3, n_types = 1)

  mc <- define_model_component(
    name = "test", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_, NA_real_),
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  # Check component structure
  # 3 factors -> 3*(3-1)/2 = 3 interactions: f1*f2, f1*f3, f2*f3
  expect_equal(mc$n_interaction_loadings, 3)
  expect_equal(mc$n_quadratic_loadings, 0)

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check all interaction parameters exist
  expect_true("test_loading_inter_1_2" %in% init$param_names)
  expect_true("test_loading_inter_1_3" %in% init$param_names)
  expect_true("test_loading_inter_2_3" %in% init$param_names)
})

# ==============================================================================
# Test 4: factor_spec='interactions' downgrades for single factor
# ==============================================================================

test_that("factor_spec='interactions' downgrades to 'linear' for single factor", {
  skip_on_cran()

  set.seed(304)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Single factor model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Should warn and downgrade
  expect_warning(
    mc <- define_model_component(
      name = "test", data = dat, outcome = "Y", factor = fm,
      covariates = c("intercept", "x"), model_type = "linear",
      loading_normalization = NA_real_,
      factor_spec = "interactions",
      evaluation_indicator = "eval"
    ),
    "requires n_factors >= 2"
  )

  # Should have been downgraded to "linear"
  expect_equal(mc$factor_spec, "linear")
  expect_equal(mc$n_interaction_loadings, 0)
})

# ==============================================================================
# Test 5: factor_spec='full' downgrades to 'quadratic' for single factor
# ==============================================================================

test_that("factor_spec='full' downgrades to 'quadratic' for single factor", {
  skip_on_cran()

  set.seed(305)
  n <- 100
  dat <- data.frame(
    intercept = 1,
    x = rnorm(n),
    Y = rnorm(n),
    eval = 1
  )

  # Single factor model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Should warn and downgrade
  expect_warning(
    mc <- define_model_component(
      name = "test", data = dat, outcome = "Y", factor = fm,
      covariates = c("intercept", "x"), model_type = "linear",
      loading_normalization = NA_real_,
      factor_spec = "full",
      evaluation_indicator = "eval"
    ),
    "requires n_factors >= 2"
  )

  # Should have been downgraded to "quadratic" (keeps quadratic, drops interactions)
  expect_equal(mc$factor_spec, "quadratic")
  expect_equal(mc$n_quadratic_loadings, 1)
  expect_equal(mc$n_interaction_loadings, 0)
})

# ==============================================================================
# Test 6: Gradient accuracy for linear model with interaction factor terms
# ==============================================================================

test_that("Gradient accuracy for linear model with interaction factor terms", {
  skip_on_cran()

  set.seed(306)

  # Generate data with 2 factors and interaction effect
  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  # True model: Y = intercept + x*beta + lambda1*f1 + lambda2*f2 + lambda_inter*f1*f2 + error
  true_intercept <- 2.0
  true_beta <- 0.5
  true_lambda1 <- 1.0  # Fixed for identification
  true_lambda2 <- 0.8
  true_lambda_inter <- 0.3
  true_sigma <- 0.5

  x <- rnorm(n)
  Y <- true_intercept + true_beta * x + true_lambda1 * f1 + true_lambda2 * f2 +
       true_lambda_inter * f1 * f2 + rnorm(n, 0, true_sigma)

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  # Create model
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_),  # Fix first loading
    factor_spec = "interactions",
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
              info = sprintf("Gradient check failed for interaction linear model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 7: Hessian accuracy for linear model with interaction factor terms
# ==============================================================================

test_that("Hessian accuracy for linear model with interaction factor terms", {
  skip_on_cran()

  set.seed(307)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  Y <- 2.0 + 0.5 * x + 1.0 * f1 + 0.8 * f2 + 0.3 * f1 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Hessian check failed for interaction linear model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 8: Gradient accuracy for probit model with interaction factor terms
# ==============================================================================

test_that("Gradient accuracy for probit model with interaction factor terms", {
  skip_on_cran()

  set.seed(308)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  z <- 0.5 + 0.7 * x + 1.0 * f1 + 0.6 * f2 + 0.25 * f1 * f2
  Y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Gradient check failed for interaction probit model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 9: Hessian accuracy for probit model with interaction factor terms
# ==============================================================================

test_that("Hessian accuracy for probit model with interaction factor terms", {
  skip_on_cran()

  set.seed(309)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  z <- 0.5 + 0.7 * x + 1.0 * f1 + 0.6 * f2 + 0.25 * f1 * f2
  Y <- as.numeric(runif(n) < pnorm(z))

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Hessian check failed for interaction probit model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 10: Gradient accuracy for ordered probit with interaction factor terms
# ==============================================================================

test_that("Gradient accuracy for ordered probit with interaction factor terms", {
  skip_on_cran()

  set.seed(310)

  n <- 400
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  z <- 0.6 * x + 1.0 * f1 + 0.7 * f2 + 0.2 * f1 * f2
  u <- rnorm(n)
  Y <- ifelse(z + u < -0.5, 1, ifelse(z + u < 0.5, 2, 3))

  dat <- data.frame(x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = "x",
    model_type = "oprobit",
    num_choices = 3,
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Gradient check failed for interaction oprobit model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 11: Hessian accuracy for ordered probit with interaction factor terms
# ==============================================================================

test_that("Hessian accuracy for ordered probit with interaction factor terms", {
  skip_on_cran()

  set.seed(311)

  n <- 400
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  z <- 0.6 * x + 1.0 * f1 + 0.7 * f2 + 0.2 * f1 * f2
  u <- rnorm(n)
  Y <- ifelse(z + u < -0.5, 1, ifelse(z + u < 0.5, 2, 3))

  dat <- data.frame(x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = "x",
    model_type = "oprobit",
    num_choices = 3,
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Hessian check failed for interaction oprobit model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 12: Gradient accuracy for logit model with interaction factor terms
# ==============================================================================

test_that("Gradient accuracy for logit model with interaction factor terms", {
  skip_on_cran()

  set.seed(312)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  x <- rnorm(n)

  # Multinomial logit with 3 choices
  v1 <- rep(0, n)  # Reference
  v2 <- 0.3 + 0.5 * x + 1.0 * f1 + 0.6 * f2 + 0.2 * f1 * f2
  v3 <- -0.2 + 0.3 * x + 0.8 * f1 + 0.5 * f2 + 0.15 * f1 * f2

  exp_v <- cbind(exp(v1), exp(v2), exp(v3))
  probs <- exp_v / rowSums(exp_v)

  Y <- integer(n)
  for (i in 1:n) {
    Y[i] <- sample(1:3, 1, prob = probs[i, ])
  }

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "logit",
    num_choices = 3,
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Gradient check failed for interaction logit model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 13: Hessian accuracy for logit model with interaction factor terms
# ==============================================================================

test_that("Hessian accuracy for logit model with interaction factor terms", {
  skip_on_cran()

  set.seed(313)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  x <- rnorm(n)

  v1 <- rep(0, n)
  v2 <- 0.3 + 0.5 * x + 1.0 * f1 + 0.6 * f2 + 0.2 * f1 * f2
  v3 <- -0.2 + 0.3 * x + 0.8 * f1 + 0.5 * f2 + 0.15 * f1 * f2

  exp_v <- cbind(exp(v1), exp(v2), exp(v3))
  probs <- exp_v / rowSums(exp_v)

  Y <- integer(n)
  for (i in 1:n) {
    Y[i] <- sample(1:3, 1, prob = probs[i, ])
  }

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "logit",
    num_choices = 3,
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "interactions",
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
              info = sprintf("Hessian check failed for interaction logit model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 14: Gradient accuracy for "full" spec (quadratic + interactions)
# ==============================================================================

test_that("Gradient accuracy for 'full' spec (quadratic + interactions)", {
  skip_on_cran()

  set.seed(314)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  # Full model: linear + quadratic + interaction terms
  Y <- 2.0 + 0.5 * x + 1.0 * f1 + 0.8 * f2 +
       0.2 * f1^2 + 0.15 * f2^2 + 0.25 * f1 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check all second-order parameters exist
  expect_true("Y_loading_quad_1" %in% init$param_names)
  expect_true("Y_loading_quad_2" %in% init$param_names)
  expect_true("Y_loading_inter_1_2" %in% init$param_names)

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for 'full' spec model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 15: Hessian accuracy for "full" spec (quadratic + interactions)
# ==============================================================================

test_that("Hessian accuracy for 'full' spec (quadratic + interactions)", {
  skip_on_cran()

  set.seed(315)

  n <- 300
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  x <- rnorm(n)
  Y <- 2.0 + 0.5 * x + 1.0 * f1 + 0.8 * f2 +
       0.2 * f1^2 + 0.15 * f2^2 + 0.25 * f1 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept = 1, x = x, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(1.0, NA_real_),
    factor_spec = "full",
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
              info = sprintf("Hessian check failed for 'full' spec model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 16: Model estimation with interaction factor terms
# ==============================================================================

test_that("Model estimation with interaction factor terms recovers parameters", {
  skip_on_cran()

  set.seed(316)

  # Generate data with known parameters
  n <- 1000
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  # True parameters
  true_intercept <- 2.0
  true_beta <- 0.5
  true_lambda1 <- 1.0  # Fixed for identification
  true_lambda2 <- 0.8
  true_lambda_inter <- 0.3
  true_sigma <- 0.5

  x <- rnorm(n)
  Y <- true_intercept + true_beta * x + true_lambda1 * f1 + true_lambda2 * f2 +
       true_lambda_inter * f1 * f2 + rnorm(n, 0, true_sigma)

  # Add measurement equations to identify both factors
  T1 <- 1.0 + 1.0 * f1 + rnorm(n, 0, 0.3)  # Loading fixed to 1.0 for factor 1
  T2 <- 1.5 + 0.9 * f1 + rnorm(n, 0, 0.4)
  T3 <- 0.5 + 1.0 * f2 + rnorm(n, 0, 0.3)  # Loading fixed to 1.0 for factor 2
  T4 <- 1.2 + 1.1 * f2 + rnorm(n, 0, 0.4)

  dat <- data.frame(
    intercept = 1, x = x, Y = Y,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4,
    eval = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Measurement equations for factor 1
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0),  # Fixed for identification
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0),
    evaluation_indicator = "eval"
  )

  # Measurement equations for factor 2
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0),  # Fixed for identification
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_),
    evaluation_indicator = "eval"
  )

  # Outcome with interaction effect
  mc_Y <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_Y),
    factor = fm
  )

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
  expect_equal(result$convergence, 0, info = "Model with interaction factors failed to converge")

  # Check that interaction loading was estimated
  expect_true("Y_loading_inter_1_2" %in% names(result$estimates),
              info = "Interaction loading parameter not found in estimates")

  # Check reasonable estimate (allowing for wider tolerance in complex model)
  est_inter <- result$estimates["Y_loading_inter_1_2"]
  expect_true(abs(est_inter - true_lambda_inter) < 0.4,
              info = sprintf("Interaction loading estimate (%.3f) far from true value (%.3f)",
                            est_inter, true_lambda_inter))
})

# ==============================================================================
# Test 17: Mixed components with different factor_specs
# ==============================================================================

test_that("Mixed components with different factor_specs work together", {
  skip_on_cran()

  set.seed(317)

  n <- 500
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  x <- rnorm(n)

  # Linear outcome with interaction only
  Y_lin <- 2.0 + 0.5 * x + 1.0 * f1 + 0.8 * f2 + 0.3 * f1 * f2 + rnorm(n, 0, 0.5)

  # Probit outcome with full spec (quadratic + interaction)
  z_prob <- 0.5 + 0.7 * x + 0.6 * f1 + 0.5 * f2 + 0.1 * f1^2 + 0.1 * f2^2 + 0.15 * f1 * f2
  Y_prob <- as.numeric(runif(n) < pnorm(z_prob))

  # Standard linear measurement (no second-order terms)
  T1 <- 1.0 + 1.0 * f1 + rnorm(n, 0, 0.3)
  T2 <- 0.8 + 1.0 * f2 + rnorm(n, 0, 0.3)

  dat <- data.frame(
    intercept = 1, x = x,
    Y_lin = Y_lin, Y_prob = Y_prob, T1 = T1, T2 = T2,
    eval = 1
  )

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Standard measurements (linear factor_spec)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0),
    evaluation_indicator = "eval"
  )

  # Linear with interactions only
  mc_lin <- define_model_component(
    name = "Y_lin", data = dat, outcome = "Y_lin", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  # Probit with full spec
  mc_prob <- define_model_component(
    name = "Y_prob", data = dat, outcome = "Y_prob", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc_T1, mc_T2, mc_lin, mc_prob), factor = fm)
  init <- initialize_parameters(ms, dat)

  # Check appropriate parameters exist for each component
  expect_true("Y_lin_loading_inter_1_2" %in% init$param_names)
  expect_false(any(grepl("Y_lin_loading_quad", init$param_names)))

  expect_true("Y_prob_loading_quad_1" %in% init$param_names)
  expect_true("Y_prob_loading_quad_2" %in% init$param_names)
  expect_true("Y_prob_loading_inter_1_2" %in% init$param_names)

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 6
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for mixed factor_specs (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 18: Probit estimation with interaction factor terms
# ==============================================================================

test_that("Probit model estimation with interaction factor terms recovers parameters", {
  skip_on_cran()

  set.seed(318)

  # Generate data with known parameters
  n <- 1500
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  # True parameters for probit outcome
  true_intercept <- 0.5
  true_beta <- 0.7
  true_lambda1 <- 0.8
  true_lambda2 <- 0.6
  true_lambda_inter <- 0.4

  x <- rnorm(n)
  z <- true_intercept + true_beta * x + true_lambda1 * f1 + true_lambda2 * f2 +
       true_lambda_inter * f1 * f2
  Y_prob <- as.numeric(runif(n) < pnorm(z))

  # Linear measurement equations to identify both factors
  T1 <- 1.0 + 1.0 * f1 + rnorm(n, 0, 0.3)
  T2 <- 1.5 + 0.9 * f1 + rnorm(n, 0, 0.4)
  T3 <- 0.5 + 1.0 * f2 + rnorm(n, 0, 0.3)
  T4 <- 1.2 + 1.1 * f2 + rnorm(n, 0, 0.4)

  dat <- data.frame(
    intercept = 1, x = x, Y_prob = Y_prob,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4,
    eval = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Measurement equations for factor 1 (linear)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0),
    evaluation_indicator = "eval"
  )

  # Measurement equations for factor 2 (linear)
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0),
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_),
    evaluation_indicator = "eval"
  )

  # Probit outcome with interaction effect
  mc_Y <- define_model_component(
    name = "Y_prob", data = dat, outcome = "Y_prob", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_Y),
    factor = fm
  )

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
  expect_equal(result$convergence, 0, info = "Probit model with interaction factors failed to converge")

  # Check that interaction loading was estimated
  expect_true("Y_prob_loading_inter_1_2" %in% names(result$estimates),
              info = "Interaction loading parameter not found in probit estimates")

  # Check reasonable estimate (wider tolerance for probit)
  est_inter <- result$estimates["Y_prob_loading_inter_1_2"]
  expect_true(abs(est_inter - true_lambda_inter) < 0.5,
              info = sprintf("Probit interaction loading estimate (%.3f) far from true value (%.3f)",
                            est_inter, true_lambda_inter))
})

# ==============================================================================
# Test 19: Logit estimation with interaction factor terms
# ==============================================================================

test_that("Logit model estimation with interaction factor terms recovers parameters", {
  skip_on_cran()

  set.seed(319)

  # Generate data with known parameters
  n <- 1500
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  # True parameters for multinomial logit outcome (3 choices)
  # Choice 1 is reference, we model z1 (choice 2 vs 1) and z2 (choice 3 vs 1)
  true_int1 <- 0.3
  true_beta1 <- 0.6
  true_lambda1_ch1 <- 0.7
  true_lambda2_ch1 <- 0.5
  true_lambda_inter_ch1 <- 0.35  # Interaction for choice 1 (first non-reference)

  true_int2 <- -0.2
  true_beta2 <- 0.4
  true_lambda1_ch2 <- 0.5
  true_lambda2_ch2 <- 0.3
  true_lambda_inter_ch2 <- 0.25

  x <- rnorm(n)
  # Multinomial logit with 3 choices (reference = choice 1)
  z1 <- true_int1 + true_beta1 * x + true_lambda1_ch1 * f1 + true_lambda2_ch1 * f2 +
        true_lambda_inter_ch1 * f1 * f2
  z2 <- true_int2 + true_beta2 * x + true_lambda1_ch2 * f1 + true_lambda2_ch2 * f2 +
        true_lambda_inter_ch2 * f1 * f2

  exp_z0 <- 1
  exp_z1 <- exp(z1)
  exp_z2 <- exp(z2)
  denom <- exp_z0 + exp_z1 + exp_z2

  p0 <- exp_z0 / denom
  p1 <- exp_z1 / denom
  p2 <- exp_z2 / denom

  # C++ expects multinomial choices coded as 1, 2, 3
  Y_logit <- numeric(n)
  for (i in seq_len(n)) {
    Y_logit[i] <- sample(1:3, 1, prob = c(p0[i], p1[i], p2[i]))
  }

  # Linear measurement equations to identify both factors
  T1 <- 1.0 + 1.0 * f1 + rnorm(n, 0, 0.3)
  T2 <- 1.5 + 0.9 * f1 + rnorm(n, 0, 0.4)
  T3 <- 0.5 + 1.0 * f2 + rnorm(n, 0, 0.3)
  T4 <- 1.2 + 1.1 * f2 + rnorm(n, 0, 0.4)

  dat <- data.frame(
    intercept = 1, x = x, Y_logit = Y_logit,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4,
    eval = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Measurement equations for factor 1 (linear)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0),
    evaluation_indicator = "eval"
  )

  # Measurement equations for factor 2 (linear)
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0),
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_),
    evaluation_indicator = "eval"
  )

  # Logit outcome with interaction effect (multinomial with 3 choices)
  mc_Y <- define_model_component(
    name = "Y_logit", data = dat, outcome = "Y_logit", factor = fm,
    covariates = c("intercept", "x"), model_type = "logit",
    num_choices = 3,  # 3 choices: 1, 2, 3 (reference = 1)
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "interactions",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_Y),
    factor = fm
  )

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
  expect_equal(result$convergence, 0, info = "Logit model with interaction factors failed to converge")

  # Check that interaction loading was estimated (for first non-reference choice)
  # Parameter naming: Y_logit_c1_loading_inter_1_2 for choice 1
  expect_true("Y_logit_c1_loading_inter_1_2" %in% names(result$estimates),
              info = "Interaction loading parameter not found in logit estimates")

  # Check reasonable estimate (wider tolerance for multinomial logit)
  est_inter <- result$estimates["Y_logit_c1_loading_inter_1_2"]
  expect_true(abs(est_inter - true_lambda_inter_ch1) < 0.6,
              info = sprintf("Logit interaction loading estimate (%.3f) far from true value (%.3f)",
                            est_inter, true_lambda_inter_ch1))
})

# ==============================================================================
# Test 20: Ordered probit estimation with interaction factor terms
# ==============================================================================

test_that("Ordered probit model estimation with interaction factor terms recovers parameters", {
  skip_on_cran()

  set.seed(320)

  # Generate data with known parameters
  n <- 1500
  f1 <- rnorm(n)
  f2 <- rnorm(n)

  # True parameters for ordered probit outcome
  true_intercept <- 0.0  # Usually fixed to 0 for identification
  true_beta <- 0.5
  true_lambda1 <- 0.6
  true_lambda2 <- 0.5
  true_lambda_inter <- 0.3
  # Thresholds (3 categories = 2 thresholds)
  true_thresh <- c(-0.5, 0.8)

  x <- rnorm(n)
  # Ordered probit: latent variable
  z_star <- true_intercept + true_beta * x + true_lambda1 * f1 + true_lambda2 * f2 +
            true_lambda_inter * f1 * f2 + rnorm(n)

  # Categorize based on thresholds, convert to integer for C++ compatibility
  Y_oprobit_fac <- cut(z_star,
                       breaks = c(-Inf, true_thresh, Inf),
                       ordered_result = TRUE)
  Y_oprobit <- as.integer(Y_oprobit_fac)  # Convert to 1, 2, 3

  # Linear measurement equations to identify both factors
  T1 <- 1.0 + 1.0 * f1 + rnorm(n, 0, 0.3)
  T2 <- 1.5 + 0.9 * f1 + rnorm(n, 0, 0.4)
  T3 <- 0.5 + 1.0 * f2 + rnorm(n, 0, 0.3)
  T4 <- 1.2 + 1.1 * f2 + rnorm(n, 0, 0.4)

  dat <- data.frame(
    intercept = 1, x = x, Y_oprobit = Y_oprobit,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4,
    eval = 1
  )

  # Create model system
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Measurement equations for factor 1 (linear)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0),
    evaluation_indicator = "eval"
  )

  # Measurement equations for factor 2 (linear)
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0),
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_),
    evaluation_indicator = "eval"
  )

  # Ordered probit outcome with interaction effect.
  # Note: oprobit has no intercept covariate — the intercept is absorbed into
  # the cut points.
  mc_Y <- define_model_component(
    name = "Y_oprobit", data = dat, outcome = "Y_oprobit", factor = fm,
    covariates = "x", model_type = "oprobit",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "interactions",
    evaluation_indicator = "eval",
    num_choices = 3
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_Y),
    factor = fm
  )

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
  expect_equal(result$convergence, 0, info = "Ordered probit model with interaction factors failed to converge")

  # Check that interaction loading was estimated
  expect_true("Y_oprobit_loading_inter_1_2" %in% names(result$estimates),
              info = "Interaction loading parameter not found in ordered probit estimates")

  # Check reasonable estimate (wider tolerance for oprobit)
  est_inter <- result$estimates["Y_oprobit_loading_inter_1_2"]
  expect_true(abs(est_inter - true_lambda_inter) < 0.5,
              info = sprintf("Ordered probit interaction loading estimate (%.3f) far from true value (%.3f)",
                            est_inter, true_lambda_inter))
})

# ==============================================================================
# Test 21: Gradient accuracy for 3-factor model with full second order terms
# ==============================================================================

test_that("Gradient accuracy for 3-factor linear model with full second order terms", {
  # This test runs by default (moderate computational cost)
  set.seed(321)

  # Generate data with 3 factors
  n <- 400
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  f3 <- rnorm(n)

  # Six measurement equations (2 per factor for identification)
  # Factor 1 indicators
  T1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 0.5)  # Loading fixed to 1.0
  T2 <- 1.5 + 1.2 * f1 + rnorm(n, 0, 0.6)  # Free loading

  # Factor 2 indicators
  T3 <- 1.0 + 1.0 * f2 + rnorm(n, 0, 0.4)  # Loading fixed to 1.0
  T4 <- 0.8 + 0.9 * f2 + rnorm(n, 0, 0.5)  # Free loading

  # Factor 3 indicators
  T5 <- 1.2 + 1.0 * f3 + rnorm(n, 0, 0.5)  # Loading fixed to 1.0
  T6 <- 0.9 + 1.1 * f3 + rnorm(n, 0, 0.4)  # Free loading

  # Outcome with full second order terms
  # Y = intercept + beta*x + lambda1*f1 + lambda2*f2 + lambda3*f3
  #   + lambda_quad_1*f1^2 + lambda_quad_2*f2^2 + lambda_quad_3*f3^2
  #   + lambda_inter_1_2*f1*f2 + lambda_inter_1_3*f1*f3 + lambda_inter_2_3*f2*f3 + error
  true_intercept <- 3.0
  true_beta <- 0.5
  true_lambda <- c(0.8, 0.6, 0.7)  # Linear loadings
  true_lambda_quad <- c(0.25, 0.20, 0.15)  # Quadratic loadings
  true_lambda_inter <- c(0.3, 0.25, 0.2)  # Interaction loadings (1-2, 1-3, 2-3)
  true_sigma <- 0.5

  x <- rnorm(n)
  Y <- true_intercept + true_beta * x +
       true_lambda[1] * f1 + true_lambda[2] * f2 + true_lambda[3] * f3 +
       true_lambda_quad[1] * f1^2 + true_lambda_quad[2] * f2^2 + true_lambda_quad[3] * f3^2 +
       true_lambda_inter[1] * f1 * f2 + true_lambda_inter[2] * f1 * f3 + true_lambda_inter[3] * f2 * f3 +
       rnorm(n, 0, true_sigma)

  dat <- data.frame(
    intercept = 1, x = x,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4, T5 = T5, T6 = T6,
    Y = Y, eval = 1
  )

  # Create 3-factor model
  fm <- define_factor_model(n_factors = 3, n_types = 1)

  # Define measurement components
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0, 0.0),  # Factor 1 only, fixed to 1
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0, 0.0),  # Factor 1 only, free
    evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0, 0.0),  # Factor 2 only, fixed to 1
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_, 0.0),  # Factor 2 only, free
    evaluation_indicator = "eval"
  )
  mc_T5 <- define_model_component(
    name = "T5", data = dat, outcome = "T5", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 0.0, 1.0),  # Factor 3 only, fixed to 1
    evaluation_indicator = "eval"
  )
  mc_T6 <- define_model_component(
    name = "T6", data = dat, outcome = "T6", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 0.0, NA_real_),  # Factor 3 only, free
    evaluation_indicator = "eval"
  )

  # Outcome with full second order terms
  mc_Y <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_, NA_real_),  # All linear loadings free
    factor_spec = "full",  # Quadratic + interactions
    evaluation_indicator = "eval"
  )

  # Check component structure
  expect_equal(mc_Y$n_quadratic_loadings, 3)  # 3 factors -> 3 quadratic loadings
  expect_equal(mc_Y$n_interaction_loadings, 3)  # 3 factors -> 3 interactions (1-2, 1-3, 2-3)

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_T5, mc_T6, mc_Y),
    factor = fm
  )
  init <- initialize_parameters(ms, dat)

  # Check that all parameter names exist
  expect_true("Y_loading_quad_1" %in% init$param_names)
  expect_true("Y_loading_quad_2" %in% init$param_names)
  expect_true("Y_loading_quad_3" %in% init$param_names)
  expect_true("Y_loading_inter_1_2" %in% init$param_names)
  expect_true("Y_loading_inter_1_3" %in% init$param_names)
  expect_true("Y_loading_inter_2_3" %in% init$param_names)

  # Run gradient check
  grad_check <- check_gradient_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = GRAD_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(grad_check$pass,
              info = sprintf("Gradient check failed for 3-factor full model (max error: %.2e)",
                            grad_check$max_error))
})

# ==============================================================================
# Test 22: Hessian accuracy for 3-factor model with full second order terms
# ==============================================================================

test_that("Hessian accuracy for 3-factor linear model with full second order terms", {
  # This test runs by default (moderate computational cost)
  set.seed(322)

  # Generate data with 3 factors
  n <- 400
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  f3 <- rnorm(n)

  # Six measurement equations (2 per factor)
  T1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 0.5)
  T2 <- 1.5 + 1.2 * f1 + rnorm(n, 0, 0.6)
  T3 <- 1.0 + 1.0 * f2 + rnorm(n, 0, 0.4)
  T4 <- 0.8 + 0.9 * f2 + rnorm(n, 0, 0.5)
  T5 <- 1.2 + 1.0 * f3 + rnorm(n, 0, 0.5)
  T6 <- 0.9 + 1.1 * f3 + rnorm(n, 0, 0.4)

  # Outcome with full second order terms
  x <- rnorm(n)
  Y <- 3.0 + 0.5 * x +
       0.8 * f1 + 0.6 * f2 + 0.7 * f3 +
       0.25 * f1^2 + 0.20 * f2^2 + 0.15 * f3^2 +
       0.3 * f1 * f2 + 0.25 * f1 * f3 + 0.2 * f2 * f3 +
       rnorm(n, 0, 0.5)

  dat <- data.frame(
    intercept = 1, x = x,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4, T5 = T5, T6 = T6,
    Y = Y, eval = 1
  )

  fm <- define_factor_model(n_factors = 3, n_types = 1)

  # Measurement components
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T5 <- define_model_component(
    name = "T5", data = dat, outcome = "T5", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 0.0, 1.0),
    evaluation_indicator = "eval"
  )
  mc_T6 <- define_model_component(
    name = "T6", data = dat, outcome = "T6", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 0.0, NA_real_),
    evaluation_indicator = "eval"
  )

  mc_Y <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_T5, mc_T6, mc_Y),
    factor = fm
  )
  init <- initialize_parameters(ms, dat)

  # Run Hessian check
  hess_check <- check_hessian_accuracy(
    ms, dat, init$init_params,
    param_fixed = init$param_fixed,
    tol = HESS_TOL, verbose = VERBOSE, n_quad = 8
  )

  expect_true(hess_check$pass,
              info = sprintf("Hessian check failed for 3-factor full model (max error: %.2e)",
                            hess_check$max_error))
})

# ==============================================================================
# Test 23: Parameter recovery for 3-factor model with full second order terms
# ==============================================================================

test_that("3-factor linear estimation with full second order terms recovers parameters", {
  # Skip unless NOT_CRAN=true (extra tests only - computationally expensive with 8^3=512 quadrature points)
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"),
              "Skipping 3-factor estimation test (run with NOT_CRAN=true)")

  set.seed(323)

  # Generate data with 3 factors (moderate sample for parameter recovery)
  # Note: Uses n=800 to balance test speed with 8^3=512 quadrature points
  n <- 800
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  f3 <- rnorm(n)

  # True parameter values
  true_lambda_T2 <- 1.2
  true_lambda_T4 <- 0.9
  true_lambda_T6 <- 1.1

  true_Y_intercept <- 3.0
  true_Y_beta <- 0.5
  true_Y_lambda <- c(0.8, 0.6, 0.7)  # Linear loadings
  true_Y_quad <- c(0.25, 0.20, 0.15)  # Quadratic loadings
  true_Y_inter <- c(0.30, 0.25, 0.20)  # Interaction loadings (1-2, 1-3, 2-3)
  true_Y_sigma <- 0.5

  # Six measurement equations (2 per factor)
  T1 <- 2.0 + 1.0 * f1 + rnorm(n, 0, 0.5)
  T2 <- 1.5 + true_lambda_T2 * f1 + rnorm(n, 0, 0.6)
  T3 <- 1.0 + 1.0 * f2 + rnorm(n, 0, 0.4)
  T4 <- 0.8 + true_lambda_T4 * f2 + rnorm(n, 0, 0.5)
  T5 <- 1.2 + 1.0 * f3 + rnorm(n, 0, 0.5)
  T6 <- 0.9 + true_lambda_T6 * f3 + rnorm(n, 0, 0.4)

  # Outcome with full second order terms
  x <- rnorm(n)
  Y <- true_Y_intercept + true_Y_beta * x +
       true_Y_lambda[1] * f1 + true_Y_lambda[2] * f2 + true_Y_lambda[3] * f3 +
       true_Y_quad[1] * f1^2 + true_Y_quad[2] * f2^2 + true_Y_quad[3] * f3^2 +
       true_Y_inter[1] * f1 * f2 + true_Y_inter[2] * f1 * f3 + true_Y_inter[3] * f2 * f3 +
       rnorm(n, 0, true_Y_sigma)

  dat <- data.frame(
    intercept = 1, x = x,
    T1 = T1, T2 = T2, T3 = T3, T4 = T4, T5 = T5, T6 = T6,
    Y = Y, eval = 1
  )

  fm <- define_factor_model(n_factors = 3, n_types = 1)

  # Measurement components
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1.0, 0.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 1.0, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, NA_real_, 0.0),
    evaluation_indicator = "eval"
  )
  mc_T5 <- define_model_component(
    name = "T5", data = dat, outcome = "T5", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 0.0, 1.0),
    evaluation_indicator = "eval"
  )
  mc_T6 <- define_model_component(
    name = "T6", data = dat, outcome = "T6", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0.0, 0.0, NA_real_),
    evaluation_indicator = "eval"
  )

  mc_Y <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_T5, mc_T6, mc_Y),
    factor = fm
  )

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
  expect_equal(result$convergence, 0,
              info = "3-factor full model failed to converge")

  # Check that all second-order parameters exist in estimates
  expect_true("Y_loading_quad_1" %in% names(result$estimates))
  expect_true("Y_loading_quad_2" %in% names(result$estimates))
  expect_true("Y_loading_quad_3" %in% names(result$estimates))
  expect_true("Y_loading_inter_1_2" %in% names(result$estimates))
  expect_true("Y_loading_inter_1_3" %in% names(result$estimates))
  expect_true("Y_loading_inter_2_3" %in% names(result$estimates))

  # Check parameter recovery (tolerance of 0.4 for reasonable recovery with n=800)
  tol <- 0.4

  # Linear loadings
  est_lambda1 <- result$estimates["Y_loading_1"]
  est_lambda2 <- result$estimates["Y_loading_2"]
  est_lambda3 <- result$estimates["Y_loading_3"]
  expect_true(abs(est_lambda1 - true_Y_lambda[1]) < tol,
              info = sprintf("Linear loading 1 estimate (%.3f) far from true (%.3f)",
                            est_lambda1, true_Y_lambda[1]))
  expect_true(abs(est_lambda2 - true_Y_lambda[2]) < tol,
              info = sprintf("Linear loading 2 estimate (%.3f) far from true (%.3f)",
                            est_lambda2, true_Y_lambda[2]))
  expect_true(abs(est_lambda3 - true_Y_lambda[3]) < tol,
              info = sprintf("Linear loading 3 estimate (%.3f) far from true (%.3f)",
                            est_lambda3, true_Y_lambda[3]))

  # Quadratic loadings
  est_quad1 <- result$estimates["Y_loading_quad_1"]
  est_quad2 <- result$estimates["Y_loading_quad_2"]
  est_quad3 <- result$estimates["Y_loading_quad_3"]
  expect_true(abs(est_quad1 - true_Y_quad[1]) < tol,
              info = sprintf("Quadratic loading 1 estimate (%.3f) far from true (%.3f)",
                            est_quad1, true_Y_quad[1]))
  expect_true(abs(est_quad2 - true_Y_quad[2]) < tol,
              info = sprintf("Quadratic loading 2 estimate (%.3f) far from true (%.3f)",
                            est_quad2, true_Y_quad[2]))
  expect_true(abs(est_quad3 - true_Y_quad[3]) < tol,
              info = sprintf("Quadratic loading 3 estimate (%.3f) far from true (%.3f)",
                            est_quad3, true_Y_quad[3]))

  # Interaction loadings
  est_inter12 <- result$estimates["Y_loading_inter_1_2"]
  est_inter13 <- result$estimates["Y_loading_inter_1_3"]
  est_inter23 <- result$estimates["Y_loading_inter_2_3"]
  expect_true(abs(est_inter12 - true_Y_inter[1]) < tol,
              info = sprintf("Interaction loading 1-2 estimate (%.3f) far from true (%.3f)",
                            est_inter12, true_Y_inter[1]))
  expect_true(abs(est_inter13 - true_Y_inter[2]) < tol,
              info = sprintf("Interaction loading 1-3 estimate (%.3f) far from true (%.3f)",
                            est_inter13, true_Y_inter[2]))
  expect_true(abs(est_inter23 - true_Y_inter[3]) < tol,
              info = sprintf("Interaction loading 2-3 estimate (%.3f) far from true (%.3f)",
                            est_inter23, true_Y_inter[3]))

  if (VERBOSE) {
    cat("\n=== 3-Factor Full Model Parameter Recovery ===\n")
    cat("\nLinear loadings:\n")
    cat(sprintf("  Lambda 1: true=%.3f, est=%.3f\n", true_Y_lambda[1], est_lambda1))
    cat(sprintf("  Lambda 2: true=%.3f, est=%.3f\n", true_Y_lambda[2], est_lambda2))
    cat(sprintf("  Lambda 3: true=%.3f, est=%.3f\n", true_Y_lambda[3], est_lambda3))
    cat("\nQuadratic loadings:\n")
    cat(sprintf("  Quad 1: true=%.3f, est=%.3f\n", true_Y_quad[1], est_quad1))
    cat(sprintf("  Quad 2: true=%.3f, est=%.3f\n", true_Y_quad[2], est_quad2))
    cat(sprintf("  Quad 3: true=%.3f, est=%.3f\n", true_Y_quad[3], est_quad3))
    cat("\nInteraction loadings:\n")
    cat(sprintf("  Inter 1-2: true=%.3f, est=%.3f\n", true_Y_inter[1], est_inter12))
    cat(sprintf("  Inter 1-3: true=%.3f, est=%.3f\n", true_Y_inter[2], est_inter13))
    cat(sprintf("  Inter 2-3: true=%.3f, est=%.3f\n", true_Y_inter[3], est_inter23))
  }
})
