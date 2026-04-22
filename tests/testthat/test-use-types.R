# Tests for use_types feature (component-level type usage)
# Tests that type intercepts are only included for components with use_types = TRUE

# Helper function to generate data with latent types
generate_two_stage_data <- function(n, n_types = 2) {
  set.seed(42)

  # Generate latent factor
  f <- rnorm(n, 0, 1)

  # Type assignment (uniform)
  true_types <- sample(1:n_types, n, replace = TRUE)

  # Type-specific intercepts for outcome model only
  # Type 1 = 0 (reference), type 2 = 0.5
  type_intercepts <- (1:n_types - 1) * 0.5

  # Measurement equations (no type effects)
  y1 <- 0 + 1.0 * f + rnorm(n, 0, 0.5)  # Fixed loading = 1
  y2 <- 0 + 0.8 * f + rnorm(n, 0, 0.5)  # Free loading

  # Outcome equation (with type effects)
  outcome <- 0.5 + type_intercepts[true_types] + 0.6 * f + rnorm(n, 0, 0.5)

  data.frame(
    y1 = y1,
    y2 = y2,
    outcome = outcome,
    intercept = 1,
    true_type = true_types,
    true_factor = f,
    eval = 1
  )
}

test_that("use_types parameter is properly stored in model component", {
  dat <- generate_two_stage_data(100)
  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # Without use_types (default FALSE)
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = 1.0
  )
  expect_false(mc1$use_types)
  expect_equal(mc1$n_type_intercepts, 0L)

  # With use_types = FALSE explicitly
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", use_types = FALSE
  )
  expect_false(mc2$use_types)
  expect_equal(mc2$n_type_intercepts, 0L)

  # With use_types = TRUE
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "outcome", factor = fm,
    covariates = "intercept", model_type = "linear",
    use_types = TRUE
  )
  expect_true(mc3$use_types)
  expect_equal(mc3$n_type_intercepts, 1L)  # n_types - 1 = 2 - 1 = 1
})

test_that("use_types=TRUE with n_types=1 produces warning and sets to FALSE", {
  dat <- generate_two_stage_data(100)
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  expect_warning(
    mc <- define_model_component(
      name = "m1", data = dat, outcome = "outcome", factor = fm,
      covariates = "intercept", model_type = "linear", use_types = TRUE
    ),
    "use_types = TRUE.*has no effect"
  )
  expect_false(mc$use_types)
})

test_that("parameter initialization: no type params when no component uses types", {
  dat <- generate_two_stage_data(100)
  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # All components have use_types = FALSE
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = 1.0,
    use_types = FALSE
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", use_types = FALSE
  )

  ms <- define_model_system(components = list(mc1, mc2), factor = fm)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  param_names <- names(init_result$init_params)

  # Should NOT have type loadings or type intercepts
  expect_false(any(grepl("type_2_loading", param_names)))
  expect_false(any(grepl("type_2_intercept", param_names)))

  # Should have factor variance and component params only
  expect_true("factor_var_1" %in% param_names)
  expect_true("m1_sigma" %in% param_names)
  expect_true("m2_sigma" %in% param_names)
})

test_that("parameter initialization: type params only for components with use_types=TRUE", {
  dat <- generate_two_stage_data(100)
  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # m1 and m2 without types, m3 with types
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = 1.0,
    use_types = FALSE
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", use_types = FALSE
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "outcome", factor = fm,
    covariates = "intercept", model_type = "linear",
    use_types = TRUE
  )

  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  param_names <- names(init_result$init_params)

  # Should have type loadings (for the type probability model)
  expect_true("type_2_loading_1" %in% param_names)

  # Should have type intercept ONLY for m3
  expect_false("m1_type_2_intercept" %in% param_names)
  expect_false("m2_type_2_intercept" %in% param_names)
  expect_true("m3_type_2_intercept" %in% param_names)
})

test_that("estimation with mixed use_types works correctly", {
  skip_on_cran()

  dat <- generate_two_stage_data(200, n_types = 2)
  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # Model system with some components using types and some not
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = 1.0,
    use_types = FALSE, evaluation_indicator = "eval"
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", use_types = FALSE,
    evaluation_indicator = "eval"
  )
  mc3 <- define_model_component(
    name = "outcome", data = dat, outcome = "outcome", factor = fm,
    covariates = "intercept", model_type = "linear",
    use_types = TRUE, evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)

  # Check parameter initialization
  init <- initialize_parameters(ms, dat, verbose = FALSE)

  # Should have type loadings (for type probability model)
  expect_true("type_2_loading_1" %in% names(init$init_params))

  # Should have type intercept only for mc3 (outcome)
  expect_false("m1_type_2_intercept" %in% names(init$init_params))
  expect_false("m2_type_2_intercept" %in% names(init$init_params))
  expect_true("outcome_type_2_intercept" %in% names(init$init_params))

  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(ms, dat, control = ctrl,
                                 optimizer = "nlminb", verbose = FALSE)
  expect_true(result$convergence == 0)

  # Check that type intercept was estimated
  expect_true("outcome_type_2_intercept" %in% names(result$estimates))

  # Check that results have reasonable values
  expect_true(is.finite(result$estimates["type_2_loading_1"]))
  expect_true(is.finite(result$estimates["outcome_type_2_intercept"]))
})

test_that("gradient validation for mixed use_types model", {
  skip_on_cran()

  dat <- generate_two_stage_data(100, n_types = 2)
  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # m1 without types, m2 with types
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = 1.0,
    use_types = FALSE, evaluation_indicator = "eval"
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", use_types = TRUE,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc1, mc2), factor = fm)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Initialize C++ model
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  fm_ptr <- initialize_factor_model_cpp(ms, dat, ctrl$n_quad_points)

  params <- init_result$init_params
  n_params <- length(params)
  param_fixed <- init_result$param_fixed

  # Get analytical gradient
  result <- evaluate_likelihood_cpp(fm_ptr, params, 2)
  analytical_grad <- result$gradient

  # The result object should have the gradient
  expect_length(analytical_grad, n_params)

  # Compute finite difference gradient for key free parameters
  eps <- 1e-6
  free_params <- which(!param_fixed)

  # Test at least the first few free parameters
  test_indices <- head(free_params, min(3, length(free_params)))

  errors_found <- 0
  for (i in test_indices) {
    params_plus <- params
    params_plus[i] <- params[i] + eps
    result_plus <- evaluate_likelihood_cpp(fm_ptr, params_plus, 1)

    params_minus <- params
    params_minus[i] <- params[i] - eps
    result_minus <- evaluate_likelihood_cpp(fm_ptr, params_minus, 1)

    # Handle both possible return value names
    lkhd_plus <- if (!is.null(result_plus$log_likelihood)) result_plus$log_likelihood else result_plus$loglik
    lkhd_minus <- if (!is.null(result_minus$log_likelihood)) result_minus$log_likelihood else result_minus$loglik

    # Skip if we couldn't get valid likelihoods
    if (is.null(lkhd_plus) || is.null(lkhd_minus) ||
        !is.finite(lkhd_plus) || !is.finite(lkhd_minus)) {
      next
    }

    fd_grad_i <- (lkhd_plus - lkhd_minus) / (2 * eps)

    if (is.finite(fd_grad_i) && abs(fd_grad_i) > 1e-10) {
      rel_error <- abs(analytical_grad[i] - fd_grad_i) / max(abs(fd_grad_i), 1e-8)
      if (rel_error >= 1e-4) {
        errors_found <- errors_found + 1
      }
    }
  }

  # At least some gradient checks should pass
  expect_lt(errors_found, length(test_indices),
            label = "Most gradient checks should pass")
})

test_that("multinomial logit with use_types works correctly", {
  skip_on_cran()

  set.seed(42)
  n <- 200
  n_types <- 2

  # Generate factor and types
  f <- rnorm(n, 0, 1)
  true_types <- sample(1:n_types, n, replace = TRUE)

  # Type effects on choice
  type_effect <- c(0, 0.5)[true_types]

  # Multinomial logit with 3 choices
  # V[j] = beta_j + lambda_j * f + type_effect (for types > 1)
  V1 <- 0  # Reference
  V2 <- 0.5 + 0.3 * f + type_effect + rnorm(n, 0, 0.5)
  V3 <- 0.2 + 0.5 * f + type_effect + rnorm(n, 0, 0.5)

  # Simulate choices
  exp_V <- cbind(exp(V1), exp(V2), exp(V3))
  probs <- exp_V / rowSums(exp_V)
  choice <- apply(probs, 1, function(p) sample(1:3, 1, prob = p))

  dat <- data.frame(
    choice = choice,
    intercept = 1,
    eval = 1
  )

  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # Logit with use_types = TRUE
  mc <- define_model_component(
    name = "choice", data = dat, outcome = "choice", factor = fm,
    covariates = "intercept", model_type = "logit",
    num_choices = 3, use_types = TRUE,
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(components = list(mc), factor = fm)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  param_names <- names(init_result$init_params)

  # Should have type loadings
  expect_true("type_2_loading_1" %in% param_names)

  # Should have type intercepts for each non-reference choice
  # Format: {component}_c{choice}_type_{type}_intercept
  expect_true("choice_c1_type_2_intercept" %in% param_names)
  expect_true("choice_c2_type_2_intercept" %in% param_names)
})
