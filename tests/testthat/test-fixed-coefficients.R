test_that("fix_coefficient validates inputs correctly", {
  skip_on_cran()

  set.seed(123)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 1)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  # Valid usage
  mc_fixed <- fix_coefficient(mc, "intercept", 0)
  expect_s3_class(mc_fixed, "model_component")
  expect_equal(length(mc_fixed$fixed_coefficients), 1)
  expect_equal(mc_fixed$fixed_coefficients[[1]]$covariate, "intercept")
  expect_equal(mc_fixed$fixed_coefficients[[1]]$value, 0)

  # Error on non-existent covariate
  expect_error(
    fix_coefficient(mc, "nonexistent", 0),
    "not found in component"
  )

  # Error on non-model_component
  expect_error(
    fix_coefficient(list(a = 1), "x", 0),
    "must be an object of class"
  )

  # Error on non-numeric value
  expect_error(
    fix_coefficient(mc, "intercept", "zero"),
    "must be a single finite numeric"
  )
})

test_that("Fixed coefficients work in linear model", {
  skip_on_cran()

  set.seed(456)
  n <- 300

  # True DGP: Y = 2.0 + 0*x1 + 1.5*x2 + f + eps
  # We will fix x1 coefficient to 0
  f <- rnorm(n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  Y <- 2.0 + 0*x1 + 1.5*x2 + f + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept = 1, x1 = x1, x2 = x2, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 1)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1", "x2"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  # Count parameters before fixing
  # Linear: fac_var(1) + intercept(1) + x1(1) + x2(1) + loading(1) + sigma(1) = 6
  ms_before <- define_model_system(components = list(mc), factor = fm)
  n_params_before <- sum(sapply(ms_before$components, function(c) c$nparam_model)) + fm$n_factors

  # Fix x1 coefficient to 0
  mc_fixed <- fix_coefficient(mc, "x1", 0)

  ms <- define_model_system(components = list(mc_fixed), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result <- estimate_model_rcpp(
    ms, dat,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0)
  expect_true(!is.null(result))
  expect_true(!is.null(result$estimates))

  # Verify one fewer parameter is estimated (x1 is fixed)
  # nparam_model should be reduced by 1 in the fixed component
  expect_equal(mc_fixed$nparam_model, mc$nparam_model - 1)

  # Result estimates includes ALL parameters (including fixed ones)
  # Fixed params appear at their fixed values
  expect_equal(length(result$estimates), n_params_before)

  # Check the fixed param (Y_x1) is at its fixed value
  expect_equal(result$estimates["Y_x1"], c(Y_x1 = 0), tolerance = 1e-10)
})

test_that("Fixed coefficients work in probit model", {
  skip_on_cran()

  set.seed(789)
  n <- 300

  # True DGP: two binary indicators of a common factor.
  #   Y  = I(0.5 + 0*x1 + f + eps > 0)   -- loading fixed to 1 for identification
  #   Y2 = I(-0.2 + 0.8*f + eps > 0)      -- second indicator so the factor
  #                                           variance is actually identified
  # A single probit indicator with fixed loading = 1 and free factor variance
  # is marginally identified only in the ratio intercept / sqrt(1 + factor_var),
  # which produces a flat likelihood and nlminb stops on relative (not strict)
  # convergence. Adding Y2 fixes this.
  f <- rnorm(n)
  x1 <- rnorm(n)
  Y <- as.integer(0.5 + 0*x1 + f + rnorm(n) > 0)
  Y2 <- as.integer(-0.2 + 0.8*f + rnorm(n) > 0)

  dat <- data.frame(intercept = 1, x1 = x1, Y = Y, Y2 = Y2, eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 1)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "probit",
    loading_normalization = 1, evaluation_indicator = "eval"
  )
  mc_helper <- define_model_component(
    name = "Y2", data = dat, outcome = "Y2", factor = fm,
    covariates = "intercept", model_type = "probit",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  # Fix x1 coefficient to 0 on the primary component
  mc_fixed <- fix_coefficient(mc, "x1", 0)

  ms <- define_model_system(components = list(mc_fixed, mc_helper), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(
    ms, dat,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0)
  expect_true(!is.null(result))

  # Check the fixed param is at its fixed value
  expect_equal(result$estimates["Y_x1"], c(Y_x1 = 0), tolerance = 1e-10)
})

test_that("Multiple fixed coefficients work", {
  skip_on_cran()

  set.seed(111)
  n <- 200

  f <- rnorm(n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  Y <- 1.0 + 0*x1 + 0*x2 + f + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept = 1, x1 = x1, x2 = x2, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 1)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1", "x2"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  # Count parameters before fixing
  ms_before <- define_model_system(components = list(mc), factor = fm)
  n_params_before <- sum(sapply(ms_before$components, function(c) c$nparam_model)) + fm$n_factors

  # Fix both x1 and x2 to 0
  mc_fixed <- fix_coefficient(mc, "x1", 0)
  mc_fixed <- fix_coefficient(mc_fixed, "x2", 0)

  expect_equal(length(mc_fixed$fixed_coefficients), 2)
  expect_equal(mc_fixed$nparam_model, mc$nparam_model - 2)

  ms <- define_model_system(components = list(mc_fixed), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result <- estimate_model_rcpp(
    ms, dat,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0)
  expect_true(!is.null(result))

  # Verify parameter count includes all params (fixed ones at their fixed values)
  expect_equal(length(result$estimates), n_params_before)

  # Check the fixed params are at their fixed values
  expect_equal(result$estimates["Y_x1"], c(Y_x1 = 0), tolerance = 1e-10)
  expect_equal(result$estimates["Y_x2"], c(Y_x2 = 0), tolerance = 1e-10)
})

test_that("Fixed coefficients with choice for multinomial logit", {
  skip_on_cran()

  set.seed(222)
  n <- 300

  # Multinomial logit with 3 choices
  f <- rnorm(n)
  x1 <- rnorm(n)

  # Choice probabilities
  V1 <- 0  # Reference
  V2 <- 0.5 + 0*x1 + 0.8*f  # x1 coef fixed to 0 for choice 2
  V3 <- 1.0 + 0.3*x1 + 1.2*f

  exp_V <- cbind(exp(V1), exp(V2), exp(V3))
  probs <- exp_V / rowSums(exp_V)

  # Draw choices
  Y <- sapply(1:n, function(i) sample(1:3, 1, prob = probs[i, ]))

  dat <- data.frame(intercept = 1, x1 = x1, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 1)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "logit",
    num_choices = 3,
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  # Count parameters before fixing
  # Mlogit 3 choices: fac_var(1) + 2 alternatives * (intercept + x1 + loading) = 1 + 2*3 = 7
  ms_before <- define_model_system(components = list(mc), factor = fm)
  n_params_before <- sum(sapply(ms_before$components, function(c) c$nparam_model)) + fm$n_factors

  # Fix x1 coefficient to 0 for choice 2 (first non-reference choice)
  mc_fixed <- fix_coefficient(mc, "x1", 0, choice = 1)

  expect_equal(mc_fixed$fixed_coefficients[[1]]$choice, 1)
  expect_equal(mc_fixed$nparam_model, mc$nparam_model - 1)

  ms <- define_model_system(components = list(mc_fixed), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result <- estimate_model_rcpp(
    ms, dat,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0)
  expect_true(!is.null(result))

  # Verify parameter count includes all params (fixed ones at their fixed values)
  expect_equal(length(result$estimates), n_params_before)

  # Check the fixed param is at its fixed value
  expect_equal(result$estimates["Y_c1_x1"], c(Y_c1_x1 = 0), tolerance = 1e-10)
})

test_that("Fixed coefficient value is used correctly", {
  skip_on_cran()

  set.seed(333)
  n <- 500

  # Generate data where the true intercept is 2.5
  f <- rnorm(n)
  x1 <- rnorm(n)
  Y <- 2.5 + 1.0*x1 + 0.8*f + rnorm(n, 0, 0.3)

  dat <- data.frame(intercept = 1, x1 = x1, Y = Y, eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 1)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval"
  )

  # Fix intercept to 2.5 (the true value)
  mc_fixed <- fix_coefficient(mc, "intercept", 2.5)

  ms <- define_model_system(components = list(mc_fixed), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result <- estimate_model_rcpp(
    ms, dat,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  expect_true(!is.null(result))

  # Should have converged successfully
  expect_equal(result$convergence, 0)

  # The fixed intercept should be exactly at its fixed value
  expect_equal(result$estimates["Y_intercept"], c(Y_intercept = 2.5), tolerance = 1e-10)

  # Estimated x1 coefficient should be close to 1.0
  expect_true(abs(result$estimates["Y_x1"] - 1.0) < 0.3)
})
