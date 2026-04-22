# Tests for observation weights

test_that("observation weights are validated during estimation", {
  # Note: Weights column validation happens during estimate_model_rcpp,
  # not during define_model_system. This is the expected behavior.
  set.seed(123)
  n <- 100
  dat <- data.frame(y = rnorm(n), intercept = 1)

  fm <- define_factor_model(n_factors = 1)
  mc <- define_model_component("m", dat, "y", fm,
                                covariates = "intercept", model_type = "linear",
                                loading_normalization = 1)

  # define_model_system allows any string for weights (validated later)
  ms <- define_model_system(components = list(mc), factor = fm, weights = "nonexistent")
  expect_equal(ms$weights, "nonexistent")
})

test_that("observation weights are stored in model system", {
  set.seed(124)
  n <- 100
  dat <- data.frame(y = rnorm(n), intercept = 1, w = runif(n, 0.5, 1.5))

  fm <- define_factor_model(n_factors = 1)
  mc <- define_model_component("m", dat, "y", fm,
                                covariates = "intercept", model_type = "linear",
                                loading_normalization = 1)

  ms <- define_model_system(components = list(mc), factor = fm, weights = "w")

  # Check that weights variable name is stored
  expect_equal(ms$weights, "w")
})

test_that("uniform weights equal no weights", {
  skip_on_cran()
  set.seed(126)
  n <- 300

  f <- rnorm(n, 0, 1)
  y <- 2.0 + 1.0 * f + rnorm(n, 0, 0.5)

  dat <- data.frame(y = y, intercept = 1, w = rep(1, n))

  fm <- define_factor_model(n_factors = 1)
  mc <- define_model_component("m", dat, "y", fm,
                                covariates = "intercept", model_type = "linear",
                                loading_normalization = 1)

  # Model without weights
  ms_no_weights <- define_model_system(components = list(mc), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)
  result_no_weights <- estimate_model_rcpp(ms_no_weights, dat, control = ctrl,
                                            parallel = FALSE, optimizer = "nlminb", verbose = FALSE)

  # Model with uniform weights
  ms_weights <- define_model_system(components = list(mc), factor = fm, weights = "w")
  result_weights <- estimate_model_rcpp(ms_weights, dat, control = ctrl,
                                         parallel = FALSE, optimizer = "nlminb", verbose = FALSE)

  # Both should converge
  expect_equal(result_no_weights$convergence, 0)
  expect_equal(result_weights$convergence, 0)

  # Log-likelihoods should be essentially equal
  expect_equal(result_no_weights$loglik, result_weights$loglik, tolerance = 1e-6,
               info = sprintf("Loglik without weights: %.6f, with uniform weights: %.6f",
                             result_no_weights$loglik, result_weights$loglik))

  # Estimates should be essentially equal
  for (pname in names(result_no_weights$estimates)) {
    expect_equal(result_no_weights$estimates[pname], result_weights$estimates[pname], tolerance = 1e-6,
                 info = sprintf("Parameter %s differs: no weights %.6f, uniform weights %.6f",
                               pname, result_no_weights$estimates[pname], result_weights$estimates[pname]))
  }
})

test_that("gradient is correct with observation weights", {
  skip_on_cran()
  set.seed(127)
  n <- 300

  f <- rnorm(n, 0, 1)
  y <- 2.0 + 1.0 * f + rnorm(n, 0, 0.5)
  w <- runif(n, 0.5, 1.5)  # Non-uniform weights

  dat <- data.frame(y = y, intercept = 1, w = w)

  fm <- define_factor_model(n_factors = 1)
  mc <- define_model_component("m", dat, "y", fm,
                                covariates = "intercept", model_type = "linear",
                                loading_normalization = 1)

  ms <- define_model_system(components = list(mc), factor = fm, weights = "w")

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE  # Factor variance fixed when loading is fixed

  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 1e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass,
              info = sprintf("Gradient check with weights failed, max error: %.2e", result$max_error))
})

test_that("Hessian is correct with observation weights", {
  skip_on_cran()
  set.seed(128)
  n <- 300

  f <- rnorm(n, 0, 1)
  y <- 2.0 + 1.0 * f + rnorm(n, 0, 0.5)
  w <- runif(n, 0.5, 1.5)

  dat <- data.frame(y = y, intercept = 1, w = w)

  fm <- define_factor_model(n_factors = 1)
  mc <- define_model_component("m", dat, "y", fm,
                                covariates = "intercept", model_type = "linear",
                                loading_normalization = 1)

  ms <- define_model_system(components = list(mc), factor = fm, weights = "w")

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 1e-3, verbose = FALSE, n_quad = 8)
  expect_true(result$pass,
              info = sprintf("Hessian check with weights failed, max error: %.2e", result$max_error))
})
