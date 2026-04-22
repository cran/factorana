# Tests for equality constraints (measurement invariance)

test_that("equality constraints are stored in model system", {
  set.seed(123)
  n <- 100
  dat <- data.frame(y1 = rnorm(n), y2 = rnorm(n), intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  # Define with equality constraints
  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  # Check that constraints are stored
  expect_equal(length(ms$equality_constraints), 1)
  expect_equal(ms$equality_constraints[[1]], c("m1_sigma", "m2_sigma"))
})

test_that("equality constraints require at least 2 parameters per group", {
  set.seed(124)
  n <- 100
  dat <- data.frame(y1 = rnorm(n), y2 = rnorm(n), intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  # Single parameter constraint should error
  expect_error(
    define_model_system(components = list(mc1, mc2), factor = fm,
                       equality_constraints = list(c("m1_sigma"))),
    regexp = "at least 2"
  )
})

test_that("tied parameters have identical estimates", {
  skip_on_cran()
  set.seed(125)
  n <- 500

  # Generate data with same true sigma for both components
  f <- rnorm(n, 0, 1)
  true_sigma <- 0.6
  y1 <- 1.0 * f + rnorm(n, 0, true_sigma)
  y2 <- 0.8 * f + rnorm(n, 0, true_sigma)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  # Constrain sigmas to be equal
  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)
  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                 optimizer = "nlminb", verbose = FALSE)

  expect_equal(result$convergence, 0,
               info = sprintf("Model with equality constraints did not converge, code: %d", result$convergence))

  # Check that tied parameters are exactly equal
  sigma1 <- unname(result$estimates["m1_sigma"])
  sigma2 <- unname(result$estimates["m2_sigma"])

  expect_equal(sigma1, sigma2, tolerance = 1e-10,
               info = sprintf("Tied sigmas should be exactly equal: %.6f vs %.6f", sigma1, sigma2))

  # Check that estimated sigma is close to true value
  expect_true(abs(sigma1 - true_sigma) < 0.15,
              info = sprintf("Estimated sigma %.3f differs from true %.3f by more than 0.15",
                            sigma1, true_sigma))
})

test_that("gradient is correct with equality constraints", {
  skip_on_cran()
  set.seed(126)
  n <- 300

  f <- rnorm(n, 0, 1)
  y1 <- 1.0 * f + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  # First param (factor_var) is fixed when loading is fixed
  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 1e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass,
              info = sprintf("Gradient check with equality constraints failed, max error: %.2e", result$max_error))
})

test_that("Hessian is correct with equality constraints", {
  skip_on_cran()
  set.seed(127)
  n <- 300

  f <- rnorm(n, 0, 1)
  y1 <- 1.0 * f + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 1.0, verbose = FALSE, n_quad = 8)
  expect_true(result$pass,
              info = sprintf("Hessian check with equality constraints failed, max error: %.2e", result$max_error))
})
