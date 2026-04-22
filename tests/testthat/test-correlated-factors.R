# Tests for correlated factor models

# Helper function to isolate tests by clearing C++ external pointers
isolate_test <- function() {
  rm(list = ls(pattern = "^(fm|mc|ms|result|ctrl|init)", envir = parent.frame()),
     envir = parent.frame())
  gc(verbose = FALSE, full = TRUE)
  invisible(NULL)
}

test_that("correlated two-factor model initializes correctly", {
  isolate_test()
  set.seed(123)
  n <- 100

  # Generate simple data
  dat <- data.frame(
    y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), y4 = rnorm(n),
    intercept = 1
  )

  # Define 2-factor model WITH correlation

  fm <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "correlation")

  # Check factor model has correlation flag set
  expect_true(fm$correlation)
  expect_equal(fm$n_factors, 2)

  # Define components
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0)
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0)
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1)
  )
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA)
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4),
    factor = fm
  )

  # Initialize parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Check that correlation parameter was added
  expect_true("factor_corr_1_2" %in% names(init_result$init_params))

  # Check initial correlation value is reasonable (should be 0 or small)
  corr_init <- init_result$init_params["factor_corr_1_2"]
  expect_true(abs(corr_init) < 0.5)
})

test_that("correlated two-factor linear model converges and recovers correlation", {
  isolate_test()
  skip_on_cran()
  set.seed(42)
  n <- 300

  # True correlation between factors
  rho_true <- 0.5

  # Generate correlated factors using Cholesky decomposition
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  f1 <- z1
  f2 <- rho_true * z1 + sqrt(1 - rho_true^2) * z2

  # Generate linear measures
  sigma <- 0.5
  y1 <- 1.0 * f1 + rnorm(n, 0, sigma)
  y2 <- 0.8 * f1 + rnorm(n, 0, sigma)
  y3 <- 1.0 * f2 + rnorm(n, 0, sigma)
  y4 <- 0.9 * f2 + rnorm(n, 0, sigma)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  # Define 2-factor model WITH correlation
  fm <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "correlation")

  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0)
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0)
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1)
  )
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA)
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4),
    factor = fm
  )

  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check convergence (0 = success)
  # After fixing Hessian computation for correlated factors, convergence code 0 is expected
  expect_equal(result$convergence, 0,
    info = sprintf("Expected convergence code 0, got: %d", result$convergence))

  # Check correlation estimate is close to true value (within 0.2)
  est_corr <- result$estimates["factor_corr_1_2"]
  expect_true(abs(est_corr - rho_true) < 0.2,
    info = sprintf("Estimated correlation %.3f differs from true %.3f by more than 0.2",
                   est_corr, rho_true))

  # Check factor variances (true = 1.0)
  est_var1 <- result$estimates["factor_var_1"]
  est_var2 <- result$estimates["factor_var_2"]
  expect_true(abs(est_var1 - 1.0) < 0.20,
    info = sprintf("Factor 1 variance %.3f should be close to 1.0", est_var1))
  expect_true(abs(est_var2 - 1.0) < 0.20,
    info = sprintf("Factor 2 variance %.3f should be close to 1.0", est_var2))

  # Check that loading estimates are reasonable
  est_loading_m2 <- result$estimates["m2_loading_1"]
  expect_true(abs(est_loading_m2 - 0.8) < 0.15,
    info = sprintf("Loading m2 estimate %.3f differs from true 0.8", est_loading_m2))

  est_loading_m4 <- result$estimates["m4_loading_2"]
  expect_true(abs(est_loading_m4 - 0.9) < 0.15,
    info = sprintf("Loading m4 estimate %.3f differs from true 0.9", est_loading_m4))

  # Check sigma estimates (true = 0.5)
  est_sigma_m1 <- result$estimates["m1_sigma"]
  est_sigma_m2 <- result$estimates["m2_sigma"]
  est_sigma_m3 <- result$estimates["m3_sigma"]
  est_sigma_m4 <- result$estimates["m4_sigma"]
  expect_true(abs(est_sigma_m1 - sigma) < 0.10,
    info = sprintf("m1 sigma %.3f should be close to %.1f", est_sigma_m1, sigma))
  expect_true(abs(est_sigma_m2 - sigma) < 0.10,
    info = sprintf("m2 sigma %.3f should be close to %.1f", est_sigma_m2, sigma))
  expect_true(abs(est_sigma_m3 - sigma) < 0.10,
    info = sprintf("m3 sigma %.3f should be close to %.1f", est_sigma_m3, sigma))
  expect_true(abs(est_sigma_m4 - sigma) < 0.10,
    info = sprintf("m4 sigma %.3f should be close to %.1f", est_sigma_m4, sigma))
})

test_that("correlated two-factor model recovers factor variances with sufficient quadrature", {
  isolate_test()
  skip_on_cran()
  set.seed(123)
  n <- 500

  # True parameters
  rho_true <- 0.6
  true_var1 <- 1.0
  true_var2 <- 1.0

  # Generate correlated factors using Cholesky decomposition
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  f1 <- z1  # Var = 1
  f2 <- rho_true * z1 + sqrt(1 - rho_true^2) * z2  # Var = 1, Cor = rho

  # 3 indicators per factor for reliable variance estimation
  sigma <- 0.5
  y1 <- 1.0 * f1 + rnorm(n, 0, sigma)
  y2 <- 0.8 * f1 + rnorm(n, 0, sigma)
  y3 <- 1.2 * f1 + rnorm(n, 0, sigma)
  y4 <- 1.0 * f2 + rnorm(n, 0, sigma)
  y5 <- 0.9 * f2 + rnorm(n, 0, sigma)
  y6 <- 1.1 * f2 + rnorm(n, 0, sigma)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, y5 = y5, y6 = y6, intercept = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "correlation")

  mc1 <- define_model_component(name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = c(1, 0))
  mc2 <- define_model_component(name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = c(NA, 0))
  mc3 <- define_model_component(name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = c(NA, 0))
  mc4 <- define_model_component(name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = c(0, 1))
  mc5 <- define_model_component(name = "m5", data = dat, outcome = "y5", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = c(0, NA))
  mc6 <- define_model_component(name = "m6", data = dat, outcome = "y6", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4, mc5, mc6), factor = fm)

  # Use 16 quadrature points for accurate variance estimation with correlated factors
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                optimizer = "nlminb", verbose = FALSE)

  expect_equal(result$convergence, 0)

  # Check factor variances are close to 1.0 (tolerance 0.25 for n=500)
  est_var1 <- result$estimates["factor_var_1"]
  est_var2 <- result$estimates["factor_var_2"]
  est_corr <- result$estimates["factor_corr_1_2"]

  expect_true(abs(est_var1 - true_var1) < 0.25,
    info = sprintf("Factor 1 variance %.3f should be close to %.1f", est_var1, true_var1))
  expect_true(abs(est_var2 - true_var2) < 0.25,
    info = sprintf("Factor 2 variance %.3f should be close to %.1f", est_var2, true_var2))
  expect_true(abs(est_corr - rho_true) < 0.15,
    info = sprintf("Correlation %.3f should be close to %.1f", est_corr, rho_true))
})

# NOTE: Three-factor correlated model test removed.
# Correlated factor models are only supported for n_factors = 2.
# The define_factor_model() function now throws an error if
# correlation = TRUE and n_factors > 2.

test_that("uncorrelated model (factor_structure = 'independent') has no correlation parameters", {
  isolate_test()
  set.seed(123)
  n <- 100

  dat <- data.frame(
    y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), y4 = rnorm(n),
    intercept = 1
  )

  # Define 2-factor model WITHOUT correlation (default)
  fm <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "independent")

  # Check factor model has correlation flag NOT set
  expect_false(fm$correlation)

  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0)
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0)
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1)
  )
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA)
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4),
    factor = fm
  )

  # Initialize parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Check that NO correlation parameters were added
  corr_params <- grep("factor_corr", names(init_result$init_params), value = TRUE)
  expect_equal(length(corr_params), 0)
})

test_that("negative correlation is recovered correctly", {
  isolate_test()
  skip_on_cran()
  set.seed(456)
  n <- 300

  # True NEGATIVE correlation between factors
  rho_true <- -0.4

  # Generate correlated factors
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  f1 <- z1
  f2 <- rho_true * z1 + sqrt(1 - rho_true^2) * z2

  # Generate linear measures
  sigma <- 0.5
  y1 <- 1.0 * f1 + rnorm(n, 0, sigma)
  y2 <- 0.8 * f1 + rnorm(n, 0, sigma)
  y3 <- 1.0 * f2 + rnorm(n, 0, sigma)
  y4 <- 0.9 * f2 + rnorm(n, 0, sigma)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "correlation")

  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0)
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0)
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1)
  )
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA)
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4),
    factor = fm
  )

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check convergence (0 = success)
  # After fixing Hessian computation for correlated factors, convergence code 0 is expected
  expect_equal(result$convergence, 0,
    info = sprintf("Expected convergence code 0, got: %d", result$convergence))

  # Check that estimated correlation is negative and close to true value
  est_corr <- result$estimates["factor_corr_1_2"]
  expect_true(est_corr < 0,
    info = sprintf("Expected negative correlation, got %.3f", est_corr))
  expect_true(abs(est_corr - rho_true) < 0.2,
    info = sprintf("Estimated correlation %.3f differs from true %.3f by more than 0.2",
                   est_corr, rho_true))
})
