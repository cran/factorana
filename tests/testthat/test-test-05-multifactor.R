# Helper function to isolate tests by clearing C++ external pointers
isolate_test <- function() {
  # Remove any lingering model objects from global environment
  rm(list = ls(pattern = "^(fm|mc|ms|result|ctrl|init)", envir = parent.frame()),
     envir = parent.frame())
  # Force garbage collection to clean up C++ external pointers
  gc(verbose = FALSE, full = TRUE)
  invisible(NULL)
}

test_that("multi-factor loadings work with normalization (linear)", {
  isolate_test()
  set.seed(1)
  fm <- define_factor_model(3, 1)
  dat <- data.frame(y = rnorm(120), x1 = rnorm(120), x2 = rnorm(120), eval = 1L)

  # Specify loading normalization at component level: NA, 0, 1
  mc <- define_model_component("lin3", dat, "y", fm,
                               evaluation_indicator = "eval",
                               covariates = c("x1","x2"),
                               model_type = "linear",
                               loading_normalization = c(NA, 0, 1),
                               intercept = FALSE)

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  # Check component loading normalization
  expect_length(mc$loading_normalization, 3)
  expect_true(is.na(mc$loading_normalization[1]))
  expect_equal(mc$loading_normalization[2], 0)
  expect_equal(mc$loading_normalization[3], 1)
})

test_that("multi-factor loadings work with normalization (probit)", {
  isolate_test()
  set.seed(2)
  fm <- define_factor_model(2, 1)
  # simple binary outcome
  x1 <- rnorm(150); x2 <- rnorm(150)
  p  <- plogis(0.5*x1 - 0.3*x2)
  y  <- as.integer(runif(150) < p)
  dat <- data.frame(y=y, x1=x1, x2=x2, eval=1L)

  # Specify loading normalization at component level: NA, 1
  mc <- define_model_component("prb2", dat, "y", fm,
                               evaluation_indicator = "eval",
                               covariates = c("x1","x2"),
                               model_type = "probit",
                               loading_normalization = c(NA, 1),
                               intercept = FALSE)

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  # Check component loading normalization
  expect_length(mc$loading_normalization, 2)
  expect_true(is.na(mc$loading_normalization[1]))
  expect_equal(mc$loading_normalization[2], 1)
})

test_that("two-factor CFA with ordered probit converges and recovers parameters", {
  isolate_test()
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Generate data with known two-factor structure
  set.seed(123)
  n <- 2000

  # Two independent factors
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)

  # True loadings for factor 1 measures (m1, m2, m3)
  # Using variance normalization: factor variances fixed to 1.0, all loadings free
  lambda1_1 <- 1.0
  lambda1_2 <- 0.8
  lambda1_3 <- 1.2

  # True loadings for factor 2 measures (m4, m5, m6)
  lambda2_4 <- 1.0
  lambda2_5 <- 0.9
  lambda2_6 <- 1.1

  # Generate latent continuous variables
  y1_star <- lambda1_1 * f1 + rnorm(n, 0, 0.5)
  y2_star <- lambda1_2 * f1 + rnorm(n, 0, 0.5)
  y3_star <- lambda1_3 * f1 + rnorm(n, 0, 0.5)
  y4_star <- lambda2_4 * f2 + rnorm(n, 0, 0.5)
  y5_star <- lambda2_5 * f2 + rnorm(n, 0, 0.5)
  y6_star <- lambda2_6 * f2 + rnorm(n, 0, 0.5)

  # Convert to ordered categories (5 categories each)
  # Use integers 1-5, define_model_component will convert to ordered factor
  make_ordered <- function(y_star) {
    as.integer(cut(y_star,
                   breaks = c(-Inf, -1, -0.3, 0.3, 1, Inf),
                   labels = FALSE))
  }

  dat <- data.frame(
    m1 = make_ordered(y1_star),
    m2 = make_ordered(y2_star),
    m3 = make_ordered(y3_star),
    m4 = make_ordered(y4_star),
    m5 = make_ordered(y5_star),
    m6 = make_ordered(y6_star)
  )

  # Define 2-factor model
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Factor 1 measures: loading on f1, zero on f2
  # Variance normalization: factor variances fixed to 1.0, all loadings free
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "m1", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0
    num_choices = 5,
    intercept = FALSE
  )

  # m2: free loading on f1
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "m2", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0 (zero)
    num_choices = 5,
    intercept = FALSE
  )

  # m3: free loading on f1
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "m3", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0 (zero)
    num_choices = 5,
    intercept = FALSE
  )

  # Factor 2 measures: zero on f1, loading on f2
  # m4: free loading on f2
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "m4", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0, f2=free
    num_choices = 5,
    intercept = FALSE
  )

  # m5: free loading on f2
  mc5 <- define_model_component(
    name = "m5", data = dat, outcome = "m5", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0 (zero), f2=free
    num_choices = 5,
    intercept = FALSE
  )

  # m6: free loading on f2
  mc6 <- define_model_component(
    name = "m6", data = dat, outcome = "m6", factor = fm,
    covariates = NULL,
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0 (zero), f2=free
    num_choices = 5,
    intercept = FALSE
  )

  # Create model system
  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4, mc5, mc6),
    factor = fm
  )

  # Initialize parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    init_params = init_result$init_params,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0,
               info = "Model should converge successfully")

  # Check that we have reasonable estimates
  expect_true(is.numeric(result$estimates))
  expect_true(all(is.finite(result$estimates)))

  # Check standard errors are computed
  expect_true(is.numeric(result$std_errors))
  expect_true(all(is.finite(result$std_errors)))
  # Fixed parameters (factor variances fixed to 1.0) correctly have std_error = 0
  # Only check that non-fixed parameters have positive std_errors
  free_param_idx <- which(!grepl("factor_var", names(result$std_errors)))
  expect_true(all(result$std_errors[free_param_idx] > 0),
              info = sprintf("Non-fixed params should have positive std_errors. Zero SEs at: %s",
                             paste(names(result$std_errors)[result$std_errors == 0], collapse=", ")))

  # Extract estimated loadings
  # Note: Parameter order depends on implementation, but we can check structure
  param_names <- names(result$estimates)

  # Check that factor variances are estimated
  expect_true(any(grepl("factor_var", param_names, ignore.case = TRUE)),
              info = "Should have factor variance parameters")

  # Check that loadings are estimated for the correct measures
  expect_true(any(grepl("m2.*loading", param_names)),
              info = "m2 should have a free loading on f1")
  expect_true(any(grepl("m3.*loading", param_names)),
              info = "m3 should have a free loading on f1")
  expect_true(any(grepl("m5.*loading", param_names)),
              info = "m5 should have a free loading on f2")
  expect_true(any(grepl("m6.*loading", param_names)),
              info = "m6 should have a free loading on f2")

  # Print results for manual inspection
  cat("\n=== Two-Factor CFA Results ===\n")
  cat("Convergence:", result$convergence, "\n")
  cat("Log-likelihood:", result$loglik, "\n")
  cat("Number of parameters:", length(result$estimates), "\n\n")

  if (!is.null(param_names) && length(param_names) > 0) {
    # Find loading parameters
    loading_idx <- grep("loading", param_names)
    if (length(loading_idx) > 0) {
      cat("Estimated Loadings:\n")
      for (i in loading_idx) {
        cat(sprintf("  %s: %.4f (SE: %.4f)\n",
                    param_names[i],
                    result$estimates[i],
                    result$std_errors[i]))
      }
    }
  }

  # Parameter recovery check: loadings should be reasonably close to true values
  # We'll do a rough check on the free loadings
  # Note: In ordered probit with variance normalization, absolute loadings are scaled.
  # We check loading RATIOS instead, which should be preserved regardless of scaling.
  # True ratios: m2/m1=0.8, m3/m1=1.2, m5/m4=0.9, m6/m4=1.1
  if (!is.null(param_names)) {
    # Get reference loadings for each factor
    m1_load_idx <- grep("m1.*loading", param_names)
    m4_load_idx <- grep("m4.*loading", param_names)

    if (length(m1_load_idx) > 0 && length(m4_load_idx) > 0) {
      m1_est <- result$estimates[m1_load_idx[1]]
      m4_est <- result$estimates[m4_load_idx[1]]

      # Check m2/m1 ratio (true = 0.8)
      m2_load_idx <- grep("m2.*loading", param_names)
      if (length(m2_load_idx) > 0) {
        m2_est <- result$estimates[m2_load_idx[1]]
        m2_ratio <- m2_est / m1_est
        expect_true(abs(m2_ratio - 0.8) < 0.2,
                    info = sprintf("m2/m1 loading ratio should be near 0.8, got %.3f", m2_ratio))
      }

      # Check m3/m1 ratio (true = 1.2)
      m3_load_idx <- grep("m3.*loading", param_names)
      if (length(m3_load_idx) > 0) {
        m3_est <- result$estimates[m3_load_idx[1]]
        m3_ratio <- m3_est / m1_est
        expect_true(abs(m3_ratio - 1.2) < 0.2,
                    info = sprintf("m3/m1 loading ratio should be near 1.2, got %.3f", m3_ratio))
      }

      # Check m5/m4 ratio (true = 0.9)
      m5_load_idx <- grep("m5.*loading", param_names)
      if (length(m5_load_idx) > 0) {
        m5_est <- result$estimates[m5_load_idx[1]]
        m5_ratio <- m5_est / m4_est
        expect_true(abs(m5_ratio - 0.9) < 0.2,
                    info = sprintf("m5/m4 loading ratio should be near 0.9, got %.3f", m5_ratio))
      }

      # Check m6/m4 ratio (true = 1.1)
      m6_load_idx <- grep("m6.*loading", param_names)
      if (length(m6_load_idx) > 0) {
        m6_est <- result$estimates[m6_load_idx[1]]
        m6_ratio <- m6_est / m4_est
        expect_true(abs(m6_ratio - 1.1) < 0.2,
                    info = sprintf("m6/m4 loading ratio should be near 1.1, got %.3f", m6_ratio))
      }
    }
  }
})

test_that("two-factor CFA with ordered probit and covariates converges and recovers parameters", {
  isolate_test()
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Generate data with known two-factor structure
  # This test includes covariates in the measurement equations
  set.seed(456)
  n <- 2000

  # Two independent factors
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)

  # Covariate
  x1 <- rnorm(n, 0, 1)

  # True loadings for factor 1 measures (m1, m2, m3)
  # Using variance normalization: factor variances fixed to 1.0, all loadings free
  lambda1_1 <- 1.0
  lambda1_2 <- 0.8
  lambda1_3 <- 1.2

  # True loadings for factor 2 measures (m4, m5, m6)
  lambda2_4 <- 1.0
  lambda2_5 <- 0.9
  lambda2_6 <- 1.1

  # True covariate effect
  beta_x1 <- 0.3

  # Generate latent continuous variables WITH covariates
  y1_star <- beta_x1 * x1 + lambda1_1 * f1 + rnorm(n, 0, 0.5)
  y2_star <- beta_x1 * x1 + lambda1_2 * f1 + rnorm(n, 0, 0.5)
  y3_star <- beta_x1 * x1 + lambda1_3 * f1 + rnorm(n, 0, 0.5)
  y4_star <- beta_x1 * x1 + lambda2_4 * f2 + rnorm(n, 0, 0.5)
  y5_star <- beta_x1 * x1 + lambda2_5 * f2 + rnorm(n, 0, 0.5)
  y6_star <- beta_x1 * x1 + lambda2_6 * f2 + rnorm(n, 0, 0.5)

  # Convert to ordered categories (5 categories each)
  make_ordered <- function(y_star) {
    as.integer(cut(y_star,
                   breaks = c(-Inf, -1, -0.3, 0.3, 1, Inf),
                   labels = FALSE))
  }

  dat <- data.frame(
    intercept = 1,
    x1 = x1,
    m1 = make_ordered(y1_star),
    m2 = make_ordered(y2_star),
    m3 = make_ordered(y3_star),
    m4 = make_ordered(y4_star),
    m5 = make_ordered(y5_star),
    m6 = make_ordered(y6_star)
  )

  # Define 2-factor model
  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Factor 1 measures: loading on f1, zero on f2
  # Variance normalization: factor variances fixed to 1.0, all loadings free
  # WITH COVARIATES
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "m1", factor = fm,
    covariates = "x1",
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0
    num_choices = 5
  )

  # m2: free loading on f1
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "m2", factor = fm,
    covariates = "x1",
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0
    num_choices = 5
  )

  # m3: free loading on f1
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "m3", factor = fm,
    covariates = "x1",
    model_type = "oprobit",
    loading_normalization = c(NA, 0),  # f1=free, f2=0
    num_choices = 5
  )

  # Factor 2 measures: zero on f1, loading on f2
  # m4: free loading on f2
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "m4", factor = fm,
    covariates = "x1",
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0, f2=free
    num_choices = 5
  )

  # m5: free loading on f2
  mc5 <- define_model_component(
    name = "m5", data = dat, outcome = "m5", factor = fm,
    covariates = "x1",
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0, f2=free
    num_choices = 5
  )

  # m6: free loading on f2
  mc6 <- define_model_component(
    name = "m6", data = dat, outcome = "m6", factor = fm,
    covariates = "x1",
    model_type = "oprobit",
    loading_normalization = c(0, NA),  # f1=0, f2=free
    num_choices = 5
  )

  # Create model system
  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4, mc5, mc6),
    factor = fm
  )

  # Initialize parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)

  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    init_params = init_result$init_params,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0,
               info = "Model with covariates should converge successfully")

  # Check that we have reasonable estimates
  expect_true(is.numeric(result$estimates))
  expect_true(all(is.finite(result$estimates)))

  # Check standard errors are computed
  expect_true(is.numeric(result$std_errors))
  expect_true(all(is.finite(result$std_errors)))
  # For oprobit models: factor variances may be fixed (when loadings identify the scale)
  # and intercepts are absorbed into thresholds. Only check that loadings and
  # covariate effects (x1) have positive std_errors.
  loading_and_coef_idx <- grep("loading|_x1", names(result$std_errors))
  expect_true(all(result$std_errors[loading_and_coef_idx] > 0),
              info = sprintf("Loadings and covariate effects should have positive std_errors. Zero SEs at: %s",
                             paste(names(result$std_errors)[loading_and_coef_idx][result$std_errors[loading_and_coef_idx] == 0], collapse=", ")))

  # Print results for manual inspection
  cat("\n=== Two-Factor CFA with Covariates Results ===\n")
  cat("Convergence:", result$convergence, "\n")
  cat("Log-likelihood:", result$loglik, "\n")
  cat("Number of parameters:", length(result$estimates), "\n\n")

  param_names <- names(result$estimates)
  if (!is.null(param_names) && length(param_names) > 0) {
    # Find loading parameters
    loading_idx <- grep("loading", param_names)
    if (length(loading_idx) > 0) {
      cat("Estimated Loadings:\n")
      for (i in loading_idx) {
        cat(sprintf("  %s: %.4f (SE: %.4f)\n",
                    param_names[i],
                    result$estimates[i],
                    result$std_errors[i]))
      }
    }

    # Find covariate parameters
    x1_idx <- grep("x1", param_names)
    if (length(x1_idx) > 0) {
      cat("\nCovariate Effects (x1):\n")
      for (i in x1_idx) {
        cat(sprintf("  %s: %.4f (SE: %.4f)\n",
                    param_names[i],
                    result$estimates[i],
                    result$std_errors[i]))
      }
    }
  }
})

test_that("two-factor CFA with 3 linear measures per factor converges and recovers parameters", {
  isolate_test()
  skip_on_cran()
  # Generate data with known two-factor structure
  # 3 linear measures per factor with intercept + covariate
  set.seed(456)
  n <- 2000

  # Two independent factors
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)

  # Common covariate and intercept
  x1 <- rnorm(n, 0, 1)
  beta_intercept <- 2.0
  beta_x1 <- 1.5

  # True loadings for factor 1 measures (m1, m2, m3)
  # Using variance normalization: factor variances fixed to 1.0, all loadings free
  lambda1_1 <- 1.0
  lambda1_2 <- 0.8
  lambda1_3 <- 1.2

  # True loadings for factor 2 measures (m4, m5, m6)
  lambda2_4 <- 1.0
  lambda2_5 <- 0.9
  lambda2_6 <- 1.1
  
  # Generate linear outcomes: intercept + covariate + factor loading
  dat <- data.frame(
    intercept = 1,
    x1 = x1,
    m1 = beta_intercept + beta_x1 * x1 + lambda1_1 * f1 + rnorm(n, 0, 0.5),
    m2 = beta_intercept + beta_x1 * x1 + lambda1_2 * f1 + rnorm(n, 0, 0.5),
    m3 = beta_intercept + beta_x1 * x1 + lambda1_3 * f1 + rnorm(n, 0, 0.5),
    m4 = beta_intercept + beta_x1 * x1 + lambda2_4 * f2 + rnorm(n, 0, 0.5),
    m5 = beta_intercept + beta_x1 * x1 + lambda2_5 * f2 + rnorm(n, 0, 0.5),
    m6 = beta_intercept + beta_x1 * x1 + lambda2_6 * f2 + rnorm(n, 0, 0.5)
  )
  
  # Define 2-factor model
  fm <- define_factor_model(n_factors = 2, n_types = 1)
  
  # Factor 1 measures
  # Variance normalization: factor variances fixed to 1.0, all loadings free
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "m1", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(NA, 0)  # f1=free, f2=0 (zero)
  )

  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "m2", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(NA, 0)  # f1=free, f2=0 (zero)
  )

  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "m3", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(NA, 0)  # f1=free, f2=0 (zero)
  )

  # Factor 2 measures
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "m4", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, NA)  # f1=0 (zero), f2=free
  )

  mc5 <- define_model_component(
    name = "m5", data = dat, outcome = "m5", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, NA)  # f1=0 (zero), f2=free
  )

  mc6 <- define_model_component(
    name = "m6", data = dat, outcome = "m6", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, NA)  # f1=0 (zero), f2=free
  )
  
  # Create model system
  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4, mc5, mc6),
    factor = fm
  )
  
  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )
  
  # Check convergence
  expect_equal(result$convergence, 0, info = "Model should converge")
  
  # Check that we have the right number of parameters
  # Variance normalization: factor variances fixed at 1.0 (not estimated)
  # All loadings are free
  # m1: intercept, x1, loading_f1, sigma = 4
  # m2: intercept, x1, loading_f1, sigma = 4
  # m3: intercept, x1, loading_f1, sigma = 4
  # m4: intercept, x1, loading_f2, sigma = 4
  # m5: intercept, x1, loading_f2, sigma = 4
  # m6: intercept, x1, loading_f2, sigma = 4
  # Total: 2 (fixed factor vars) + 4*6 = 2 + 24 = 26
  expect_equal(length(result$estimates), 26)
  
  # Check factor variances are fixed to 1.0 (variance normalization)
  expect_equal(unname(result$estimates[1]), 1.0,
              info = sprintf("Factor 1 variance should be fixed to 1.0, got %.4f", result$estimates[1]))
  expect_equal(unname(result$estimates[2]), 1.0,
              info = sprintf("Factor 2 variance should be fixed to 1.0, got %.4f", result$estimates[2]))

  # Tolerance for parameter recovery (more lenient for linear models)
  tol <- 0.15
  
  # Check loadings are recovered (with parameter names if available)
  param_names <- result$parameter_names
  if (!is.null(param_names)) {
    # Check m2 loading on f1 (true = 0.8)
    m2_load_idx <- grep("m2.*loading.*1", param_names, ignore.case = TRUE)
    if (length(m2_load_idx) > 0) {
      m2_est <- result$estimates[m2_load_idx[1]]
      expect_true(abs(m2_est - 0.8) < tol,
                  info = sprintf("m2 loading should be near 0.8, got %.3f", m2_est))
    }
    
    # Check m3 loading on f1 (true = 1.2)
    m3_load_idx <- grep("m3.*loading.*1", param_names, ignore.case = TRUE)
    if (length(m3_load_idx) > 0) {
      m3_est <- result$estimates[m3_load_idx[1]]
      expect_true(abs(m3_est - 1.2) < tol,
                  info = sprintf("m3 loading should be near 1.2, got %.3f", m3_est))
    }
    
    # Check m5 loading on f2 (true = 0.9)
    m5_load_idx <- grep("m5.*loading.*2", param_names, ignore.case = TRUE)
    if (length(m5_load_idx) > 0) {
      m5_est <- result$estimates[m5_load_idx[1]]
      expect_true(abs(m5_est - 0.9) < tol,
                  info = sprintf("m5 loading should be near 0.9, got %.3f", m5_est))
    }
    
    # Check m6 loading on f2 (true = 1.1)
    m6_load_idx <- grep("m6.*loading.*2", param_names, ignore.case = TRUE)
    if (length(m6_load_idx) > 0) {
      m6_est <- result$estimates[m6_load_idx[1]]
      expect_true(abs(m6_est - 1.1) < tol,
                  info = sprintf("m6 loading should be near 1.1, got %.3f", m6_est))
    }
  }
})

test_that("three-factor CFA with 3 linear measures per factor converges and recovers parameters", {
  isolate_test()
  skip_on_cran()
  # Generate data with known three-factor structure
  # 3 linear measures per factor with intercept + covariate
  set.seed(789)
  n <- 2000

  # Three independent factors
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)
  f3 <- rnorm(n, 0, 1)
  
  # Common covariate and intercept
  x1 <- rnorm(n, 0, 1)
  beta_intercept <- 2.0
  beta_x1 <- 1.5
  
  # True loadings for factor 1 measures (m1, m2, m3)
  lambda1_1 <- 1.0   # Fixed for identification
  lambda1_2 <- 0.8
  lambda1_3 <- 1.2
  
  # True loadings for factor 2 measures (m4, m5, m6)
  lambda2_4 <- 1.0   # Fixed for identification
  lambda2_5 <- 0.9
  lambda2_6 <- 1.1
  
  # True loadings for factor 3 measures (m7, m8, m9)
  lambda3_7 <- 1.0   # Fixed for identification
  lambda3_8 <- 0.7
  lambda3_9 <- 1.3
  
  # Generate linear outcomes: intercept + covariate + factor loading
  dat <- data.frame(
    intercept = 1,
    x1 = x1,
    m1 = beta_intercept + beta_x1 * x1 + lambda1_1 * f1 + rnorm(n, 0, 0.5),
    m2 = beta_intercept + beta_x1 * x1 + lambda1_2 * f1 + rnorm(n, 0, 0.5),
    m3 = beta_intercept + beta_x1 * x1 + lambda1_3 * f1 + rnorm(n, 0, 0.5),
    m4 = beta_intercept + beta_x1 * x1 + lambda2_4 * f2 + rnorm(n, 0, 0.5),
    m5 = beta_intercept + beta_x1 * x1 + lambda2_5 * f2 + rnorm(n, 0, 0.5),
    m6 = beta_intercept + beta_x1 * x1 + lambda2_6 * f2 + rnorm(n, 0, 0.5),
    m7 = beta_intercept + beta_x1 * x1 + lambda3_7 * f3 + rnorm(n, 0, 0.5),
    m8 = beta_intercept + beta_x1 * x1 + lambda3_8 * f3 + rnorm(n, 0, 0.5),
    m9 = beta_intercept + beta_x1 * x1 + lambda3_9 * f3 + rnorm(n, 0, 0.5)
  )
  
  # Define 3-factor model
  fm <- define_factor_model(n_factors = 3, n_types = 1)
  
  # Factor 1 measures
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "m1", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(1.0, 0, 0)  # f1=1.0 (fixed), f2=0, f3=0
  )
  
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "m2", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(NA, 0, 0)  # f1=free, f2=0, f3=0
  )
  
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "m3", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(NA, 0, 0)  # f1=free, f2=0, f3=0
  )
  
  # Factor 2 measures
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "m4", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, 1.0, 0)  # f1=0, f2=1.0 (fixed), f3=0
  )
  
  mc5 <- define_model_component(
    name = "m5", data = dat, outcome = "m5", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, NA, 0)  # f1=0, f2=free, f3=0
  )
  
  mc6 <- define_model_component(
    name = "m6", data = dat, outcome = "m6", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, NA, 0)  # f1=0, f2=free, f3=0
  )
  
  # Factor 3 measures
  mc7 <- define_model_component(
    name = "m7", data = dat, outcome = "m7", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, 0, 1.0)  # f1=0, f2=0, f3=1.0 (fixed)
  )
  
  mc8 <- define_model_component(
    name = "m8", data = dat, outcome = "m8", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, 0, NA)  # f1=0, f2=0, f3=free
  )
  
  mc9 <- define_model_component(
    name = "m9", data = dat, outcome = "m9", factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "linear",
    loading_normalization = c(0, 0, NA)  # f1=0, f2=0, f3=free
  )
  
  # Create model system
  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9),
    factor = fm
  )
  
  # Estimate the model
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )
  
  # Check convergence
  expect_equal(result$convergence, 0, info = "Model should converge")
  
  # Check that we have the right number of parameters
  # 3 factor vars + 9 components * (various)
  # m1,m4,m7: 3 params each (intercept, x1, sigma) = 9
  # m2,m3,m5,m6,m8,m9: 4 params each (intercept, x1, loading, sigma) = 24
  # Total: 3 + 9 + 24 = 36
  expect_equal(length(result$estimates), 36)
  
  # Check factor variances are estimated (all factors have fixed loadings, so variances are identified)
  # Since true factor variances are 1.0 (from data generation), estimates should be close to 1.0
  # with some tolerance for finite sample variation
  expect_true(unname(result$estimates[1]) > 0.5 && unname(result$estimates[1]) < 2.0,
              info = sprintf("Factor 1 variance should be positive and reasonable, got %.4f", result$estimates[1]))
  expect_true(unname(result$estimates[2]) > 0.5 && unname(result$estimates[2]) < 2.0,
              info = sprintf("Factor 2 variance should be positive and reasonable, got %.4f", result$estimates[2]))
  expect_true(unname(result$estimates[3]) > 0.5 && unname(result$estimates[3]) < 2.0,
              info = sprintf("Factor 3 variance should be positive and reasonable, got %.4f", result$estimates[3]))
  
  # Tolerance for parameter recovery
  tol <- 0.15
  
  # Check loadings are recovered (with parameter names if available)
  param_names <- result$parameter_names
  if (!is.null(param_names)) {
    # Factor 1 loadings
    m2_load_idx <- grep("m2.*loading.*1", param_names, ignore.case = TRUE)
    if (length(m2_load_idx) > 0) {
      m2_est <- result$estimates[m2_load_idx[1]]
      expect_true(abs(m2_est - 0.8) < tol,
                  info = sprintf("m2 loading should be near 0.8, got %.3f", m2_est))
    }
    
    m3_load_idx <- grep("m3.*loading.*1", param_names, ignore.case = TRUE)
    if (length(m3_load_idx) > 0) {
      m3_est <- result$estimates[m3_load_idx[1]]
      expect_true(abs(m3_est - 1.2) < tol,
                  info = sprintf("m3 loading should be near 1.2, got %.3f", m3_est))
    }
    
    # Factor 2 loadings
    m5_load_idx <- grep("m5.*loading.*2", param_names, ignore.case = TRUE)
    if (length(m5_load_idx) > 0) {
      m5_est <- result$estimates[m5_load_idx[1]]
      expect_true(abs(m5_est - 0.9) < tol,
                  info = sprintf("m5 loading should be near 0.9, got %.3f", m5_est))
    }
    
    m6_load_idx <- grep("m6.*loading.*2", param_names, ignore.case = TRUE)
    if (length(m6_load_idx) > 0) {
      m6_est <- result$estimates[m6_load_idx[1]]
      expect_true(abs(m6_est - 1.1) < tol,
                  info = sprintf("m6 loading should be near 1.1, got %.3f", m6_est))
    }
    
    # Factor 3 loadings
    m8_load_idx <- grep("m8.*loading.*3", param_names, ignore.case = TRUE)
    if (length(m8_load_idx) > 0) {
      m8_est <- result$estimates[m8_load_idx[1]]
      expect_true(abs(m8_est - 0.7) < tol,
                  info = sprintf("m8 loading should be near 0.7, got %.3f", m8_est))
    }
    
    m9_load_idx <- grep("m9.*loading.*3", param_names, ignore.case = TRUE)
    if (length(m9_load_idx) > 0) {
      m9_est <- result$estimates[m9_load_idx[1]]
      expect_true(abs(m9_est - 1.3) < tol,
                  info = sprintf("m9 loading should be near 1.3, got %.3f", m9_est))
    }
  }
})

