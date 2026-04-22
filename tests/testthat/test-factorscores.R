# Test factor score estimation
# This test validates factor score recovery by comparing estimated scores
# to the true simulated factors.

test_that("Factor score estimation recovers true factors in linear model", {
  skip_on_cran()

  # Generate data with known factor structure
  set.seed(42)
  n <- 200

  # True factor values (save for comparison)
  true_factors <- rnorm(n, 0, 1)

  # Generate 3 linear measures with known loadings
  lambda1 <- 1.0  # Normalized to 1 for identification
  lambda2 <- 0.8
  lambda3 <- 1.2
  sigma <- 0.5    # Error SD for each measure

  y1 <- lambda1 * true_factors + rnorm(n, 0, sigma)
  y2 <- lambda2 * true_factors + rnorm(n, 0, sigma)
  y3 <- lambda3 * true_factors + rnorm(n, 0, sigma)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, intercept = 1)

  # Define 1-factor model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Define measurement components
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1  # Fixed at 1 for identification
  )

  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear"
  )

  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear"
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3),
    factor = fm
  )

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  # Estimate the model
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  expect_true(result$convergence == 0,
              info = "Model estimation should converge")

  # Estimate factor scores
  factor_scores <- estimate_factorscores_rcpp(
    result = result,
    data = dat,
    verbose = FALSE
  )

  expect_s3_class(factor_scores, "factorana_factorscores")
  expect_equal(nrow(factor_scores), n)

  # Get converged observations
  converged_idx <- which(factor_scores$converged)
  n_converged <- length(converged_idx)

  # Validation statistics
  estimated <- factor_scores$factor_1[converged_idx]
  true_vals <- true_factors[converged_idx]
  std_errors <- factor_scores$se_factor_1[converged_idx]

  # 1. Correlation: should be reasonably high (> 0.5 for 3 measures with noise)
  correlation <- cor(estimated, true_vals)
  cat("\n  Correlation between estimated and true factors:", round(correlation, 3), "\n")
  expect_gt(correlation, 0.5,
            label = "Correlation between estimated and true factors")

  # 2. Mean absolute standardized error: |error|/SE
  #    If SEs are well-calibrated, this should be ~0.8 (E[|Z|] for standard normal)
  abs_std_errors <- abs(estimated - true_vals) / std_errors
  mean_abs_std_error <- mean(abs_std_errors, na.rm = TRUE)
  cat("  Mean |error|/SE:", round(mean_abs_std_error, 3), "\n")

  # 3. Coverage: % of true values within 95% CI
  lower <- estimated - 1.96 * std_errors
  upper <- estimated + 1.96 * std_errors
  covered <- (true_vals >= lower) & (true_vals <= upper)
  coverage_rate <- mean(covered, na.rm = TRUE)
  cat("  95% CI coverage rate:", round(coverage_rate * 100, 1), "%\n")
  # Coverage is limited with only 3 measures; expect reasonable but not 95%
  expect_gt(coverage_rate, 0.50,
            label = "95% CI coverage")

  # 4. Bias: mean(estimated - true) should be near 0
  bias <- mean(estimated - true_vals, na.rm = TRUE)
  cat("  Bias (mean error):", round(bias, 3), "\n")
  expect_lt(abs(bias), 0.3,
            label = "Absolute bias")

  # 5. Convergence rate - should be 100% for all observations
  conv_rate <- n_converged / n
  cat("  Convergence rate:", round(conv_rate * 100, 1), "%\n")
  expect_equal(conv_rate, 1.0,
               label = "Convergence rate (should be 100%)")
})


test_that("Factor score estimation works with probit model", {
  skip_on_cran()

  # Generate data with known factor structure
  set.seed(123)
  n <- 300

  # True factor values
  true_factors <- rnorm(n, 0, 1)

  # Generate 3 binary measures with known loadings
  lambda1 <- 1.0
  lambda2 <- 0.8
  lambda3 <- 1.2

  # Latent variables
  y1_star <- lambda1 * true_factors + rnorm(n, 0, 1)
  y2_star <- lambda2 * true_factors + rnorm(n, 0, 1)
  y3_star <- lambda3 * true_factors + rnorm(n, 0, 1)

  # Binary outcomes
  dat <- data.frame(
    y1 = as.integer(y1_star > 0),
    y2 = as.integer(y2_star > 0),
    y3 = as.integer(y3_star > 0),
    intercept = 1
  )

  # Define model
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "probit",
    loading_normalization = 1
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "probit"
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "probit"
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3),
    factor = fm
  )

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  # Estimate the model
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Estimate factor scores
  factor_scores <- estimate_factorscores_rcpp(
    result = result,
    data = dat,
    verbose = FALSE
  )

  # Validation
  converged_idx <- which(factor_scores$converged)
  estimated <- factor_scores$factor_1[converged_idx]
  true_vals <- true_factors[converged_idx]
  std_errors <- factor_scores$se_factor_1[converged_idx]

  correlation <- cor(estimated, true_vals)
  cat("\n  Probit model - Correlation:", round(correlation, 3), "\n")
  cat("  Probit model - Convergence rate:",
      round(length(converged_idx) / n * 100, 1), "%\n")

  # Check 95% CI coverage
  lower <- estimated - 1.96 * std_errors
  upper <- estimated + 1.96 * std_errors
  covered <- (true_vals >= lower) & (true_vals <= upper)
  coverage_rate <- mean(covered, na.rm = TRUE)
  cat("  Probit model - 95% CI coverage:", round(coverage_rate * 100, 1), "%\n")

  expect_gt(correlation, 0.4,
            label = "Probit model correlation")
  expect_equal(length(converged_idx) / n, 1.0,
               label = "Probit model convergence rate (should be 100%)")
  expect_gt(coverage_rate, 0.50,
            label = "Probit model coverage rate")
})


test_that("Factor score estimation handles two-factor model", {
  skip_on_cran()

  set.seed(456)
  n <- 300

  # Two independent factors
  f1 <- rnorm(n, 0, 1)
  f2 <- rnorm(n, 0, 1)

  sigma <- 0.5

  # Factor 1 measures
  y1 <- 1.0 * f1 + rnorm(n, 0, sigma)
  y2 <- 0.8 * f1 + rnorm(n, 0, sigma)

  # Factor 2 measures
  y3 <- 1.0 * f2 + rnorm(n, 0, sigma)
  y4 <- 0.9 * f2 + rnorm(n, 0, sigma)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, n_types = 1)

  # Factor 1 measures (load on factor 1 only)
  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0)  # Fixed at 1 for f1, 0 for f2
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0)  # Free for f1, 0 for f2
  )

  # Factor 2 measures (load on factor 2 only)
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1)  # 0 for f1, fixed at 1 for f2
  )
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA)  # 0 for f1, free for f2
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3, mc4),
    factor = fm
  )

  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  # Estimate
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Factor scores
  factor_scores <- estimate_factorscores_rcpp(
    result = result,
    data = dat,
    verbose = FALSE
  )

  # obs_id, 2 factors, 2 SEs, converged, log_posterior = 7 columns
  expect_equal(ncol(factor_scores), 7)

  converged_idx <- which(factor_scores$converged)
  n_converged <- length(converged_idx)

  # Convergence rate - should be 100%
  conv_rate <- n_converged / n
  cat("\n  Two-factor model - Convergence rate:", round(conv_rate * 100, 1), "%\n")
  expect_equal(conv_rate, 1.0,
               label = "Two-factor model convergence rate (should be 100%)")

  # Factor 1 validation
  est_f1 <- factor_scores$factor_1[converged_idx]
  true_f1 <- f1[converged_idx]
  se_f1 <- factor_scores$se_factor_1[converged_idx]

  cor_f1 <- cor(est_f1, true_f1)
  cat("  Two-factor model - Factor 1 correlation:", round(cor_f1, 3), "\n")

  # Check 95% CI coverage for factor 1
  lower_f1 <- est_f1 - 1.96 * se_f1
  upper_f1 <- est_f1 + 1.96 * se_f1
  covered_f1 <- (true_f1 >= lower_f1) & (true_f1 <= upper_f1)
  coverage_f1 <- mean(covered_f1, na.rm = TRUE)
  cat("  Two-factor model - Factor 1 95% CI coverage:", round(coverage_f1 * 100, 1), "%\n")

  # Factor 2 validation
  est_f2 <- factor_scores$factor_2[converged_idx]
  true_f2 <- f2[converged_idx]
  se_f2 <- factor_scores$se_factor_2[converged_idx]

  cor_f2 <- cor(est_f2, true_f2)
  cat("  Two-factor model - Factor 2 correlation:", round(cor_f2, 3), "\n")

  # Check 95% CI coverage for factor 2
  lower_f2 <- est_f2 - 1.96 * se_f2
  upper_f2 <- est_f2 + 1.96 * se_f2
  covered_f2 <- (true_f2 >= lower_f2) & (true_f2 <= upper_f2)
  coverage_f2 <- mean(covered_f2, na.rm = TRUE)
  cat("  Two-factor model - Factor 2 95% CI coverage:", round(coverage_f2 * 100, 1), "%\n")

  # Cross-correlations should be low (factors are orthogonal)
  cross_cor_12 <- cor(est_f1, true_f2)
  cross_cor_21 <- cor(est_f2, true_f1)
  cat("  Cross-correlation (est_f1 vs true_f2):", round(cross_cor_12, 3), "\n")
  cat("  Cross-correlation (est_f2 vs true_f1):", round(cross_cor_21, 3), "\n")

  expect_gt(cor_f1, 0.4,
            label = "Factor 1 correlation")
  expect_gt(cor_f2, 0.4,
            label = "Factor 2 correlation")
  expect_gt(coverage_f1, 0.50,
            label = "Factor 1 95% CI coverage")
  expect_gt(coverage_f2, 0.50,
            label = "Factor 2 95% CI coverage")
})


test_that("Factor score estimation works with ordered probit model", {
  skip_on_cran()

  set.seed(789)
  n <- 300
  true_factors <- rnorm(n, 0, 1)

  # Generate ordinal outcomes with 4 categories
  lambda1 <- 1.0
  lambda2 <- 0.8
  lambda3 <- 1.2

  # Latent continuous variables
  y1_star <- lambda1 * true_factors + rnorm(n, 0, 1)
  y2_star <- lambda2 * true_factors + rnorm(n, 0, 1)
  y3_star <- lambda3 * true_factors + rnorm(n, 0, 1)

  # Convert to ordered categories (1, 2, 3, 4)
  cut_to_ordered <- function(y_star) {
    cuts <- c(-Inf, -0.5, 0.5, 1.5, Inf)
    as.integer(cut(y_star, breaks = cuts, labels = FALSE))
  }

  # Create data frame with integer outcomes
  dat <- data.frame(
    y1 = cut_to_ordered(y1_star),
    y2 = cut_to_ordered(y2_star),
    y3 = cut_to_ordered(y3_star)
  )

  # For ordered probit, keep outcome as integer (1, 2, 3, 4)
  # The C++ code expects numeric values, not factors

  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = NULL, model_type = "oprobit",
    loading_normalization = 1, num_choices = 4,
    intercept = FALSE
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = NULL, model_type = "oprobit",
    num_choices = 4,
    intercept = FALSE
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = NULL, model_type = "oprobit",
    num_choices = 4,
    intercept = FALSE
  )

  ms <- define_model_system(
    components = list(mc1, mc2, mc3),
    factor = fm
  )

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  # Estimate the model
  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = ctrl,
    parallel = FALSE,
    optimizer = "nlminb",
    verbose = FALSE
  )

  # Estimate factor scores
  factor_scores <- estimate_factorscores_rcpp(
    result = result,
    data = dat,
    verbose = FALSE
  )

  # Validation
  converged_idx <- which(factor_scores$converged)
  estimated <- factor_scores$factor_1[converged_idx]
  true_vals <- true_factors[converged_idx]
  std_errors <- factor_scores$se_factor_1[converged_idx]

  correlation <- cor(estimated, true_vals)
  cat("\n  Ordered probit model - Correlation:", round(correlation, 3), "\n")
  cat("  Ordered probit model - Convergence rate:",
      round(length(converged_idx) / n * 100, 1), "%\n")

  # Check 95% CI coverage
  lower <- estimated - 1.96 * std_errors
  upper <- estimated + 1.96 * std_errors
  covered <- (true_vals >= lower) & (true_vals <= upper)
  coverage_rate <- mean(covered, na.rm = TRUE)
  cat("  Ordered probit model - 95% CI coverage:", round(coverage_rate * 100, 1), "%\n")

  expect_gt(correlation, 0.4,
            label = "Ordered probit model correlation")
  expect_equal(length(converged_idx) / n, 1.0,
               label = "Ordered probit model convergence rate (should be 100%)")
  expect_gt(coverage_rate, 0.50,
            label = "Ordered probit model coverage rate")
})

test_that("Zero-observation individuals get NA factor scores and SEs", {
  skip_on_cran()

  set.seed(99)
  n <- 50
  true_f <- rnorm(n)

  # 3 linear measures with an evaluation indicator.
  # First 5 observations have eval=0 for all items → no data.
  eval_ind <- rep(1, n)
  eval_ind[1:5] <- 0

  dat <- data.frame(
    intercept = 1,
    y1 = 1.0 * true_f + rnorm(n, 0, 0.5),
    y2 = 0.8 * true_f + rnorm(n, 0, 0.5),
    y3 = 0.9 * true_f + rnorm(n, 0, 0.5),
    eval = eval_ind
  )

  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = 1,
                                 evaluation_indicator = "eval")
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = NA_real_,
                                 evaluation_indicator = "eval")
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = NA_real_,
                                 evaluation_indicator = "eval")

  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb",
                                parallel = FALSE, verbose = FALSE)
  expect_equal(result$convergence, 0)

  # Estimate factor scores
  fscores <- estimate_factorscores_rcpp(result, dat, control = ctrl,
                                         parallel = FALSE, verbose = FALSE,
                                         include_prior = TRUE)

  # Zero-obs individuals (rows 1-5) should have NA scores and SEs
  for (i in 1:5) {
    expect_true(is.na(fscores$factor_1[i]),
                info = sprintf("Obs %d should have NA factor score", i))
    expect_true(is.na(fscores$se_factor_1[i]),
                info = sprintf("Obs %d should have NA SE", i))
  }

  # Observed individuals (rows 6-50) should have finite scores and SEs
  for (i in 6:50) {
    expect_true(is.finite(fscores$factor_1[i]),
                info = sprintf("Obs %d should have finite factor score", i))
    expect_true(is.finite(fscores$se_factor_1[i]),
                info = sprintf("Obs %d should have finite SE", i))
  }

  # SEs for observed individuals should be bounded by prior SD when include_prior=TRUE
  # (posterior precision >= prior precision for log-concave likelihoods)
  prior_sd <- sqrt(result$estimates["factor_var_1"])
  for (i in 6:50) {
    expect_lte(fscores$se_factor_1[i], prior_sd * 1.01,  # tiny tolerance for numerics
               label = sprintf("Obs %d SE (%.4f) vs prior SD (%.4f)",
                              i, fscores$se_factor_1[i], prior_sd))
  }
})
