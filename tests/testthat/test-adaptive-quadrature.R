# =============================================================================
# Test: Adaptive Quadrature for Two-Stage Estimation
# =============================================================================

test_that("Adaptive quadrature works for two-stage estimation", {
  skip_on_cran()

  set.seed(12345)
  n <- 400

  # ==========================================================================
  # Simulate Data
  # ==========================================================================
  # Design: 4 linear measures where y3 and y4 are missing for 25% of individuals.
  # For missing individuals, only y1 and y2 (with weak loadings and high noise)
  # are used, creating high uncertainty in factor scores.

  true_factor_var <- 1.0
  # Weaker loadings and higher noise to make factors harder to identify
  true_loadings <- c(0.6, 0.5, 0.9, 0.8)  # y1, y2 are weak; y3, y4 are strong
  true_sigmas <- c(0.9, 1.0, 0.4, 0.45)   # y1, y2 have high noise; y3, y4 are precise

  true_linear_beta <- 0.5
  true_linear_loading <- 0.6
  true_linear_sigma <- 0.4

  true_probit_beta <- 0.3
  true_probit_loading <- 0.5

  f <- rnorm(n, 0, sqrt(true_factor_var))
  x <- rnorm(n, 0, 1)

  y1 <- true_loadings[1] * f + rnorm(n, 0, true_sigmas[1])
  y2 <- true_loadings[2] * f + rnorm(n, 0, true_sigmas[2])
  y3 <- true_loadings[3] * f + rnorm(n, 0, true_sigmas[3])
  y4 <- true_loadings[4] * f + rnorm(n, 0, true_sigmas[4])

  missing_idx <- sample(1:n, size = floor(0.25 * n), replace = FALSE)
  eval_y3 <- rep(1, n)
  eval_y3[missing_idx] <- 0
  eval_y4 <- rep(1, n)
  eval_y4[missing_idx] <- 0

  y_linear <- true_linear_beta * x + true_linear_loading * f +
              rnorm(n, 0, true_linear_sigma)

  latent_probit <- true_probit_beta * x + true_probit_loading * f + rnorm(n, 0, 1)
  y_probit <- as.integer(latent_probit > 0)

  dat <- data.frame(
    y1 = y1, y2 = y2, y3 = y3, y4 = y4,
    eval_y3 = eval_y3, eval_y4 = eval_y4,
    x = x,
    y_linear = y_linear,
    y_probit = y_probit,
    intercept = 1
  )

  # ==========================================================================
  # Stage 1: Estimate Measurement System
  # ==========================================================================

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component(
    name = "m1", data = dat, outcome = "y1", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = 1
  )
  mc2 <- define_model_component(
    name = "m2", data = dat, outcome = "y2", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = NA_real_
  )
  mc3 <- define_model_component(
    name = "m3", data = dat, outcome = "y3", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = NA_real_,
    evaluation_indicator = "eval_y3"
  )
  mc4 <- define_model_component(
    name = "m4", data = dat, outcome = "y4", factor = fm,
    covariates = "intercept", model_type = "linear", loading_normalization = NA_real_,
    evaluation_indicator = "eval_y4"
  )

  ms_stage1 <- define_model_system(
    components = list(mc1, mc2, mc3, mc4),
    factor = fm
  )

  control_stage1 <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    model_system = ms_stage1,
    data = dat,
    control = control_stage1,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  expect_equal(result_stage1$convergence, 0, label = "Stage 1 converged")

  # ==========================================================================
  # Compute Factor Scores and Check SE Pattern
  # ==========================================================================

  fscores <- estimate_factorscores_rcpp(
    result = result_stage1,
    data = dat,
    control = control_stage1,
    verbose = FALSE
  )

  factor_scores_mat <- as.matrix(fscores[, "factor_1", drop = FALSE])
  factor_ses_mat <- as.matrix(fscores[, "se_factor_1", drop = FALSE])
  factor_var <- as.numeric(result_stage1$estimates["factor_var_1"])

  se_complete <- factor_ses_mat[-missing_idx, 1]
  se_missing <- factor_ses_mat[missing_idx, 1]
  se_ratio <- mean(se_missing) / mean(se_complete)

  expect_true(se_ratio > 1.2,
              label = sprintf("SE ratio (missing/complete) = %.2f > 1.2", se_ratio))

  # ==========================================================================
  # Check Adaptive Quadrature Point Distribution
  # ==========================================================================

  # Use a lower threshold to get more quadrature points and better accuracy
  # Higher threshold = fewer points = faster but less accurate
  threshold <- 0.15
  max_quad <- 8

  compute_nquad <- function(se, fvar, thresh, max_q) {
    nq <- 1 + 2 * floor(se / fvar / thresh)
    nq <- pmin(nq, max_q)
    nq[se > sqrt(fvar)] <- max_q
    return(nq)
  }

  nquad_obs <- compute_nquad(factor_ses_mat[, 1], factor_var, threshold, max_quad)
  avg_nquad <- mean(nquad_obs)

  # Verify that adaptive quadrature reduces computational cost
  # (some observations should use fewer than max quadrature points)
  expect_true(avg_nquad < max_quad,
              label = sprintf("Average quad points = %.2f < %d", avg_nquad, max_quad))

  # Verify that observations with missing data need more quadrature points on average
  avg_nquad_missing <- mean(nquad_obs[missing_idx])
  avg_nquad_complete <- mean(nquad_obs[-missing_idx])
  expect_true(avg_nquad_missing >= avg_nquad_complete,
              label = sprintf("Missing obs use %.2f pts, complete use %.2f pts",
                              avg_nquad_missing, avg_nquad_complete))

  # ==========================================================================
  # Stage 2: Define Model with All Model Types
  # ==========================================================================

  mc_linear <- define_model_component(
    name = "outcome_linear", data = dat, outcome = "y_linear", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = NA_real_
  )

  mc_probit <- define_model_component(
    name = "outcome_probit", data = dat, outcome = "y_probit", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = NA_real_
  )

  ms_stage2 <- define_model_system(
    components = list(mc_linear, mc_probit),
    factor = fm,
    previous_stage = result_stage1
  )

  # ==========================================================================
  # Test 1: Standard Quadrature (Baseline)
  # ==========================================================================

  control_standard <- define_estimation_control(
    n_quad_points = max_quad,
    num_cores = 1,
    adaptive_integration = FALSE
  )

  result_standard <- estimate_model_rcpp(
    model_system = ms_stage2,
    data = dat,
    control = control_standard,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  expect_equal(result_standard$convergence, 0, label = "Standard quadrature converged")

  # ==========================================================================
  # Test 2: Adaptive Quadrature (Single Worker)
  # ==========================================================================

  control_adaptive <- define_estimation_control(
    n_quad_points = max_quad,
    num_cores = 1,
    adaptive_integration = TRUE,
    adapt_int_thresh = threshold
  )

  factor_vars_vec <- c(factor_var_1 = factor_var)

  result_adaptive <- estimate_model_rcpp(
    model_system = ms_stage2,
    data = dat,
    control = control_adaptive,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE,
    factor_scores = factor_scores_mat,
    factor_ses = factor_ses_mat,
    factor_vars = factor_vars_vec
  )

  expect_equal(result_adaptive$convergence, 0, label = "Adaptive quadrature converged")

  # ==========================================================================
  # Test 3: Adaptive Quadrature (Parallel)
  # ==========================================================================

  control_parallel <- define_estimation_control(
    n_quad_points = max_quad,
    num_cores = 2,
    adaptive_integration = TRUE,
    adapt_int_thresh = threshold
  )

  result_parallel <- estimate_model_rcpp(
    model_system = ms_stage2,
    data = dat,
    control = control_parallel,
    optimizer = "nlminb",
    parallel = TRUE,
    verbose = FALSE,
    factor_scores = factor_scores_mat,
    factor_ses = factor_ses_mat,
    factor_vars = factor_vars_vec
  )

  expect_equal(result_parallel$convergence, 0, label = "Parallel adaptive converged")

  # ==========================================================================
  # Compare Results Across Methods
  # ==========================================================================

  n_stage1_params <- length(result_stage1$estimates)
  stage2_param_idx <- (n_stage1_params + 1):length(result_standard$estimates)

  est_standard <- result_standard$estimates[stage2_param_idx]
  est_adaptive <- result_adaptive$estimates[stage2_param_idx]
  est_parallel <- result_parallel$estimates[stage2_param_idx]

  # Verify adaptive quadrature is actually being used (LL values should differ)
  # If factor scores weren't being transmitted, adaptive would equal standard
  ll_diff_adaptive <- abs(result_adaptive$loglik - result_standard$loglik)
  ll_diff_parallel <- abs(result_parallel$loglik - result_standard$loglik)

  expect_true(ll_diff_adaptive > 1.0,
              label = sprintf("Adaptive LL differs from standard (diff = %.2f > 1.0)", ll_diff_adaptive))
  expect_true(ll_diff_parallel > 1.0,
              label = sprintf("Parallel LL differs from standard (diff = %.2f > 1.0)", ll_diff_parallel))

  # Log-likelihoods: adaptive quadrature is an approximation, so we allow larger
  # differences. The key validation is parameter estimates, not raw LL values.
  # With threshold 0.15, typical LL difference is 5-15% of |LL|.
  ll_tol <- abs(result_standard$loglik) * 0.15
  expect_true(ll_diff_adaptive < ll_tol,
              label = sprintf("Adaptive LL diff = %.4f < %.4f (15%% of |LL|)", ll_diff_adaptive, ll_tol))
  expect_true(ll_diff_parallel < ll_tol,
              label = sprintf("Parallel LL diff = %.4f < %.4f (15%% of |LL|)", ll_diff_parallel, ll_tol))

  # Parameter estimates should be close (< 15% relative difference)
  # Adaptive quadrature optimizes an approximate likelihood surface, so
  # parameter estimates may differ somewhat from standard quadrature.
  rel_diff_adaptive <- abs(est_adaptive - est_standard) / (abs(est_standard) + 0.1)
  rel_diff_parallel <- abs(est_parallel - est_standard) / (abs(est_standard) + 0.1)

  expect_true(max(rel_diff_adaptive) < 0.15,
              label = sprintf("Max adaptive diff = %.3f < 0.15", max(rel_diff_adaptive)))
  expect_true(max(rel_diff_parallel) < 0.15,
              label = sprintf("Max parallel diff = %.3f < 0.15", max(rel_diff_parallel)))

  # Adaptive single and parallel must give identical results
  # (they compute the exact same likelihood function, just with different parallelization)
  # Allow machine precision tolerance due to floating-point summation order
  expect_equal(est_adaptive, est_parallel, tolerance = .Machine$double.eps^0.5)
})


test_that("Adaptive quadrature validation checks work", {
  skip_on_cran()

  set.seed(54321)
  n <- 100

  f <- rnorm(n)
  dat <- data.frame(
    y1 = f + rnorm(n, 0, 0.5),
    y2 = 0.8 * f + rnorm(n, 0, 0.6),
    intercept = 1
  )

  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                model_type = "linear", loading_normalization = NA_real_)
  ms <- define_model_system(list(mc1, mc2), fm)

  control <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = control,
                                parallel = FALSE, verbose = FALSE)

  fscores <- estimate_factorscores_rcpp(result, dat, control = control, verbose = FALSE)
  factor_scores_mat <- as.matrix(fscores[, "factor_1", drop = FALSE])
  factor_ses_mat <- as.matrix(fscores[, "se_factor_1", drop = FALSE])
  factor_var <- as.numeric(result$estimates["factor_var_1"])

  # Test with factor_scores but NULL factor_ses (should use standard quadrature)
  result2 <- estimate_model_rcpp(
    ms, dat, control = control,
    parallel = FALSE, verbose = FALSE,
    factor_scores = factor_scores_mat,
    factor_ses = NULL,
    factor_vars = c(factor_var_1 = factor_var)
  )

  expect_equal(result2$convergence, 0)
  expect_equal(result2$loglik, result$loglik, tolerance = 1e-6)
})
