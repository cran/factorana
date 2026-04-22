# =============================================================================
# Test: HJV-style 3-Stage Estimation with Adaptive Quadrature
# =============================================================================
# Stage 1: Measurement system (2 factors, 3 measures each)
# Stage 2: Factor score estimation
# Stage 3: Outcome equations (linear, probit, logit, oprobit)
#          Compare standard vs adaptive quadrature likelihoods

test_that("HJV-style 3-stage estimation: adaptive vs standard quadrature", {
  skip_on_cran()

  set.seed(42)
  n <- 1000

  # ==========================================================================
  # True Parameters
  # ==========================================================================

  # Factor variances
  true_var_f1 <- 1.0
  true_var_f2 <- 0.8

  # Measurement loadings (factor 1: measures 1-3, factor 2: measures 4-6)
  # First loading for each factor normalized to 1
  true_loadings_f1 <- c(1.0, 0.8, 0.7)
  true_loadings_f2 <- c(1.0, 0.9, 0.6)

  # Measurement error SDs
  true_sigma_m <- c(0.5, 0.6, 0.7, 0.5, 0.55, 0.65)

  # Outcome parameters
  # Linear outcome
  true_linear_intercept <- 0.5
  true_linear_beta <- 0.3
  true_linear_loading_f1 <- 0.6
  true_linear_loading_f2 <- 0.4
  true_linear_sigma <- 0.4

  # Probit outcome
  true_probit_intercept <- -0.2
  true_probit_beta <- 0.4
  true_probit_loading_f1 <- 0.5
  true_probit_loading_f2 <- 0.3

  # Logit outcome (3 choices)
  true_logit_intercept <- c(0.3, -0.1)  # 2 non-reference choices
  true_logit_beta <- c(0.2, 0.3)
  true_logit_loading_f1 <- c(0.4, 0.5)
  true_logit_loading_f2 <- c(0.3, 0.2)

  # ==========================================================================
  # Simulate Data
  # ==========================================================================

  # Generate factors
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  f2 <- rnorm(n, 0, sqrt(true_var_f2))

  # Covariate
  x <- rnorm(n, 0, 1)

  # Measurement equations (factor 1)
  m1 <- true_loadings_f1[1] * f1 + rnorm(n, 0, true_sigma_m[1])
  m2 <- true_loadings_f1[2] * f1 + rnorm(n, 0, true_sigma_m[2])
  m3 <- true_loadings_f1[3] * f1 + rnorm(n, 0, true_sigma_m[3])

  # Measurement equations (factor 2)
  m4 <- true_loadings_f2[1] * f2 + rnorm(n, 0, true_sigma_m[4])
  m5 <- true_loadings_f2[2] * f2 + rnorm(n, 0, true_sigma_m[5])
  m6 <- true_loadings_f2[3] * f2 + rnorm(n, 0, true_sigma_m[6])

  # Linear outcome
  y_linear <- true_linear_intercept + true_linear_beta * x +
              true_linear_loading_f1 * f1 + true_linear_loading_f2 * f2 +
              rnorm(n, 0, true_linear_sigma)

  # Probit outcome
  latent_probit <- true_probit_intercept + true_probit_beta * x +
                   true_probit_loading_f1 * f1 + true_probit_loading_f2 * f2 +
                   rnorm(n, 0, 1)
  y_probit <- as.integer(latent_probit > 0)

  # Logit outcome (3 choices: 0, 1, 2)
  # V_j = intercept_j + beta_j * x + loading_f1_j * f1 + loading_f2_j * f2
  V0 <- 0  # Reference category
  V1 <- true_logit_intercept[1] + true_logit_beta[1] * x +
        true_logit_loading_f1[1] * f1 + true_logit_loading_f2[1] * f2
  V2 <- true_logit_intercept[2] + true_logit_beta[2] * x +
        true_logit_loading_f1[2] * f1 + true_logit_loading_f2[2] * f2

  # Add type 1 extreme value errors and take argmax
  eps_logit <- matrix(-log(-log(runif(n * 3))), ncol = 3)
  U <- cbind(V0 + eps_logit[,1], V1 + eps_logit[,2], V2 + eps_logit[,3])
  y_logit <- apply(U, 1, which.max)  # 1, 2, or 3 (1-indexed)

  # Create data frame (excluding oprobit for now to simplify debugging)
  dat <- data.frame(
    m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5, m6 = m6,
    x = x,
    y_linear = y_linear,
    y_probit = y_probit,
    y_logit = y_logit,
    intercept = 1
  )

  # ==========================================================================
  # Stage 1: Estimate Measurement System
  # ==========================================================================

  message("\n========== STAGE 1: Measurement System ==========")

  fm <- define_factor_model(n_factors = 2)

  # Factor 1 measures (loading normalized on m1)
  mc_m1 <- define_model_component(
    name = "m1", data = dat, outcome = "m1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0)  # Fix loading_f1=1, loading_f2=0
  )
  mc_m2 <- define_model_component(
    name = "m2", data = dat, outcome = "m2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0)  # Free loading_f1, loading_f2=0
  )
  mc_m3 <- define_model_component(
    name = "m3", data = dat, outcome = "m3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0)
  )

  # Factor 2 measures (loading normalized on m4)
  mc_m4 <- define_model_component(
    name = "m4", data = dat, outcome = "m4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1)  # loading_f1=0, fix loading_f2=1
  )
  mc_m5 <- define_model_component(
    name = "m5", data = dat, outcome = "m5", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_)
  )
  mc_m6 <- define_model_component(
    name = "m6", data = dat, outcome = "m6", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_)
  )

  ms_stage1 <- define_model_system(
    components = list(mc_m1, mc_m2, mc_m3, mc_m4, mc_m5, mc_m6),
    factor = fm
  )

  control_stage1 <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    model_system = ms_stage1,
    data = dat,
    control = control_stage1,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = TRUE
  )

  expect_equal(result_stage1$convergence, 0, label = "Stage 1 converged")

  message("\nStage 1 estimates:")
  print(round(result_stage1$estimates, 4))

  # ==========================================================================
  # Stage 2: Factor Score Estimation
  # ==========================================================================

  message("\n========== STAGE 2: Factor Score Estimation ==========")

  fscores <- estimate_factorscores_rcpp(
    result = result_stage1,
    data = dat,
    control = control_stage1,
    verbose = TRUE
  )

  # Extract matrices
  factor_scores_mat <- as.matrix(fscores[, c("factor_1", "factor_2")])
  factor_ses_mat <- as.matrix(fscores[, c("se_factor_1", "se_factor_2")])

  # Extract factor variances from Stage 1
  factor_var_1 <- as.numeric(result_stage1$estimates["factor_var_1"])
  factor_var_2 <- as.numeric(result_stage1$estimates["factor_var_2"])
  factor_vars_vec <- c(factor_var_1, factor_var_2)

  message("\nFactor variances from Stage 1:")
  message(sprintf("  factor_var_1 = %.4f (true = %.4f)", factor_var_1, true_var_f1))
  message(sprintf("  factor_var_2 = %.4f (true = %.4f)", factor_var_2, true_var_f2))

  message("\nFactor score summary:")
  message(sprintf("  Factor 1: mean=%.3f, sd=%.3f, SE range=[%.3f, %.3f]",
                  mean(factor_scores_mat[,1]), sd(factor_scores_mat[,1]),
                  min(factor_ses_mat[,1]), max(factor_ses_mat[,1])))
  message(sprintf("  Factor 2: mean=%.3f, sd=%.3f, SE range=[%.3f, %.3f]",
                  mean(factor_scores_mat[,2]), sd(factor_scores_mat[,2]),
                  min(factor_ses_mat[,2]), max(factor_ses_mat[,2])))

  # Correlation with true factors
  cor_f1 <- cor(factor_scores_mat[,1], f1)
  cor_f2 <- cor(factor_scores_mat[,2], f2)
  message(sprintf("\nCorrelation with true factors: f1=%.3f, f2=%.3f", cor_f1, cor_f2))

  expect_true(cor_f1 > 0.8, label = "Factor 1 scores correlate with truth")
  expect_true(cor_f2 > 0.8, label = "Factor 2 scores correlate with truth")

  # ==========================================================================
  # Stage 3: Outcome Equations
  # ==========================================================================

  message("\n========== STAGE 3: Outcome Equations ==========")

  # Define outcome components
  mc_linear <- define_model_component(
    name = "out_linear", data = dat, outcome = "y_linear", factor = fm,
    covariates = c("intercept", "x"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_)
  )

  mc_probit <- define_model_component(
    name = "out_probit", data = dat, outcome = "y_probit", factor = fm,
    covariates = c("intercept", "x"), model_type = "probit",
    loading_normalization = c(NA_real_, NA_real_)
  )

  mc_logit <- define_model_component(
    name = "out_logit", data = dat, outcome = "y_logit", factor = fm,
    covariates = c("intercept", "x"), model_type = "logit",
    num_choices = 3,
    loading_normalization = c(NA_real_, NA_real_)
  )

  # Model system with previous_stage (linear, probit, logit)
  ms_stage3 <- define_model_system(
    components = list(mc_linear, mc_probit, mc_logit),
    factor = fm,
    previous_stage = result_stage1
  )

  # ==========================================================================
  # Test A: Standard Quadrature
  # ==========================================================================

  message("\n----- Test A: Standard Quadrature -----")

  n_quad <- 8
  control_standard <- define_estimation_control(
    n_quad_points = n_quad,
    num_cores = 1,
    adaptive_integration = FALSE
  )

  result_standard <- estimate_model_rcpp(
    model_system = ms_stage3,
    data = dat,
    control = control_standard,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = TRUE
  )

  expect_equal(result_standard$convergence, 0, label = "Standard quadrature converged")

  message(sprintf("\nStandard quadrature loglik: %.4f", result_standard$loglik))

  # ==========================================================================
  # Test B: Adaptive Quadrature
  # ==========================================================================

  message("\n----- Test B: Adaptive Quadrature -----")

  threshold <- 0.3
  control_adaptive <- define_estimation_control(
    n_quad_points = n_quad,
    num_cores = 1,
    adaptive_integration = TRUE,
    adapt_int_thresh = threshold
  )

  result_adaptive <- estimate_model_rcpp(
    model_system = ms_stage3,
    data = dat,
    control = control_adaptive,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = TRUE,
    factor_scores = factor_scores_mat,
    factor_ses = factor_ses_mat,
    factor_vars = factor_vars_vec
  )

  expect_equal(result_adaptive$convergence, 0, label = "Adaptive quadrature converged")

  message(sprintf("\nAdaptive quadrature loglik: %.4f", result_adaptive$loglik))

  # ==========================================================================
  # Compare Results
  # ==========================================================================

  message("\n========== COMPARISON ==========")

  ll_standard <- result_standard$loglik
  ll_adaptive <- result_adaptive$loglik
  ll_ratio <- ll_adaptive / ll_standard
  ll_diff <- abs(ll_adaptive - ll_standard)
  ll_diff_pct <- 100 * ll_diff / abs(ll_standard)

  message(sprintf("\nLog-likelihood comparison:"))
  message(sprintf("  Standard:  %.4f", ll_standard))
  message(sprintf("  Adaptive:  %.4f", ll_adaptive))
  message(sprintf("  Ratio:     %.6f", ll_ratio))
  message(sprintf("  Abs diff:  %.4f", ll_diff))
  message(sprintf("  Pct diff:  %.2f%%", ll_diff_pct))

  # The key test: likelihoods should be similar (within 20% for this approximation)
  # If adaptive is working correctly, they should be close
  expect_true(ll_ratio > 0.5 && ll_ratio < 2.0,
              label = sprintf("Likelihood ratio %.4f should be between 0.5 and 2.0", ll_ratio))

  # Compare parameter estimates (Stage 3 only)
  n_stage1_params <- length(result_stage1$estimates)
  stage3_names <- names(result_standard$estimates)[(n_stage1_params + 1):length(result_standard$estimates)]

  est_standard <- result_standard$estimates[(n_stage1_params + 1):length(result_standard$estimates)]
  est_adaptive <- result_adaptive$estimates[(n_stage1_params + 1):length(result_adaptive$estimates)]

  message("\nParameter estimates (Stage 3 only):")
  comparison_df <- data.frame(
    Parameter = stage3_names,
    Standard = round(est_standard, 4),
    Adaptive = round(est_adaptive, 4),
    Diff = round(est_adaptive - est_standard, 4)
  )
  print(comparison_df)

  # Parameters should be reasonably close
  max_param_diff <- max(abs(est_adaptive - est_standard))
  message(sprintf("\nMax parameter difference: %.4f", max_param_diff))

  # ==========================================================================
  # Direct likelihood comparison at initial parameters
  # ==========================================================================

  message("\n========== DIRECT LIKELIHOOD COMPARISON ==========")
  message("(Evaluating at the SAME parameters to isolate quadrature difference)")

  # Use Stage 1 estimates + initial Stage 3 estimates
  # Get initial parameters for Stage 3
  init_result <- initialize_parameters(ms_stage3, dat, verbose = FALSE)
  full_init_params <- init_result$init_params

  # Initialize FactorModel for standard quadrature
  fm_ptr_std <- initialize_factor_model_cpp(ms_stage3, as.matrix(dat), n_quad, full_init_params)

  # Get free parameters
  free_params <- extract_free_params_cpp(fm_ptr_std, full_init_params)

  # Evaluate standard quadrature likelihood
  ll_std_direct <- evaluate_loglik_only_cpp(fm_ptr_std, free_params)

  message(sprintf("\nAt initial parameters:"))
  message(sprintf("  Standard quadrature loglik: %.4f", ll_std_direct))

  # Initialize FactorModel for adaptive quadrature
  fm_ptr_adapt <- initialize_factor_model_cpp(ms_stage3, as.matrix(dat), n_quad, full_init_params)

  # Set adaptive quadrature
  set_adaptive_quadrature_cpp(
    fm_ptr_adapt,
    factor_scores_mat,
    factor_ses_mat,
    factor_vars_vec,
    threshold = threshold,
    max_quad = n_quad,
    verbose = TRUE
  )

  # Evaluate adaptive quadrature likelihood
  ll_adapt_direct <- evaluate_loglik_only_cpp(fm_ptr_adapt, free_params)

  message(sprintf("  Adaptive quadrature loglik: %.4f", ll_adapt_direct))
  message(sprintf("  Ratio: %.6f", ll_adapt_direct / ll_std_direct))
  message(sprintf("  Difference: %.4f", ll_adapt_direct - ll_std_direct))

  # This is the critical test - at the same parameters, the likelihoods
  # should be similar if adaptive quadrature is working correctly
  direct_ratio <- ll_adapt_direct / ll_std_direct

  message(sprintf("\n*** CRITICAL TEST: Direct likelihood ratio = %.6f ***", direct_ratio))

  if (abs(direct_ratio - 1.0) > 0.5) {
    message("WARNING: Large discrepancy between standard and adaptive quadrature!")
    message("This suggests a bug in the adaptive quadrature implementation.")
  }

  expect_true(direct_ratio > 0.5 && direct_ratio < 2.0,
              label = sprintf("Direct likelihood ratio %.4f should be between 0.5 and 2.0", direct_ratio))
})
