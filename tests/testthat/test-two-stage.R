# Two-Stage Estimation Test Suite
#
# This test suite verifies that multi-stage estimation works correctly,
# where early stage components and parameters are fixed while later stages
# are optimized.

test_that("Two-stage estimation: Roy model measurement then selection/wages", {
  skip_on_cran()

  set.seed(42)
  n <- 500  # Smaller sample for faster testing

  # Generate Roy model data
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  f <- rnorm(n)  # Latent ability

  # Test scores (measure ability)
  T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)
  T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)
  T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)

  # Potential wages
  wage0 <- 2.0 + 0.5*x1 + 0.3*x2 + 0.5*f + rnorm(n, 0, 0.6)
  wage1 <- 2.5 + 0.6*x1 + 1.0*f + rnorm(n, 0, 0.7)

  # Sector choice
  z_sector <- 0.0 + 0.4*x2 + 0.8*f
  sector <- as.numeric(runif(n) < pnorm(z_sector))

  # Observed wage (Roy selection)
  wage <- ifelse(sector == 1, wage1, wage0)

  dat <- data.frame(
    intercept = 1,
    x1 = x1, x2 = x2,
    T1 = T1, T2 = T2, T3 = T3,
    wage = wage,
    sector = sector,
    eval_tests = 1,
    eval_wage0 = 1 - sector,
    eval_wage1 = sector,
    eval_sector = 1
  )

  # ======================================================================
  # STAGE 1: Estimate measurement system (test scores only)
  # ======================================================================

  fm <- define_factor_model(n_factors = 1, n_types = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = 1.0,  # Fix for identification
    evaluation_indicator = "eval_tests"
  )

  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Free
    evaluation_indicator = "eval_tests"
  )

  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = NA_real_,  # Free
    evaluation_indicator = "eval_tests"
  )

  ms_stage1 <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3),
    factor = fm
  )

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)
  result_stage1 <- estimate_model_rcpp(
    model_system = ms_stage1,
    data = dat,
    init_params = NULL,
    control = ctrl,
    optimizer = "nlminb",
    verbose = FALSE
  )

  expect_equal(result_stage1$convergence, 0, info = "Stage 1 should converge")
  expect_true(!is.null(result_stage1$estimates), info = "Stage 1 should have estimates")

  n_stage1_params <- length(result_stage1$estimates)
  stage1_loglik <- result_stage1$loglik

  cat(sprintf("\nStage 1: %d parameters, loglik = %.4f\n", n_stage1_params, stage1_loglik))

  # ======================================================================
  # STAGE 2: Add wage/sector equations, fixing stage 1 parameters
  # ======================================================================

  # Reuse same factor model
  mc_wage0 <- define_model_component(
    name = "wage0", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1", "x2"), model_type = "linear",
    loading_normalization = NA_real_,  # Free
    evaluation_indicator = "eval_wage0"
  )

  mc_wage1 <- define_model_component(
    name = "wage1", data = dat, outcome = "wage", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_,  # Free
    evaluation_indicator = "eval_wage1"
  )

  mc_sector <- define_model_component(
    name = "sector", data = dat, outcome = "sector", factor = fm,
    covariates = c("intercept", "x2"), model_type = "probit",
    loading_normalization = NA_real_,  # Free
    evaluation_indicator = "eval_sector"
  )

  ms_stage2 <- define_model_system(
    components = list(mc_wage0, mc_wage1, mc_sector),
    factor = fm,
    previous_stage = result_stage1  # Fix stage 1 parameters
  )

  # Verify previous_stage_info was set up correctly
  expect_false(is.null(ms_stage2$previous_stage_info), info = "Stage 2 should have previous_stage_info")
  expect_equal(ms_stage2$previous_stage_info$n_components, 3, info = "Should have 3 fixed components")
  expect_equal(ms_stage2$previous_stage_info$n_params_fixed, n_stage1_params,
               info = "Should have correct number of fixed parameters")

  # Estimate stage 2
  result_stage2 <- estimate_model_rcpp(
    model_system = ms_stage2,
    data = dat,
    init_params = NULL,
    control = ctrl,
    optimizer = "nlminb",
    verbose = FALSE
  )

  expect_equal(result_stage2$convergence, 0, info = "Stage 2 should converge")

  n_stage2_total <- length(result_stage2$estimates)
  n_stage2_free <- n_stage2_total - n_stage1_params

  cat(sprintf("Stage 2: %d total parameters (%d fixed + %d free), loglik = %.4f\n",
              n_stage2_total, n_stage1_params, n_stage2_free, result_stage2$loglik))

  # ======================================================================
  # VERIFICATION 1: Stage 1 parameters unchanged in stage 2
  # ======================================================================

  stage1_params <- result_stage1$estimates
  stage2_params_fixed <- result_stage2$estimates[1:n_stage1_params]

  max_diff <- max(abs(stage1_params - stage2_params_fixed))
  expect_true(max_diff < 1e-10,
              info = sprintf("Stage 1 params should be unchanged (max diff = %.2e)", max_diff))

  cat(sprintf("  Stage 1 parameters preserved: max diff = %.2e\n", max_diff))

  # ======================================================================
  # VERIFICATION 2: Standard errors preserved from stage 1
  # ======================================================================

  stage1_se <- result_stage1$std_errors
  stage2_se_fixed <- result_stage2$std_errors[1:n_stage1_params]

  max_se_diff <- max(abs(stage1_se - stage2_se_fixed))
  expect_true(max_se_diff < 1e-10,
              info = sprintf("Stage 1 SEs should be preserved (max diff = %.2e)", max_se_diff))

  cat(sprintf("  Stage 1 SEs preserved: max diff = %.2e\n", max_se_diff))

  # ======================================================================
  # VERIFICATION 3: Compare to single-stage (full model) estimation
  # ======================================================================

  # Estimate full model in one stage
  ms_full <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_wage0, mc_wage1, mc_sector),
    factor = fm
  )

  result_full <- estimate_model_rcpp(
    model_system = ms_full,
    data = dat,
    init_params = NULL,
    control = ctrl,
    optimizer = "nlminb",
    verbose = FALSE
  )

  expect_equal(result_full$convergence, 0, info = "Full model should converge")

  cat(sprintf("Full model (single-stage): %d parameters, loglik = %.4f\n",
              length(result_full$estimates), result_full$loglik))

  # Log-likelihoods should be similar (though may differ due to optimization paths)
  # With n=500, a difference of ~1 in loglik is acceptable (0.2% difference)
  loglik_diff <- abs(result_stage2$loglik - result_full$loglik)
  expect_true(loglik_diff < 2.0,
              info = sprintf("Two-stage vs single-stage loglik diff = %.4f", loglik_diff))

  cat(sprintf("  Loglik difference (two-stage vs full): %.4f\n", loglik_diff))

  # Parameters might differ slightly due to different optimization paths,
  # but should be close
  param_diff <- abs(result_stage2$estimates - result_full$estimates)
  max_param_diff <- max(param_diff)

  cat(sprintf("  Max parameter difference: %.4e\n", max_param_diff))
  cat(sprintf("  Mean parameter difference: %.4e\n", mean(param_diff)))

  # This is a soft check - optimizers may find slightly different local optima
  expect_true(max_param_diff < 0.5,
              info = sprintf("Parameters should be similar (max diff = %.4e)", max_param_diff))

  cat("\n✓ All two-stage estimation tests passed!\n")
})
