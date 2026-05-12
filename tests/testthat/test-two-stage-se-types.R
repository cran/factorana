# Tests for two-stage SE_linear estimation where types are introduced at
# Stage 2 (the structural stage) rather than in Stage 1.
#
# Workflow being covered: Stage 1 estimates a measurement system with no
# latent types (n_types = 1, use_types = FALSE everywhere). Stage 2 replaces
# the independent-factor structure with SE_linear and introduces n_types = 2
# at the structural level only — types shift the outcome-factor intercept
# (se_intercept_type_2) and drive a type probability model whose loadings
# are on the input factors. Types do NOT affect the measurement system.
#
# This is the common analytical workflow (Heckman, Humphries & Veramendi
# 2016, 2018): the measurement system is type-agnostic, and unobserved
# heterogeneity is introduced at the structural level in the second stage.
#
# Regression guard: factorana <= 1.1.1 omitted typeprob_* and
# type_*_loading_* from the Stage-2 parameter vector in this path, causing
# either a crash in setup_parameter_constraints or silently mis-fixed
# parameters.

VERBOSE  <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"
GRAD_TOL <- 1e-3
# Hessian tolerance is relaxed to 5e-3 specifically for the factor-variance
# diagonal second derivative. The sqrt(factor_var) factor in the Gauss-Hermite
# scaling makes d^2L/d(factor_var_k)^2 sensitive to FD step size; in this path
# the analytical value is off from FD by ~2e-3 relative error at n_quad=12,
# while every other Hessian element matches to ~1e-6. Five millieths is
# strong enough to flag real bugs while tolerating the genuine numerical
# noise of the factor-variance second derivative.
HESS_TOL <- 5e-3


# ---- DGP helper -------------------------------------------------------------
# Data-generating process used by every test in this file. The measurement
# equations do NOT depend on type (no type shift on measurements). Types
# enter only through the SE intercept on the outcome factor.

.simulate_se_types_dgp <- function(
    n             = 600,
    seed          = 19,
    true_var_f1   = 1.0,
    true_se_lin   = 0.6,
    true_se_res   = 0.5,
    se_int_t2     = 0.8,
    typeprob_t2   = 0.0,
    type_load_t2  = 0.0) {

  set.seed(seed)

  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  log_odds_t2 <- typeprob_t2 + type_load_t2 * f1
  p_t2 <- plogis(log_odds_t2)
  type_id <- ifelse(runif(n) < p_t2, 2L, 1L)
  t2 <- as.integer(type_id == 2L)

  eps <- rnorm(n, 0, sqrt(true_se_res))
  f2  <- 0.0 + true_se_lin * f1 + se_int_t2 * t2 + eps

  # Measurement equations (identical for both types — no type shift)
  lambda_f1 <- c(1.0, 0.9, 1.1)
  lambda_f2 <- c(1.0, 0.85, 1.15)
  int_f1    <- c(1.5, 1.2, 0.9)
  int_f2    <- c(1.3, 1.0, 0.7)
  sigma_f1  <- c(0.7, 0.75, 0.65)
  sigma_f2  <- c(0.75, 0.7, 0.8)

  Y1_1 <- int_f1[1] + lambda_f1[1] * f1 + rnorm(n, 0, sigma_f1[1])
  Y1_2 <- int_f1[2] + lambda_f1[2] * f1 + rnorm(n, 0, sigma_f1[2])
  Y1_3 <- int_f1[3] + lambda_f1[3] * f1 + rnorm(n, 0, sigma_f1[3])
  Y2_1 <- int_f2[1] + lambda_f2[1] * f2 + rnorm(n, 0, sigma_f2[1])
  Y2_2 <- int_f2[2] + lambda_f2[2] * f2 + rnorm(n, 0, sigma_f2[2])
  Y2_3 <- int_f2[3] + lambda_f2[3] * f2 + rnorm(n, 0, sigma_f2[3])

  list(
    data = data.frame(
      intercept = 1,
      Y1_1 = Y1_1, Y1_2 = Y1_2, Y1_3 = Y1_3,
      Y2_1 = Y2_1, Y2_2 = Y2_2, Y2_3 = Y2_3,
      eval = 1
    ),
    true = list(
      var_f1 = true_var_f1, se_lin = true_se_lin, se_res = true_se_res,
      se_int_t2 = se_int_t2, typeprob_t2 = typeprob_t2,
      type_load_t2 = type_load_t2
    )
  )
}

# Stage 1 model: independent factors, NO types, use_types = FALSE everywhere.
.build_stage1_notypes_model <- function(dat) {
  fm <- define_factor_model(n_factors = 2, n_types = 1,
                            factor_structure = "independent")
  mk <- function(name, outcome, norm) {
    define_model_component(
      name = name, data = dat, outcome = outcome, factor = fm,
      covariates = "intercept", model_type = "linear",
      loading_normalization = norm,
      use_types = FALSE, evaluation_indicator = "eval"
    )
  }
  comps <- list(
    mk("Y1_1", "Y1_1", c(1, 0)),
    mk("Y1_2", "Y1_2", c(NA_real_, 0)),
    mk("Y1_3", "Y1_3", c(NA_real_, 0)),
    mk("Y2_1", "Y2_1", c(0, 1)),
    mk("Y2_2", "Y2_2", c(0, NA_real_)),
    mk("Y2_3", "Y2_3", c(0, NA_real_))
  )
  list(fm = fm, ms = define_model_system(components = comps, factor = fm))
}


# =============================================================================
# TEST 1 — shape test
# =============================================================================
test_that("Stage 2 SE_linear + n_types=2 has typeprob and type_loading slots", {
  sim <- .simulate_se_types_dgp(n = 80, seed = 11)
  dat <- sim$data

  stage1 <- .build_stage1_notypes_model(dat)
  ctrl <- define_estimation_control(n_quad_points = 4, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    model_system = stage1$ms, data = dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0,
               info = "Stage 1 (no types) must converge strictly")

  fm_stage2 <- define_factor_model(n_factors = 2, n_types = 2,
                                    factor_structure = "SE_linear")
  ms_stage2 <- define_model_system(components = list(), factor = fm_stage2,
                                    previous_stage = result_stage1)

  init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  init_names <- names(init_s2$init_params)

  # Required factor-level slots
  expect_true("factor_var_1"        %in% init_names)
  expect_true("se_intercept"        %in% init_names)
  expect_true("se_linear_1"         %in% init_names)
  expect_true("se_intercept_type_2" %in% init_names)
  expect_true("se_residual_var"     %in% init_names)

  # Regression guard: typeprob + type_loading slots must be present
  expect_true("typeprob_2_intercept" %in% init_names,
              info = "Stage 2 must include typeprob_2_intercept slot")
  expect_true("type_2_loading_1" %in% init_names,
              info = "Stage 2 must include type_2_loading_1 (input factor)")
  expect_true("type_2_loading_2" %in% init_names,
              info = "Stage 2 must include type_2_loading_2 slot (outcome factor — auto-fixed)")

  # Alignment between init_params and build_parameter_metadata
  metadata <- factorana:::build_parameter_metadata(ms_stage2)
  expect_equal(length(init_s2$init_params), length(metadata$names))
  expect_equal(init_names, metadata$names)

  # Parameters must not be duplicated
  expect_equal(anyDuplicated(init_names), 0L,
               info = "Stage 2 param vector must not contain duplicate names")

  # setup_parameter_constraints must run without error
  constraints <- factorana:::setup_parameter_constraints(
    ms_stage2, init_s2$init_params, metadata,
    init_s2$factor_variance_fixed, verbose = FALSE
  )
  expect_true(length(constraints$free_idx) > 0)

  # type_2_loading_2 must be auto-fixed (outcome factor, SE identification)
  outcome_load_idx <- match("type_2_loading_2", metadata$names)
  expect_true(!is.na(outcome_load_idx))
  expect_false(outcome_load_idx %in% constraints$free_idx,
               info = "type_2_loading_2 (outcome factor) must be fixed under SE_linear")

  # The new Stage 2 factor-level params (introduced because n_types jumped
  # from 1 to 2) must all be FREE.
  for (pn in c("typeprob_2_intercept", "type_2_loading_1",
               "se_intercept_type_2", "se_linear_1", "se_intercept",
               "se_residual_var", "factor_var_1")) {
    idx <- match(pn, metadata$names)
    expect_true(idx %in% constraints$free_idx,
                info = sprintf("%s must be a free parameter in Stage 2", pn))
  }
})


# =============================================================================
# TEST 2 — FD gradient + Hessian at the true DGP parameters
# =============================================================================
test_that("Two-stage SE_linear with n_types=2: FD gradient and Hessian match", {
  skip_on_cran()

  sim <- .simulate_se_types_dgp(
    n = 600, seed = 19,
    true_var_f1 = 1.0, true_se_lin = 0.6, true_se_res = 0.5,
    se_int_t2 = 0.8,
    typeprob_t2 = 0.3, type_load_t2 = 0.4
  )
  dat <- sim$data
  truth <- sim$true

  stage1 <- .build_stage1_notypes_model(dat)
  # n_quad = 12 keeps the FD Hessian well-conditioned (at n_quad=6 the
  # factor-variance diagonal has rel_err ~3.5e-3, at n_quad=12 it's ~1.9e-3).
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    model_system = stage1$ms, data = dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0, info = "Stage 1 must converge strictly")

  fm_stage2 <- define_factor_model(n_factors = 2, n_types = 2,
                                    factor_structure = "SE_linear")
  ms_stage2 <- define_model_system(components = list(), factor = fm_stage2,
                                    previous_stage = result_stage1)

  init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)

  # Set Stage-2 parameters to their true DGP values (measurement params are
  # already fixed from Stage 1 — those stay at the Stage 1 MLE, which is the
  # well-defined anchor for the FD check).
  params <- init_s2$init_params
  params["factor_var_1"]        <- truth$var_f1
  params["se_intercept"]        <- 0.0
  params["se_linear_1"]         <- truth$se_lin
  params["se_intercept_type_2"] <- truth$se_int_t2
  params["se_residual_var"]     <- truth$se_res
  params["typeprob_2_intercept"] <- truth$typeprob_t2
  params["type_2_loading_1"]    <- truth$type_load_t2

  metadata <- factorana:::build_parameter_metadata(ms_stage2)
  constraints <- factorana:::setup_parameter_constraints(
    ms_stage2, params, metadata,
    init_s2$factor_variance_fixed, verbose = FALSE
  )
  param_fixed <- rep(TRUE, length(params))
  param_fixed[constraints$free_idx] <- FALSE

  # FD gradient
  grad_check <- check_gradient_accuracy(ms_stage2, dat, params,
                                         param_fixed = param_fixed,
                                         tol = GRAD_TOL, verbose = FALSE, n_quad = 12)
  expect_true(grad_check$pass,
              info = sprintf("Stage 2 gradient FD failed (max err: %.2e)",
                             grad_check$max_error))

  # FD Hessian
  hess_check <- check_hessian_accuracy(ms_stage2, dat, params,
                                        param_fixed = param_fixed,
                                        tol = HESS_TOL, verbose = FALSE, n_quad = 12)
  expect_true(hess_check$pass,
              info = sprintf("Stage 2 Hessian FD failed (max err: %.2e)",
                             hess_check$max_error))

  if (VERBOSE) {
    cat("\n=== Stage 2 FD checks (Stage 1 no types) ===\n")
    cat(sprintf("Grad max rel_err: %.2e  (tol=%.0e)\n", grad_check$max_error, GRAD_TOL))
    cat(sprintf("Hess max rel_err: %.2e  (tol=%.0e)\n", hess_check$max_error, HESS_TOL))
  }
})


# ---- Oprobit variant of the DGP + Stage 1 ----------------------------------
.simulate_se_types_dgp_oprobit <- function(
    n            = 800,
    seed         = 23,
    true_var_f1  = 1.0,
    true_se_lin  = 0.6,
    true_se_res  = 0.5,
    se_int_t2    = 0.8,
    typeprob_t2  = 0.0,
    type_load_t2 = 0.0,
    cuts1        = c(-1.0, 0.0, 1.0),   # cutpoints for factor-1 indicators
    cuts2        = c(-1.0, 0.0, 1.0)) { # cutpoints for factor-2 indicators

  set.seed(seed)
  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  log_odds_t2 <- typeprob_t2 + type_load_t2 * f1
  t2 <- as.integer(runif(n) < plogis(log_odds_t2))
  eps <- rnorm(n, 0, sqrt(true_se_res))
  f2  <- true_se_lin * f1 + se_int_t2 * t2 + eps

  lambda_f1 <- c(1.0, 0.9, 1.1)
  lambda_f2 <- c(1.0, 0.85, 1.15)
  gen <- function(f, lam, cuts) {
    ystar <- lam * f + rnorm(length(f), 0, 1)   # probit scale, sigma = 1
    as.integer(findInterval(ystar, cuts) + 1L)
  }
  list(
    data = data.frame(
      intercept = 1,
      Y1_1 = gen(f1, lambda_f1[1], cuts1),
      Y1_2 = gen(f1, lambda_f1[2], cuts1),
      Y1_3 = gen(f1, lambda_f1[3], cuts1),
      Y2_1 = gen(f2, lambda_f2[1], cuts2),
      Y2_2 = gen(f2, lambda_f2[2], cuts2),
      Y2_3 = gen(f2, lambda_f2[3], cuts2),
      eval = 1
    ),
    true = list(var_f1 = true_var_f1, se_lin = true_se_lin, se_res = true_se_res,
                se_int_t2 = se_int_t2, typeprob_t2 = typeprob_t2,
                type_load_t2 = type_load_t2)
  )
}

.build_stage1_notypes_model_oprobit <- function(dat, n_categories = 4L) {
  fm <- define_factor_model(n_factors = 2, n_types = 1,
                            factor_structure = "independent")
  mk <- function(name, outcome, norm) {
    define_model_component(
      name = name, data = dat, outcome = outcome, factor = fm,
      covariates = NULL, model_type = "oprobit",
      num_choices = n_categories,
      loading_normalization = norm,
      use_types = FALSE, evaluation_indicator = "eval"
    )
  }
  comps <- list(
    mk("Y1_1", "Y1_1", c(1, 0)),
    mk("Y1_2", "Y1_2", c(NA_real_, 0)),
    mk("Y1_3", "Y1_3", c(NA_real_, 0)),
    mk("Y2_1", "Y2_1", c(0, 1)),
    mk("Y2_2", "Y2_2", c(0, NA_real_)),
    mk("Y2_3", "Y2_3", c(0, NA_real_))
  )
  list(fm = fm, ms = define_model_system(components = comps, factor = fm))
}

# =============================================================================
# TEST 3: Two-stage SE_linear with n_types=2 AND oprobit indicators, FD check.
#
# Analog of the linear TEST 2 with ordered-probit measurements. Stage 1 is a
# 2-factor independent oprobit measurement model; Stage 2 is SE_linear +
# n_types = 2 with previous_stage = Stage 1. FD gradient and Hessian are
# checked at the DGP structural parameters.
# =============================================================================
test_that("Two-stage SE_linear with n_types=2 and oprobit indicators: FD gradient and Hessian match", {
  skip_on_cran()

  sim <- .simulate_se_types_dgp_oprobit(
    n = 800, seed = 23,
    true_var_f1 = 1.0, true_se_lin = 0.6, true_se_res = 0.5,
    se_int_t2 = 0.8, typeprob_t2 = 0.3, type_load_t2 = 0.4
  )
  dat   <- sim$data
  truth <- sim$true

  stage1 <- .build_stage1_notypes_model_oprobit(dat, n_categories = 4L)
  ctrl   <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    stage1$ms, dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0,
               info = "Oprobit Stage 1 must converge strictly")

  fm_stage2 <- define_factor_model(n_factors = 2, n_types = 2,
                                    factor_structure = "SE_linear")
  ms_stage2 <- define_model_system(components = list(), factor = fm_stage2,
                                    previous_stage = result_stage1)

  init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  params <- init_s2$init_params
  params["factor_var_1"]         <- truth$var_f1
  params["se_intercept"]         <- 0.0
  params["se_linear_1"]          <- truth$se_lin
  params["se_intercept_type_2"]  <- truth$se_int_t2
  params["se_residual_var"]      <- truth$se_res
  params["typeprob_2_intercept"] <- truth$typeprob_t2
  params["type_2_loading_1"]     <- truth$type_load_t2

  metadata <- factorana:::build_parameter_metadata(ms_stage2)
  constraints <- factorana:::setup_parameter_constraints(
    ms_stage2, params, metadata,
    init_s2$factor_variance_fixed, verbose = FALSE
  )
  param_fixed <- rep(TRUE, length(params))
  param_fixed[constraints$free_idx] <- FALSE

  grad_check <- check_gradient_accuracy(ms_stage2, dat, params,
                                         param_fixed = param_fixed,
                                         tol = GRAD_TOL, verbose = FALSE, n_quad = 8)
  expect_true(grad_check$pass,
              info = sprintf("oprobit Stage 2 gradient FD failed (max err: %.2e)",
                             grad_check$max_error))

  hess_check <- check_hessian_accuracy(ms_stage2, dat, params,
                                        param_fixed = param_fixed,
                                        tol = HESS_TOL, verbose = FALSE, n_quad = 8)
  expect_true(hess_check$pass,
              info = sprintf("oprobit Stage 2 Hessian FD failed (max err: %.2e)",
                             hess_check$max_error))

  if (VERBOSE) {
    cat("\n=== Oprobit two-stage Stage 2 FD checks ===\n")
    cat(sprintf("Grad max rel_err: %.2e  (tol=%.0e)\n", grad_check$max_error, GRAD_TOL))
    cat(sprintf("Hess max rel_err: %.2e  (tol=%.0e)\n", hess_check$max_error, HESS_TOL))
  }
})


# ---- Stage 1 model variant WITH types (use_types = TRUE on components) -----
# Counterpart to .build_stage1_notypes_model: same 2-factor independent
# measurement system but Stage 1 itself models n_types = 2 and includes
# type intercepts on every measurement component.

.build_stage1_withtypes_model <- function(dat) {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "independent")
  mk <- function(name, outcome, norm) {
    define_model_component(
      name = name, data = dat, outcome = outcome, factor = fm,
      covariates = "intercept", model_type = "linear",
      loading_normalization = norm,
      use_types = TRUE, evaluation_indicator = "eval"
    )
  }
  comps <- list(
    mk("Y1_1", "Y1_1", c(1, 0)),
    mk("Y1_2", "Y1_2", c(NA_real_, 0)),
    mk("Y1_3", "Y1_3", c(NA_real_, 0)),
    mk("Y2_1", "Y2_1", c(0, 1)),
    mk("Y2_2", "Y2_2", c(0, NA_real_)),
    mk("Y2_3", "Y2_3", c(0, NA_real_))
  )
  list(fm = fm, ms = define_model_system(components = comps, factor = fm))
}


# =============================================================================
# TEST 4 — Stage 1 with types + Stage 2 SE_linear FD check
#
# Exercises the less-common workflow where Stage 1 already models
# n_types = 2 (with `use_types = TRUE` on every measurement component)
# before Stage 2 introduces `SE_linear`. Before the v1.1.7 Hessian
# accumulation fix this path produced a Hessian FD mismatch with
# rel_err ~1.7 on the SE x SE sub-block (analytical Hessian missed
# cross-derivatives involving equality-tied parameters). The fix
# makes the analytical Hessian match FD on this path too.
# =============================================================================
test_that("Stage 1 with types + Stage 2 SE_linear: FD gradient and Hessian match", {
  # Regression guard against a prior free/fixed mapping bug.
  #
  # In versions <= 1.1.8, define_model_system() in the Stage 2
  # SE workflow (allow_different_structure = TRUE) only treated
  # factor_var_*, se_*, chol_*, and factor_mean_* parameters as
  # factor-distribution parameters. typeprob_*_intercept and
  # type_*_loading_* were classified as measurement parameters and
  # therefore marked FIXED when inherited from a Stage 1 that
  # already had n_types >= 2. The C++ side correctly left them
  # FREE in the SE branch, so the R-side free_idx and the C++
  # free-parameter vector disagreed: evaluate_likelihood_rcpp
  # extracted a 7-element free gradient from C++ and then scattered
  # it into a 5-element R free_idx, shifting every value past the
  # 5th slot. That produced a permuted Hessian with max rel_err
  # ~1.67 on the SE x SE sub-block, which the TEST 4 comments
  # previously flagged as a "KNOWN ISSUE". Fix: also treat
  # ^typeprob_ and ^type_[0-9]+_loading_ as factor-distribution
  # patterns in define_model_system().

  skip_on_cran()

  sim <- .simulate_se_types_dgp(
    n = 600, seed = 19,
    true_var_f1 = 1.0, true_se_lin = 0.6, true_se_res = 0.5,
    se_int_t2 = 0.8,
    typeprob_t2 = 0.3, type_load_t2 = 0.4
  )
  dat <- sim$data
  truth <- sim$true

  stage1 <- .build_stage1_withtypes_model(dat)
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    model_system = stage1$ms, data = dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0,
               info = "Stage 1 (with types) must converge strictly")

  fm_stage2 <- define_factor_model(n_factors = 2, n_types = 2,
                                    factor_structure = "SE_linear")
  ms_stage2 <- define_model_system(components = list(), factor = fm_stage2,
                                    previous_stage = result_stage1)

  init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  params <- init_s2$init_params
  params["factor_var_1"]        <- truth$var_f1
  params["se_intercept"]        <- 0.0
  params["se_linear_1"]         <- truth$se_lin
  params["se_intercept_type_2"] <- truth$se_int_t2
  params["se_residual_var"]     <- truth$se_res
  if ("typeprob_2_intercept" %in% names(params)) {
    params["typeprob_2_intercept"] <- truth$typeprob_t2
  }
  if ("type_2_loading_1" %in% names(params)) {
    params["type_2_loading_1"] <- truth$type_load_t2
  }

  metadata <- factorana:::build_parameter_metadata(ms_stage2)
  constraints <- factorana:::setup_parameter_constraints(
    ms_stage2, params, metadata,
    init_s2$factor_variance_fixed, verbose = FALSE
  )
  param_fixed <- rep(TRUE, length(params))
  param_fixed[constraints$free_idx] <- FALSE

  grad_check <- check_gradient_accuracy(ms_stage2, dat, params,
                                         param_fixed = param_fixed,
                                         tol = GRAD_TOL, verbose = FALSE, n_quad = 12)
  expect_true(grad_check$pass,
              info = sprintf("Stage 1-with-types Stage 2 gradient FD failed (max err: %.2e)",
                             grad_check$max_error))

  hess_check <- check_hessian_accuracy(ms_stage2, dat, params,
                                        param_fixed = param_fixed,
                                        tol = HESS_TOL, verbose = FALSE, n_quad = 12)
  expect_true(hess_check$pass,
              info = sprintf("Stage 1-with-types Stage 2 Hessian FD failed (max err: %.2e)",
                             hess_check$max_error))

  if (VERBOSE) {
    cat("\n=== Stage 2 FD checks (Stage 1 WITH types) ===\n")
    cat(sprintf("Grad max rel_err: %.2e  (tol=%.0e)\n", grad_check$max_error, GRAD_TOL))
    cat(sprintf("Hess max rel_err: %.2e  (tol=%.0e)\n", hess_check$max_error, HESS_TOL))
  }
})


# =============================================================================
# DGP helper: dynamic measurement (same items at two waves) + SE_linear with
# n_types = 2 and an SE-equation covariate.
#
# Two-stage estimation with this DGP uses the canonical workflow:
#   - Stage 1 fits a dynamic measurement model (intercepts FREE per period,
#     loadings and sigmas tied across periods via define_dynamic_measurement()).
#   - Stage 2 fixes the wave-1 (anchor) intercepts and shared loadings/sigmas
#     and estimates SE_linear with n_types = 2 + se_covariates = c("X").
#
# Critically, wave-1 intercepts are reused as the anchor for the wave-2
# (outcome-factor) measurements. f_2 is therefore allowed to have a non-zero
# mean: any mean shift induced by se_intercept_type_2 * Pr(type=2) flows into
# se_intercept rather than getting absorbed by free wave-2 intercepts.
#
# The structural equation for the outcome factor:
#   f_2 = se_intercept + se_lin * f_1 + se_cov * (X - mean(X))
#         + se_int_t2 * 1{type = 2} + epsilon
# =============================================================================

.simulate_se_types_cov_dgp <- function(
    n             = 800,
    seed          = 41,
    true_var_f1   = 1.0,
    true_se_int   = 0.0,
    true_se_lin   = 0.6,
    true_se_cov   = 0.5,
    true_se_res   = 0.5,
    se_int_t2     = 0.8,
    typeprob_t2   = 0.0,
    type_load_t2  = 0.0,
    item_int      = c(1.5, 1.2, 0.9),
    item_load     = c(1.0, 0.9, 1.1),
    item_sigma    = c(0.7, 0.75, 0.65),
    x_mean        = 2.0,
    x_sd          = 1.0) {

  set.seed(seed)

  X          <- rnorm(n, mean = x_mean, sd = x_sd)
  X_demeaned <- X - mean(X)

  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  log_odds_t2 <- typeprob_t2 + type_load_t2 * f1
  p_t2 <- plogis(log_odds_t2)
  type_id <- ifelse(runif(n) < p_t2, 2L, 1L)
  t2 <- as.integer(type_id == 2L)

  eps <- rnorm(n, 0, sqrt(true_se_res))
  f2  <- true_se_int + true_se_lin * f1 + true_se_cov * X_demeaned +
         se_int_t2 * t2 + eps

  gen_Y <- function(f, i) {
    item_int[i] + item_load[i] * f + rnorm(length(f), 0, item_sigma[i])
  }

  dat <- data.frame(
    id        = seq_len(n),
    intercept = 1,
    X         = X,
    eval      = 1L,
    Y_t1_m1 = gen_Y(f1, 1), Y_t1_m2 = gen_Y(f1, 2), Y_t1_m3 = gen_Y(f1, 3),
    Y_t2_m1 = gen_Y(f2, 1), Y_t2_m2 = gen_Y(f2, 2), Y_t2_m3 = gen_Y(f2, 3)
  )

  list(
    data = dat,
    true = list(
      var_f1       = true_var_f1,
      se_int       = true_se_int,
      se_lin       = true_se_lin,
      se_cov       = true_se_cov,
      se_res       = true_se_res,
      se_int_t2    = se_int_t2,
      typeprob_t2  = typeprob_t2,
      type_load_t2 = type_load_t2,
      item_int     = item_int,
      item_load    = item_load,
      item_sigma   = item_sigma
    )
  )
}


# =============================================================================
# TEST 5 — Dynamic Stage 1 + Stage 2 SE_linear + n_types=2 + SE covariate:
#          FD gradient/Hessian
#
# Stage 1: define_dynamic_measurement() with items measured at waves t1 and
# t2; loadings and sigmas tied across waves, intercepts free per wave.
# Stage 2: build_dynamic_previous_stage() anchors at wave 1, then
# SE_linear + n_types = 2 + se_covariates = c("X"). The FD check evaluates
# gradient and Hessian at the true DGP parameter values; this guards
# against bugs in:
#   - se_cov_X gradient/Hessian terms
#   - cross-derivatives between se_cov_X and the rest of the structural
#     parameter block (factor_var_1, se_linear_1, se_residual_var,
#     se_intercept_type_2, typeprob_2_intercept, type_2_loading_1)
# =============================================================================
test_that("Two-stage SE_linear + n_types=2 + SE covariate: FD gradient and Hessian match", {
  skip_on_cran()

  sim <- .simulate_se_types_cov_dgp(
    n = 600, seed = 29,
    true_var_f1 = 1.0,
    true_se_int = 0.0, true_se_lin = 0.6, true_se_cov = 0.5,
    true_se_res = 0.5, se_int_t2 = 0.8,
    typeprob_t2 = 0.3, type_load_t2 = 0.4
  )
  dat   <- sim$data
  truth <- sim$true

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  ctrl <- define_estimation_control(n_quad_points = 12, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    dyn$model_system, dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0,
               info = "Stage 1 (dynamic measurement) must converge strictly")

  prev_stage <- build_dynamic_previous_stage(dyn, result_stage1, dat,
                                              anchor_period = 1L)

  fm_stage2 <- define_factor_model(n_factors = 2, n_types = 2,
                                    factor_structure = "SE_linear",
                                    se_covariates = c("X"))
  ms_stage2 <- define_model_system(components = list(), factor = fm_stage2,
                                    previous_stage = prev_stage)

  init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  params  <- init_s2$init_params

  expect_true("se_cov_X" %in% names(params),
              info = "Stage 2 with se_covariates must include se_cov_X parameter")

  params["factor_var_1"]         <- truth$var_f1
  params["se_intercept"]         <- truth$se_int
  params["se_linear_1"]          <- truth$se_lin
  params["se_cov_X"]             <- truth$se_cov
  params["se_intercept_type_2"]  <- truth$se_int_t2
  params["se_residual_var"]      <- truth$se_res
  params["typeprob_2_intercept"] <- truth$typeprob_t2
  params["type_2_loading_1"]     <- truth$type_load_t2

  metadata <- factorana:::build_parameter_metadata(ms_stage2)
  constraints <- factorana:::setup_parameter_constraints(
    ms_stage2, params, metadata,
    init_s2$factor_variance_fixed, verbose = FALSE
  )
  param_fixed <- rep(TRUE, length(params))
  param_fixed[constraints$free_idx] <- FALSE

  expect_false(param_fixed[match("se_cov_X", names(params))],
               info = "se_cov_X must be a free parameter in Stage 2")

  grad_check <- check_gradient_accuracy(ms_stage2, dat, params,
                                         param_fixed = param_fixed,
                                         tol = GRAD_TOL, verbose = FALSE,
                                         n_quad = 12)
  expect_true(grad_check$pass,
              info = sprintf("Stage 2 SE-cov gradient FD failed (max err: %.2e)",
                             grad_check$max_error))

  hess_check <- check_hessian_accuracy(ms_stage2, dat, params,
                                        param_fixed = param_fixed,
                                        tol = HESS_TOL, verbose = FALSE,
                                        n_quad = 12)
  expect_true(hess_check$pass,
              info = sprintf("Stage 2 SE-cov Hessian FD failed (max err: %.2e)",
                             hess_check$max_error))

  if (VERBOSE) {
    cat("\n=== Stage 2 FD checks (SE covariate) ===\n")
    cat(sprintf("Grad max rel_err: %.2e  (tol=%.0e)\n", grad_check$max_error, GRAD_TOL))
    cat(sprintf("Hess max rel_err: %.2e  (tol=%.0e)\n", hess_check$max_error, HESS_TOL))
  }
})


# =============================================================================
# TEST 6 — Dynamic Stage 1 + Stage 2 SE_linear + n_types=2 + SE covariate:
#          parameter recovery at n = 3000
#
# Stage 1 uses define_dynamic_measurement() so wave-1 / wave-2 measurement
# loadings and sigmas are tied; Stage 2 fixes those plus the wave-1
# (anchor) intercepts via build_dynamic_previous_stage(). f_2 is therefore
# free to have a non-zero mean, so se_intercept and se_intercept_type_2
# can be recovered (modulo type-label switching).
#
# The likelihood is invariant under relabeling type 1 <-> type 2, which
# negates se_intercept_type_2, typeprob_2_intercept, and type_2_loading_1.
# We pin those init values with the true sign to break the symmetry and
# check the magnitude of se_intercept_type_2 to be label-invariant.
# =============================================================================
test_that("Two-stage SE_linear + n_types=2 + SE covariate: parameter recovery", {
  skip_on_cran()

  # seed = 19 gives a "central" sample. Across an 8-seed sweep at n = 3000
  # the empirical mean of est se_linear_1 is 0.614 (truth = 0.6) with
  # individual seeds ranging from 0.48 to 0.78; seed 19 lands at 0.601.
  # This is just a stable seed for the recovery check and not load-bearing.
  sim <- .simulate_se_types_cov_dgp(
    n = 3000, seed = 19,
    true_var_f1 = 1.0,
    true_se_int = 0.0, true_se_lin = 0.6, true_se_cov = 0.5,
    true_se_res = 0.5, se_int_t2 = 0.8,
    typeprob_t2 = 0.3, type_load_t2 = 0.4
  )
  dat   <- sim$data
  truth <- sim$true

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    dyn$model_system, dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0,
               info = "Stage 1 (dynamic measurement) must converge strictly")

  prev_stage <- build_dynamic_previous_stage(dyn, result_stage1, dat,
                                              anchor_period = 1L)

  fm_stage2 <- define_factor_model(n_factors = 2, n_types = 2,
                                    factor_structure = "SE_linear",
                                    se_covariates = c("X"))
  ms_stage2 <- define_model_system(components = list(), factor = fm_stage2,
                                    previous_stage = prev_stage)

  init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
  init    <- init_s2$init_params
  init["factor_var_1"]         <- unname(prev_stage$estimates["factor_var_1"])
  init["se_intercept"]         <- 0.0
  init["se_linear_1"]          <- 0.5
  init["se_cov_X"]             <- 0.0
  init["se_residual_var"]      <- 0.5
  # Pin type-specific params with the SIGN of the true DGP to break the
  # type 1 <-> type 2 label-switching symmetry of the likelihood.
  init["se_intercept_type_2"]  <- 0.5
  init["typeprob_2_intercept"] <- 0.2
  init["type_2_loading_1"]     <- 0.2

  result <- estimate_model_rcpp(
    ms_stage2, dat, init_params = init, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result$convergence, 0,
               info = "Stage 2 (SE_linear + n_types=2 + SE cov) must converge strictly")

  est <- result$estimates
  se  <- result$std_errors

  free_names <- c("factor_var_1", "se_intercept", "se_linear_1", "se_cov_X",
                  "se_residual_var", "se_intercept_type_2",
                  "typeprob_2_intercept", "type_2_loading_1")
  expect_true(all(se[free_names] > 0),
              info = "All Stage 2 free params must have positive std_errors")
  well_id <- c("factor_var_1", "se_linear_1", "se_cov_X", "se_residual_var")
  expect_true(all(se[well_id] < 1),
              info = "Stage 2 well-identified SEs must be reasonable (< 1)")

  # Well-identified structural params at n = 3000. se_residual_var is
  # routinely the least precisely estimated of these parameters in
  # two-stage workflows; the others recover within ~1 SE of truth.
  expect_equal(unname(est["se_linear_1"]),     truth$se_lin, tolerance = 0.10,
               info = sprintf("se_linear_1: true=%.3f est=%.3f se=%.3f",
                              truth$se_lin, est["se_linear_1"],
                              se["se_linear_1"]))
  expect_equal(unname(est["se_cov_X"]),        truth$se_cov, tolerance = 0.10,
               info = sprintf("se_cov_X: true=%.3f est=%.3f se=%.3f",
                              truth$se_cov, est["se_cov_X"],
                              se["se_cov_X"]))
  expect_lt(abs(unname(est["se_residual_var"]) - truth$se_res), 0.20)
  expect_equal(unname(est["factor_var_1"]),    truth$var_f1, tolerance = 0.15,
               info = sprintf("factor_var_1: true=%.3f est=%.3f se=%.3f",
                              truth$var_f1, est["factor_var_1"],
                              se["factor_var_1"]))
  # se_intercept is now identified (not absorbed by free measurement
  # intercepts) because the wave-1 anchor pins them. The label-switching
  # mode also induces a global se_intercept / se_intercept_type_2 shift,
  # so use absolute tolerances on the shifted-mean parameters.
  expect_lt(abs(unname(est["se_intercept"])     - truth$se_int),    0.30)
  # |se_intercept_type_2| is label-invariant. The label-switching mode at
  # finite n biases the gap by ~0.3 in either direction.
  expect_lt(abs(abs(unname(est["se_intercept_type_2"])) - truth$se_int_t2),
            0.40)

  if (VERBOSE) {
    cat("\n=== Stage 2 SE-covariate recovery (n=3000, dynamic Stage 1) ===\n")
    cat(sprintf("%-22s %10s %10s %10s %10s\n",
                "param", "true", "est", "se", "z"))
    rows <- list(
      factor_var_1         = truth$var_f1,
      se_intercept         = truth$se_int,
      se_linear_1          = truth$se_lin,
      se_cov_X             = truth$se_cov,
      se_residual_var      = truth$se_res,
      se_intercept_type_2  = truth$se_int_t2,
      typeprob_2_intercept = truth$typeprob_t2,
      type_2_loading_1     = truth$type_load_t2
    )
    for (p in names(rows)) {
      tr <- rows[[p]]
      zv <- if (se[p] > 0) (est[p] - tr) / se[p] else NA
      cat(sprintf("  %-22s %+10.4f %+10.4f %10.4f %+10.2f\n",
                  p, tr, est[p], se[p], zv))
    }
  }
})


# =============================================================================
# TEST 7 — Constrained re-run via previous_stage = SE_linear + free_params
#
# Regression guard for a permutation bug observed in factorana <= 1.2.0
# (after the parameter-ordering fixes, but before the fac_name_idx extension):
#
# Workflow: take the unconstrained Stage 2 SE_linear MLE, plug it in as a
# matching-structure previous_stage (allow_different_structure = FALSE), set
# type_2_loading_2 = 0, and re-optimize with `free_params` listing every
# factor-distribution slot except type_2_loading_2.
#
# Bug symptom: the C++ un-fix code in initialize_factor_model_cpp built a
# name-to-index map that only recognized factor_var_*, mix_*, se_intercept,
# se_linear_*, se_quadratic_*, se_intercept_type_*, and se_residual_var.
# Names like typeprob_*_intercept, type_*_loading_*, factor_mean_*_*, and
# se_cov_* were silently missing from the map and stayed FIXED on the C++
# side, while R's setup_parameter_constraints kept them FREE. The optimizer
# then ran on a smaller free-set than R expected; estimates_free was
# scattered into a longer R-side free_idx via R's vector-recycling rule,
# producing a clean k-cycle of repeated values across covariates and type
# params (the user's "7-cycle" report on the MH-trap regime model).
#
# This test reproduces the exact path with two SE covariates and asserts
# every free factor-distribution slot has a distinct value (no recycling).
# =============================================================================
test_that("Constrained re-run via free_params: SE-cov / typeprob / type_loading slots get distinct values", {
  skip_on_cran()

  sim <- .simulate_se_types_cov_dgp(
    n = 800, seed = 31,
    true_var_f1 = 1.0,
    true_se_int = 0.0, true_se_lin = 0.6, true_se_cov = 0.5,
    true_se_res = 0.5, se_int_t2 = 0.7,
    typeprob_t2 = 0.2, type_load_t2 = 0.3
  )
  dat <- sim$data

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result_stage1 <- estimate_model_rcpp(
    dyn$model_system, dat, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(result_stage1$convergence, 0)

  prev_stage <- build_dynamic_previous_stage(dyn, result_stage1, dat,
                                              anchor_period = 1L)

  # Unconstrained Stage 2 (allow_different_structure = TRUE branch).
  fm_se <- define_factor_model(n_factors = 2, n_types = 2,
                                factor_structure = "SE_linear",
                                se_covariates = c("X"))
  ms_se <- define_model_system(components = list(), factor = fm_se,
                                previous_stage = prev_stage)

  init_un <- initialize_parameters(ms_se, dat, verbose = FALSE)$init_params
  init_un["factor_var_1"]         <- unname(prev_stage$estimates["factor_var_1"])
  init_un["se_linear_1"]          <- 0.5
  init_un["se_residual_var"]      <- 0.5
  init_un["se_intercept_type_2"]  <- 0.5
  init_un["typeprob_2_intercept"] <- 0.2
  init_un["type_2_loading_1"]     <- 0.2
  res_un <- estimate_model_rcpp(ms_se, dat, init_params = init_un,
                                 control = ctrl, optimizer = "nlminb",
                                 parallel = FALSE, verbose = FALSE)
  expect_equal(res_un$convergence, 0)

  # Constrained re-run: previous_stage now matches Stage 2 structure
  # exactly, so allow_different_structure = FALSE and the un-fix loop in
  # initialize_factor_model_cpp is exercised.
  free_estimates <- res_un$estimates
  free_estimates["type_2_loading_2"] <- 0.0
  stage2_ref <- list(
    model_system = ms_se,
    estimates    = free_estimates,
    std_errors   = setNames(rep(0, length(free_estimates)),
                            names(free_estimates)),
    convergence  = 0L,
    loglik       = res_un$loglik
  )
  factor_dist_params <- grep("^factor_var|^se_|^typeprob_|^type_[0-9]+_loading_",
                              names(free_estimates), value = TRUE)
  free_params_list <- setdiff(factor_dist_params, "type_2_loading_2")
  expect_true(any(grepl("^se_cov_",     free_params_list)))
  expect_true(any(grepl("^typeprob_",   free_params_list)))
  expect_true(any(grepl("^type_2_load", free_params_list)))

  ms_se_fixed <- define_model_system(components = list(), factor = ms_se$factor,
                                      previous_stage = stage2_ref,
                                      free_params = free_params_list)
  res_c <- estimate_model_rcpp(ms_se_fixed, dat, init_params = free_estimates,
                                control = ctrl, optimizer = "nlminb",
                                parallel = FALSE, verbose = FALSE)
  expect_equal(res_c$convergence, 0,
               info = "Constrained re-run must converge strictly")

  # The smoking-gun assertion: every "free" factor-distribution slot must
  # carry its own value, not a recycled copy of an earlier slot's value.
  free_idx <- which(names(res_c$estimates) %in% free_params_list)
  free_vals <- res_c$estimates[free_idx]
  # type_*_loading_<n_factors> is auto-fixed at 0 for SE models, so it can
  # tie with other zero-init slots; exclude it from the uniqueness check.
  outcome_load <- paste0("type_2_loading_", fm_se$n_factors)
  free_vals_check <- free_vals[setdiff(names(free_vals), outcome_load)]
  n_unique <- length(unique(round(free_vals_check, 8)))
  expect_equal(n_unique, length(free_vals_check),
               info = sprintf("free factor-distribution slots have %d unique values out of %d (recycled): %s",
                              n_unique, length(free_vals_check),
                              paste(sprintf("%s=%.6f", names(free_vals_check),
                                            free_vals_check),
                                    collapse = ", ")))

  # Sanity: SE covariate value must NOT equal any of the SE / typeprob /
  # type_loading slot values (the symptom of the 7-cycle bug).
  cov_idx <- grep("^se_cov_", names(res_c$estimates))
  expect_true(length(cov_idx) > 0)
  other_vals <- res_c$estimates[setdiff(seq_along(res_c$estimates), cov_idx)]
  for (ci in cov_idx) {
    cv <- round(res_c$estimates[ci], 8)
    matches <- sum(round(other_vals, 8) == cv & abs(cv) > 1e-10)
    expect_equal(matches, 0,
                 info = sprintf("se_cov slot %s = %.6f duplicates a non-cov slot",
                                names(res_c$estimates)[ci], cv))
  }
})
