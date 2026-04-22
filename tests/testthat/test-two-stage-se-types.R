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
