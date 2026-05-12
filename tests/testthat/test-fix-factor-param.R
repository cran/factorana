# Tests for fix_factor_param(): factor-distribution parameter fixing
# at model-definition time.

test_that("fix_factor_param stores constraint on the factor model", {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "SE_linear")
  expect_null(fm$fixed_params)

  fm <- fix_factor_param(fm, "type_2_loading_2", 0.0)
  expect_equal(fm$fixed_params, c(type_2_loading_2 = 0.0))

  fm <- fix_factor_param(fm, "type_2_loading_1", 0.3)
  expect_equal(fm$fixed_params,
               c(type_2_loading_2 = 0.0, type_2_loading_1 = 0.3))
})

test_that("fix_factor_param accepts batch named-vector form", {
  fm <- define_factor_model(n_factors = 3, n_types = 2,
                            factor_structure = "SE_linear")
  fm <- fix_factor_param(fm, c(type_2_loading_2 = 0.0,
                                type_2_loading_3 = 0.0))
  expect_equal(fm$fixed_params,
               c(type_2_loading_2 = 0.0, type_2_loading_3 = 0.0))
})

test_that("fix_factor_param replaces existing values on second call", {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "SE_linear")
  fm <- fix_factor_param(fm, "typeprob_2_intercept", 0.5)
  fm <- fix_factor_param(fm, "typeprob_2_intercept", -0.2)
  expect_equal(unname(fm$fixed_params["typeprob_2_intercept"]), -0.2)
  expect_equal(length(fm$fixed_params), 1L)
})

test_that("fix_factor_param with NA value unfixes the parameter", {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "SE_linear")
  fm <- fix_factor_param(fm, c(typeprob_2_intercept = 0.3,
                                type_2_loading_1 = 0.4))
  fm <- fix_factor_param(fm, "typeprob_2_intercept", NA_real_)
  expect_equal(fm$fixed_params, c(type_2_loading_1 = 0.4))

  fm <- fix_factor_param(fm, "type_2_loading_1", NA_real_)
  expect_null(fm$fixed_params)
})

test_that("fix_factor_param errors on invalid parameter name", {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "SE_linear")
  expect_error(fix_factor_param(fm, "type_2_loading_5", 0),
               "is not a valid factor-distribution parameter name")
  expect_error(fix_factor_param(fm, "typeprob2_intercept", 0),
               "is not a valid factor-distribution parameter name")
  # se_cov_X is invalid because no se_covariates were specified.
  expect_error(fix_factor_param(fm, "se_cov_X", 0),
               "is not a valid factor-distribution parameter name")
})

test_that("fix_factor_param accepts se_cov_* names when se_covariates set", {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "SE_linear",
                            se_covariates = c("X", "Z"))
  fm <- fix_factor_param(fm, "se_cov_X", 0.5)
  expect_equal(unname(fm$fixed_params["se_cov_X"]), 0.5)
})

test_that("fix_factor_param accepts factor_mean_* names when factor_covariates set", {
  fm <- define_factor_model(n_factors = 2,
                            factor_covariates = c("age"))
  fm <- fix_factor_param(fm, "factor_mean_1_age", 0.0)
  expect_equal(unname(fm$fixed_params["factor_mean_1_age"]), 0.0)
})

test_that("fix_factor_param: outcome-factor type loading is idempotent at 0", {
  fm <- define_factor_model(n_factors = 3, n_types = 2,
                            factor_structure = "SE_linear")
  # n_factors = 3 => outcome factor is index 3
  expect_silent(fm <- fix_factor_param(fm, "type_2_loading_3", 0.0))
  expect_error(fix_factor_param(fm, "type_2_loading_3", 0.5),
               "auto-fixed at 0")
})

test_that("fix_factor_param errors on duplicate names in batch input", {
  fm <- define_factor_model(n_factors = 2, n_types = 2,
                            factor_structure = "SE_linear")
  expect_error(
    fix_factor_param(fm, c("type_2_loading_1", "type_2_loading_1"),
                      c(0.1, 0.2)),
    "Duplicate parameter name"
  )
})

# =============================================================================
# End-to-end integration: fixed parameters propagate through estimation
# =============================================================================
test_that("fix_factor_param: fixed parameter has SE = 0 and recovers value", {
  skip_on_cran()

  set.seed(2027)
  n <- 600
  f1 <- rnorm(n, 0, 1)
  eps <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.6 * f1 + eps
  dat <- data.frame(
    intercept = 1, eval = 1L,
    Y_t1_m1 = f1 + rnorm(n, 0, 0.6),
    Y_t1_m2 = 0.9 * f1 + rnorm(n, 0, 0.65),
    Y_t1_m3 = 1.1 * f1 + rnorm(n, 0, 0.55),
    Y_t2_m1 = f2 + rnorm(n, 0, 0.6),
    Y_t2_m2 = 0.9 * f2 + rnorm(n, 0, 0.65),
    Y_t2_m3 = 1.1 * f2 + rnorm(n, 0, 0.55)
  )

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  s1 <- estimate_model_rcpp(dyn$model_system, dat, control = ctrl,
                             optimizer = "nlminb", parallel = FALSE,
                             verbose = FALSE)
  expect_equal(s1$convergence, 0)

  prev <- build_dynamic_previous_stage(dyn, s1, dat, anchor_period = 1L)

  fm_se <- define_factor_model(n_factors = 2, n_types = 2,
                                factor_structure = "SE_linear")
  # Fix the input-factor type loading at 0 and check it stays exactly 0.
  fm_se <- fix_factor_param(fm_se, "type_2_loading_1", 0.0)
  ms_se <- define_model_system(components = list(), factor = fm_se,
                                previous_stage = prev)

  res <- estimate_model_rcpp(ms_se, dat, control = ctrl,
                              optimizer = "nlminb", parallel = FALSE,
                              verbose = FALSE)
  expect_equal(res$convergence, 0,
               info = "Estimation must converge with type_2_loading_1 fixed at 0")

  # Fixed parameter is exactly at the user-specified value with SE = 0.
  expect_equal(unname(res$estimates["type_2_loading_1"]), 0.0,
               tolerance = 1e-12)
  expect_equal(unname(res$std_errors["type_2_loading_1"]), 0.0,
               tolerance = 1e-12)
  expect_true("type_2_loading_1" %in% names(res$estimates))
})

# =============================================================================
# Conflict semantics with previous_stage / free_params
# =============================================================================
test_that("fix_factor_param wins over free_params with a single warning", {
  skip_on_cran()

  # Three-stage chain (mirrors the MH-trap regime workflow):
  #   Stage 1: dynamic measurement (independent factors, no types)
  #   Stage 2: SE_linear + n_types=2 (unconstrained) — produces all the
  #            SE / typeprob / type_loading parameter names we need.
  #   Stage 3: SE_linear + n_types=2 with previous_stage = Stage 2.
  #            Use fix_factor_param() AND list the same name in
  #            free_params; fix_factor_param must win (the warning
  #            fires) and the parameter must remain at the user value.

  set.seed(2028)
  n <- 600
  f1 <- rnorm(n, 0, 1)
  eps <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.6 * f1 + 0.3 * (runif(n) < 0.5) + eps
  dat <- data.frame(
    intercept = 1, eval = 1L,
    Y_t1_m1 = f1 + rnorm(n, 0, 0.6),
    Y_t1_m2 = 0.9 * f1 + rnorm(n, 0, 0.65),
    Y_t1_m3 = 1.1 * f1 + rnorm(n, 0, 0.55),
    Y_t2_m1 = f2 + rnorm(n, 0, 0.6),
    Y_t2_m2 = 0.9 * f2 + rnorm(n, 0, 0.65),
    Y_t2_m3 = 1.1 * f2 + rnorm(n, 0, 0.55)
  )

  dyn <- define_dynamic_measurement(
    data = dat, items = c("m1", "m2", "m3"),
    period_prefixes = c("Y_t1_", "Y_t2_"),
    model_type = "linear", evaluation_indicator = "eval"
  )
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  s1 <- estimate_model_rcpp(dyn$model_system, dat, control = ctrl,
                             optimizer = "nlminb", parallel = FALSE,
                             verbose = FALSE)
  expect_equal(s1$convergence, 0)
  prev <- build_dynamic_previous_stage(dyn, s1, dat, anchor_period = 1L)

  # Stage 2: unconstrained SE_linear + n_types=2.
  fm2 <- define_factor_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_linear")
  ms2 <- define_model_system(components = list(), factor = fm2,
                              previous_stage = prev)
  init2 <- initialize_parameters(ms2, dat, verbose = FALSE)$init_params
  init2["se_linear_1"]          <- 0.5
  init2["se_residual_var"]      <- 0.5
  init2["se_intercept_type_2"]  <- 0.3
  init2["typeprob_2_intercept"] <- 0.1
  init2["type_2_loading_1"]     <- 0.2
  s2 <- estimate_model_rcpp(ms2, dat, init_params = init2, control = ctrl,
                             optimizer = "nlminb", parallel = FALSE,
                             verbose = FALSE)
  expect_equal(s2$convergence, 0)

  # Stage 3: matching SE_linear + n_types=2 with fix_factor_param fighting
  # against free_params. The fix wins and a warning fires.
  fm3 <- define_factor_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_linear")
  fm3 <- fix_factor_param(fm3, "typeprob_2_intercept", 0.0)
  free_params_list <- c("factor_var_1", "typeprob_2_intercept",
                        "se_intercept", "se_linear_1", "se_residual_var",
                        "se_intercept_type_2", "type_2_loading_1")

  ms3 <- define_model_system(components = list(), factor = fm3,
                              previous_stage = s2,
                              free_params = free_params_list)
  init3 <- s2$estimates

  # Capture all warnings via withCallingHandlers; expect_warning only
  # matches the first warning, but two distinct warnings fire here
  # (free_params override AND value override vs previous_stage).
  warns <- character(0)
  res3 <- withCallingHandlers(
    estimate_model_rcpp(ms3, dat, init_params = init3, control = ctrl,
                        optimizer = "nlminb", parallel = FALSE,
                        verbose = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("fix_factor_param.*overrides free_params", warns)),
              info = paste("Captured warnings:", paste(warns, collapse = " | ")))
  expect_equal(res3$convergence, 0)
  expect_equal(unname(res3$estimates["typeprob_2_intercept"]), 0.0,
               tolerance = 1e-12)
  expect_equal(unname(res3$std_errors["typeprob_2_intercept"]), 0.0,
               tolerance = 1e-12)
})
