# Regression test for the previous_stage / init_params stale-name path
# (factorana 1.3.2). Mirrors the warn-and-skip philosophy applied to
# define_model_component() in 1.3.1: previous_stage anchors that carry
# parameter names not present in the current model used to detonate
# inside setup_parameter_constraints()'s per-param branch logic when
# `param_metadata$types[i]` was NA (init_params longer than the current
# model's parameter vector). Now we warn once and skip those positions.

test_that("setup_parameter_constraints warns and skips stale init_params names", {
  set.seed(1)
  n <- 100
  dat <- data.frame(intercept = 1, y1 = rnorm(n), y2 = rnorm(n), eval = 1L)
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = 1,
                                 evaluation_indicator = "eval")
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = NA_real_,
                                 evaluation_indicator = "eval")
  ms <- define_model_system(components = list(mc1, mc2), factor = fm)

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  metadata <- factorana:::build_parameter_metadata(ms)

  # Inject a stale name into init_params and extend param_metadata mismatch.
  stale_init <- c(init$init_params,
                  m_dropped_intercept = 999.0,
                  m_dropped_sigma = 0.5)

  # setup_parameter_constraints should emit one warning naming the stale slots
  # and return without erroring.
  warns <- character(0)
  out <- withCallingHandlers(
    factorana:::setup_parameter_constraints(ms, stale_init, metadata,
                                             init$factor_variance_fixed,
                                             verbose = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("Skipping \\d+ previous_stage / init_params name", warns)),
              info = paste("Captured warnings:", paste(warns, collapse = " | ")))
  expect_true(any(grepl("m_dropped_intercept", warns)))
  expect_true(any(grepl("m_dropped_sigma", warns)))
  # The valid free_idx must NOT include the stale positions (which are past
  # the end of param_metadata$names).
  expect_true(all(out$free_idx <= length(metadata$names)))
})

test_that("estimate_model_rcpp drops stale previous_stage names without C++ size mismatch", {
  # Reproducer for the 1.3.2 bug report: 1.3.2 warned about stale names
  # inside setup_parameter_constraints() but still handed the unfiltered
  # init_params to the C++ side, which errored with "Fixed values size
  # mismatch". 1.3.3 reconciles full_init_params with the current model's
  # canonical name list BEFORE the C++ call.
  skip_on_cran()

  set.seed(4)
  N <- 400
  df_full <- data.frame(
    id = seq_len(N),
    m1 = sample(1:4, N, TRUE), m2 = sample(1:4, N, TRUE),
    m3 = sample(1:4, N, TRUE), m4 = sample(1:4, N, TRUE),
    eval_m1 = 1L, eval_m2 = 1L, eval_m3 = 1L, eval_m4 = 1L
  )
  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Stage 1: 4 items, all eval = 1.
  comps1 <- lapply(1:4, function(i) define_model_component(
    name = paste0("m", i), data = df_full, outcome = paste0("m", i),
    factor = fm, covariates = NULL, model_type = "oprobit", num_choices = 4,
    loading_normalization = if (i == 1) 1 else NA_real_,
    evaluation_indicator = paste0("eval_m", i)))
  ms1 <- define_model_system(components = comps1, factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)
  res1 <- estimate_model_rcpp(ms1, df_full, control = ctrl, optimizer = "nlminb",
                               parallel = FALSE, verbose = FALSE)
  expect_equal(res1$convergence, 0)

  # Stage 2: same items, but m4 has all-zero eval -> 1.3.1 drops it from
  # the component list. The previous_stage anchor (Stage 1 estimates) still
  # carries m4_* names; those must be dropped by reconciliation, not crash
  # the C++ side.
  df_s2 <- df_full; df_s2$eval_m4 <- 0L
  comps2 <- suppressWarnings(lapply(1:4, function(i) define_model_component(
    name = paste0("m", i), data = df_s2, outcome = paste0("m", i),
    factor = fm, covariates = NULL, model_type = "oprobit", num_choices = 4,
    loading_normalization = if (i == 1) 1 else NA_real_,
    evaluation_indicator = paste0("eval_m", i))))
  ms2 <- define_model_system(components = comps2, factor = fm)
  expect_equal(length(ms2$components), 3L)

  prev <- list(model_system = ms2, estimates = res1$estimates,
                std_errors = res1$std_errors, convergence = 0L, loglik = 0)

  # Stage 3: build a model_system that anchors on `prev`. The reported
  # bug is that the C++ side errors with "Fixed values size mismatch"
  # because the unfiltered prev anchor still carries m4_* names.
  # We free `factor_var_1` so the optimizer has at least one free
  # parameter to move (otherwise nlminb errors on an empty start vector,
  # which is a separate degenerate case, not the bug under test).
  fm_s2 <- define_factor_model(n_factors = 1, n_types = 1)
  ms3 <- define_model_system(components = list(), factor = fm_s2,
                              previous_stage = prev,
                              free_params = "factor_var_1")
  init3 <- initialize_parameters(ms3, df_s2, verbose = FALSE)$init_params

  warns <- character(0)
  res3 <- withCallingHandlers(
    estimate_model_rcpp(ms3, df_s2, init_params = init3, control = ctrl,
                        optimizer = "nlminb", parallel = FALSE, verbose = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(res3$convergence, 0,
               info = "Stage 3 must converge (no 'Fixed values size mismatch')")
  expect_true(any(grepl("Dropping .* name\\(s\\) not present", warns)) ||
              any(grepl("Skipping .* previous_stage", warns)),
              info = paste("Expected stale-name warning. Got:",
                           paste(warns, collapse = " | ")))
  # The C++ side did NOT error with "Fixed values size mismatch".
  expect_false(any(grepl("Fixed values size mismatch", warns)))
  # Result estimates should NOT contain m4_* names.
  expect_false(any(grepl("^m4_", names(res3$estimates))),
               info = paste("Unexpected m4_* in estimates:",
                            paste(grep("^m4_", names(res3$estimates), value = TRUE),
                                  collapse = ", ")))
})

test_that("Stale-name path: clean init_params (no stale) does not warn", {
  set.seed(2)
  n <- 100
  dat <- data.frame(intercept = 1, y1 = rnorm(n), y2 = rnorm(n), eval = 1L)
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = 1,
                                 evaluation_indicator = "eval")
  mc2 <- define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = NA_real_,
                                 evaluation_indicator = "eval")
  ms <- define_model_system(components = list(mc1, mc2), factor = fm)

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  metadata <- factorana:::build_parameter_metadata(ms)

  warns <- character(0)
  out <- withCallingHandlers(
    factorana:::setup_parameter_constraints(ms, init$init_params, metadata,
                                             init$factor_variance_fixed,
                                             verbose = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("Skipping .* previous_stage", warns)),
               info = paste("Unexpected warnings:", paste(warns, collapse = " | ")))
})
