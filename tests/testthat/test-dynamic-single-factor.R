# Dynamic single-factor test: one latent construct measured at two time
# points, coupled by a linear structural equation.
#
# This test uses the define_dynamic_measurement() / build_dynamic_previous_stage()
# wrapper functions, which encapsulate:
#   - Stage 1: 2-factor independent measurement system with
#     equality_constraints tying loadings and sigmas across periods and
#     period-specific intercepts.
#   - Stage 2 bridge: build a dummy previous_stage that carries the
#     anchor-period (wave 1) intercepts into every factor slot.
#
# Then Stage 2 SE_linear recovers alpha, beta, sigma_eps^2, and Var(f_1).

VERBOSE <- Sys.getenv("FACTORANA_TEST_VERBOSE", "FALSE") == "TRUE"


# ---- DGP helper -------------------------------------------------------------

.simulate_dynamic_single_factor_dgp <- function(
    n            = 1500,
    seed         = 41,
    true_var_f1  = 1.0,
    true_alpha   = 0.4,
    true_beta    = 0.6,
    true_sigma_e = sqrt(0.5),
    item_int     = c(1.5, 1.0, 0.8),
    item_load    = c(1.0, 0.9, 1.1),
    item_sigma   = c(0.7, 0.75, 0.65)) {

  set.seed(seed)

  f1  <- rnorm(n, 0, sqrt(true_var_f1))
  eps <- rnorm(n, 0, true_sigma_e)
  f2  <- true_alpha + true_beta * f1 + eps

  gen_Y <- function(f, i) {
    item_int[i] + item_load[i] * f + rnorm(length(f), 0, item_sigma[i])
  }

  dat_wide <- data.frame(
    id        = seq_len(n),
    intercept = 1,
    eval      = 1L,
    Y_t1_m1 = gen_Y(f1, 1), Y_t1_m2 = gen_Y(f1, 2), Y_t1_m3 = gen_Y(f1, 3),
    Y_t2_m1 = gen_Y(f2, 1), Y_t2_m2 = gen_Y(f2, 2), Y_t2_m3 = gen_Y(f2, 3)
  )

  list(
    wide = dat_wide,
    true = list(
      var_f1     = true_var_f1,
      alpha      = true_alpha,
      beta       = true_beta,
      sigma_e2   = true_sigma_e^2,
      item_int   = item_int,
      item_load  = item_load,
      item_sigma = item_sigma
    )
  )
}


# =============================================================================
# TEST 1: Stage 1 tied measurement recovers DGP intercepts on wave 1
# =============================================================================
test_that("define_dynamic_measurement: Stage 1 recovers wave-1 tau_m, shared loadings/sigmas", {
  skip_on_cran()

  sim   <- .simulate_dynamic_single_factor_dgp(n = 1500, seed = 41)
  truth <- sim$true

  dyn <- define_dynamic_measurement(
    data                 = sim$wide,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  expect_s3_class(dyn, "dynamic_measurement")

  # Equality constraints: 2 loadings (items m2, m3) + 3 sigmas = 5 groups
  expect_equal(length(dyn$equality_constraints), 5L)

  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  s1 <- estimate_model_rcpp(
    dyn$model_system, sim$wide, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(s1$convergence, 0,
               info = "Stage 1 must converge strictly")

  est <- s1$estimates

  # Wave-1 intercepts recover DGP (E[f_1] = 0 by convention)
  for (i in 1:3) {
    nm <- paste0("Y_t1_m", i, "_intercept")
    expect_equal(unname(est[nm]), truth$item_int[i], tolerance = 0.10,
                 info = sprintf("%s: true=%.3f est=%.3f",
                                nm, truth$item_int[i], est[nm]))
  }

  # Tied loadings and sigmas recover DGP
  expect_equal(unname(est["Y_t1_m2_loading_1"]), truth$item_load[2], tolerance = 0.05)
  expect_equal(unname(est["Y_t1_m3_loading_1"]), truth$item_load[3], tolerance = 0.05)
  for (i in 1:3) {
    nm <- paste0("Y_t1_m", i, "_sigma")
    expect_equal(unname(est[nm]), truth$item_sigma[i], tolerance = 0.05)
  }

  # Equality constraints hold exactly in the estimates
  expect_equal(unname(est["Y_t1_m2_loading_1"]),
               unname(est["Y_t2_m2_loading_2"]), tolerance = 1e-10)
  expect_equal(unname(est["Y_t1_m1_sigma"]),
               unname(est["Y_t2_m1_sigma"]),     tolerance = 1e-10)
})


# =============================================================================
# TEST 2: End-to-end: wrapper plus Stage 2 SE_linear recovers structural params
# =============================================================================
test_that("define_dynamic_measurement + Stage 2 SE_linear recovers alpha, beta, sigma_eps, var_f1", {
  skip_on_cran()

  sim   <- .simulate_dynamic_single_factor_dgp(n = 1500, seed = 41)
  truth <- sim$true
  ctrl  <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  # Stage 1 via the wrapper
  dyn <- define_dynamic_measurement(
    data                 = sim$wide,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  s1 <- estimate_model_rcpp(
    dyn$model_system, sim$wide, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(s1$convergence, 0)

  # Build the Stage 2 previous_stage bridge via the wrapper
  dummy <- build_dynamic_previous_stage(dyn, s1, sim$wide, anchor_period = 1L)

  # Stage 2: SE_linear
  fm_s2 <- define_factor_model(n_factors = 2, n_types = 1,
                                factor_structure = "SE_linear")
  ms_s2 <- define_model_system(components = list(), factor = fm_s2,
                                previous_stage = dummy)

  init_s2 <- initialize_parameters(ms_s2, sim$wide, verbose = FALSE)
  init_s2$init_params["factor_var_1"]    <- unname(dummy$estimates["factor_var_1"])
  init_s2$init_params["se_intercept"]    <- 0.0
  init_s2$init_params["se_linear_1"]     <- 0.5
  init_s2$init_params["se_residual_var"] <- 0.5

  r2 <- estimate_model_rcpp(
    ms_s2, sim$wide, init_params = init_s2$init_params, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(r2$convergence, 0)

  est <- r2$estimates
  se  <- r2$std_errors

  # The load-bearing assertion: alpha recovers.
  expect_equal(unname(est["se_intercept"]), truth$alpha, tolerance = 0.12,
               info = sprintf("se_intercept (alpha): true=%.3f est=%.3f se=%.3f",
                              truth$alpha, est["se_intercept"], se["se_intercept"]))

  # Other structural params
  expect_equal(unname(est["se_linear_1"]),     truth$beta,     tolerance = 0.10)
  expect_equal(unname(est["se_residual_var"]), truth$sigma_e2, tolerance = 0.15)
  expect_equal(unname(est["factor_var_1"]),    truth$var_f1,   tolerance = 0.15)

  new_free <- c("factor_var_1", "se_intercept", "se_linear_1", "se_residual_var")
  expect_true(all(se[new_free] > 0))

  if (VERBOSE) {
    cat("\n=== Stage 2 SE_linear recovery (via wrapper) ===\n")
    cat(sprintf("%-22s %10s %10s %10s\n", "param", "true", "est", "se"))
    for (p in new_free) {
      tr <- switch(p,
                   factor_var_1    = truth$var_f1,
                   se_intercept    = truth$alpha,
                   se_linear_1     = truth$beta,
                   se_residual_var = truth$sigma_e2)
      cat(sprintf("  %-22s %+10.4f %+10.4f %10.4f\n",
                  p, tr, est[p], se[p]))
    }
  }
})


# ---- DGP helper with types (types shift f_2 mean only) ----------------------

.simulate_dynamic_single_factor_dgp_types <- function(
    n            = 3000,
    seed         = 53,
    true_var_f1  = 1.0,
    true_alpha   = 0.4,
    true_beta    = 0.6,
    true_sigma_e = sqrt(0.5),
    true_alpha_t2 = 0.8,
    typeprob_t2   = 0.3,
    type_load_t2  = 0.4,
    item_int      = c(1.5, 1.0, 0.8),
    item_load     = c(1.0, 0.9, 1.1),
    item_sigma    = c(0.7, 0.75, 0.65)) {

  set.seed(seed)

  f1 <- rnorm(n, 0, sqrt(true_var_f1))
  p_t2 <- plogis(typeprob_t2 + type_load_t2 * f1)
  t2   <- as.integer(runif(n) < p_t2)
  eps  <- rnorm(n, 0, true_sigma_e)
  f2   <- true_alpha + true_beta * f1 + true_alpha_t2 * t2 + eps

  gen_Y <- function(f, i) {
    item_int[i] + item_load[i] * f + rnorm(length(f), 0, item_sigma[i])
  }

  dat_wide <- data.frame(
    id        = seq_len(n),
    intercept = 1,
    eval      = 1L,
    Y_t1_m1 = gen_Y(f1, 1), Y_t1_m2 = gen_Y(f1, 2), Y_t1_m3 = gen_Y(f1, 3),
    Y_t2_m1 = gen_Y(f2, 1), Y_t2_m2 = gen_Y(f2, 2), Y_t2_m3 = gen_Y(f2, 3)
  )

  list(
    wide = dat_wide,
    true = list(
      var_f1       = true_var_f1,
      alpha        = true_alpha,
      beta         = true_beta,
      sigma_e2     = true_sigma_e^2,
      alpha_t2     = true_alpha_t2,
      typeprob_t2  = typeprob_t2,
      type_load_t2 = type_load_t2
    )
  )
}


# =============================================================================
# TEST 3: End-to-end with n_types = 2 at Stage 2 (types shift f_2 mean only)
#
# The wrapper is type-agnostic at Stage 1 (hardcoded n_types = 1L). Stage 2
# adds n_types = 2 via define_factor_model(); the SE_linear model picks up
# the type-specific intercept se_intercept_type_2, the typeprob intercept,
# and the type loading on the input factor. This test verifies that the
# wrapper plus types path recovers the well-identified structural
# parameters, and that all free parameters have positive standard errors.
# =============================================================================
test_that("define_dynamic_measurement + Stage 2 SE_linear + n_types=2 recovers structural params", {
  skip_on_cran()

  sim   <- .simulate_dynamic_single_factor_dgp_types(n = 3000, seed = 53)
  truth <- sim$true
  ctrl  <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  dyn <- define_dynamic_measurement(
    data                 = sim$wide,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )
  s1 <- estimate_model_rcpp(
    dyn$model_system, sim$wide, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(s1$convergence, 0)

  dummy <- build_dynamic_previous_stage(dyn, s1, sim$wide, anchor_period = 1L)

  fm_s2 <- define_factor_model(n_factors = 2, n_types = 2,
                                factor_structure = "SE_linear")
  ms_s2 <- define_model_system(components = list(), factor = fm_s2,
                                previous_stage = dummy)

  # Neutral init (not at truth). Earlier manual testing confirmed the
  # mode is robust: three different inits (neutral, alpha-absorbs-shift,
  # type-absorbs-shift) all converge to the same log-likelihood.
  init_s2 <- initialize_parameters(ms_s2, sim$wide, verbose = FALSE)
  init_s2$init_params["factor_var_1"]         <- unname(dummy$estimates["factor_var_1"])
  init_s2$init_params["se_intercept"]         <- 0.0
  init_s2$init_params["se_linear_1"]          <- 0.5
  init_s2$init_params["se_residual_var"]      <- 0.5
  init_s2$init_params["se_intercept_type_2"]  <- 0.3
  init_s2$init_params["typeprob_2_intercept"] <- 0.0
  init_s2$init_params["type_2_loading_1"]     <- 0.0

  r2 <- estimate_model_rcpp(
    ms_s2, sim$wide, init_params = init_s2$init_params, control = ctrl,
    optimizer = "nlminb", parallel = FALSE, verbose = FALSE
  )
  expect_equal(r2$convergence, 0,
               info = "Stage 2 (SE_linear + n_types=2) must converge strictly")

  est <- r2$estimates
  se  <- r2$std_errors

  # Well-identified structural parameters (tight tolerances).
  expect_equal(unname(est["se_linear_1"]),         truth$beta,     tolerance = 0.10,
               info = sprintf("se_linear_1 (beta): true=%.3f est=%.3f",
                              truth$beta, est["se_linear_1"]))
  expect_equal(unname(est["se_intercept_type_2"]), truth$alpha_t2, tolerance = 0.20,
               info = sprintf("se_intercept_type_2 (alpha_type): true=%.3f est=%.3f",
                              truth$alpha_t2, est["se_intercept_type_2"]))
  expect_equal(unname(est["se_residual_var"]),     truth$sigma_e2, tolerance = 0.15,
               info = sprintf("se_residual_var: true=%.3f est=%.3f",
                              truth$sigma_e2, est["se_residual_var"]))
  expect_equal(unname(est["factor_var_1"]),        truth$var_f1,   tolerance = 0.15,
               info = sprintf("factor_var_1: true=%.3f est=%.3f",
                              truth$var_f1, est["factor_var_1"]))

  # se_intercept (alpha) trades off with typeprob_2_intercept and
  # type_2_loading_1 through the mean decomposition
  # E[f_2] = se_intercept + se_intercept_type_2 * Pr(type=2). The SUM
  # is tightly identified; the individual components are identified but
  # noisier at moderate n. Empirically at n = 3000, SE scales roughly
  # as 1/sqrt(n) and all z-scores stay under ~2 across sample sizes.
  # Tolerances here are set at roughly 1 SE at n = 3000.
  expect_equal(unname(est["se_intercept"]),         truth$alpha,       tolerance = 0.25,
               info = sprintf("se_intercept (alpha): true=%.3f est=%.3f se=%.3f",
                              truth$alpha, est["se_intercept"], se["se_intercept"]))
  expect_equal(unname(est["typeprob_2_intercept"]), truth$typeprob_t2, tolerance = 0.30,
               info = sprintf("typeprob_2_intercept: true=%.3f est=%.3f se=%.3f",
                              truth$typeprob_t2, est["typeprob_2_intercept"],
                              se["typeprob_2_intercept"]))
  expect_equal(unname(est["type_2_loading_1"]),     truth$type_load_t2, tolerance = 0.50,
               info = sprintf("type_2_loading_1: true=%.3f est=%.3f se=%.3f",
                              truth$type_load_t2, est["type_2_loading_1"],
                              se["type_2_loading_1"]))

  # Positive SE on every new free parameter (sanity: Fisher information
  # is not singular).
  new_free <- c("factor_var_1", "se_intercept", "se_linear_1",
                "se_residual_var", "se_intercept_type_2",
                "typeprob_2_intercept", "type_2_loading_1")
  expect_true(all(se[new_free] > 0),
              info = "All Stage 2 free params must have positive std_errors")

  if (VERBOSE) {
    cat("\n=== Stage 2 SE_linear + n_types=2 recovery (via wrapper) ===\n")
    cat(sprintf("%-22s %10s %10s %10s %10s\n", "param", "true", "est", "se", "z"))
    for (p in new_free) {
      tr <- switch(p,
                   factor_var_1         = truth$var_f1,
                   se_intercept         = truth$alpha,
                   se_linear_1          = truth$beta,
                   se_residual_var      = truth$sigma_e2,
                   se_intercept_type_2  = truth$alpha_t2,
                   typeprob_2_intercept = truth$typeprob_t2,
                   type_2_loading_1     = truth$type_load_t2)
      zv <- if (se[p] > 0) (est[p] - tr) / se[p] else NA
      cat(sprintf("  %-22s %+10.4f %+10.4f %10.4f %+10.2f\n",
                  p, tr, est[p], se[p], zv))
    }
  }
})


# =============================================================================
# TEST 4: structural check for model_type = "oprobit"
#
# Ordered probit needs a different tying pattern than linear. factorana
# parameterises cutpoints as increments:
#     cutpoint_k = thresh_1 + thresh_2 + ... + thresh_k
# so thresh_1 plays the role of a linear intercept (location) and the
# later increments play the role of sigma (scale / category spacing).
# The wrapper ties only the increments across periods, leaving thresh_1
# period-specific so a wave-to-wave shift in the latent factor mean can
# be absorbed. The wrapper also silently strips the default
# `intercept` covariate for oprobit (factorana rejects an intercept
# covariate on oprobit components).
#
# This is a STRUCTURAL test: it checks the constructed model_system and
# the bridge object, but does not run optimisation. Recovery of the
# ordered-probit dynamic model is empirically fragile at moderate n
# (see the patch README for the Mental Health Trap findings); we track
# estimation-level testing separately.
# =============================================================================
test_that("oprobit wrapper ties only threshold increments and strips intercept", {
  n_items      <- 3L
  n_categories <- 4L     # => 3 thresholds per item

  # Minimal data frame: the wrapper only reads column names and types.
  set.seed(1); n <- 100
  mk_col <- function() sample.int(n_categories, n, replace = TRUE)
  dat <- data.frame(
    intercept = 1, eval = 1L,
    Y_t1_m1 = mk_col(), Y_t1_m2 = mk_col(), Y_t1_m3 = mk_col(),
    Y_t2_m1 = mk_col(), Y_t2_m2 = mk_col(), Y_t2_m3 = mk_col()
  )

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = paste0("m", 1:n_items),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "oprobit",
    n_categories         = n_categories,
    covariates           = "intercept",   # should be stripped silently
    evaluation_indicator = "eval"
  )

  # ---- covariates auto-stripped for oprobit ----
  expect_true(is.null(dyn$covariates) || length(dyn$covariates) == 0L,
              info = "oprobit default covariates should be NULL after stripping intercept")

  # ---- equality constraints ----
  # Expected: (n_items - 1) loading ties + n_items * (n_categories - 2)
  # threshold-increment ties (thresh_k for k = 2..K-1). thresh_1 is NOT tied.
  n_loading_ties <- n_items - 1L              # items 2..n, loading free on both factors
  n_thresh_ties  <- n_items * (n_categories - 2L)  # thresh_2 and thresh_3 for each item
  expected_n_eq  <- n_loading_ties + n_thresh_ties
  expect_equal(length(dyn$equality_constraints), expected_n_eq,
               info = sprintf("expected %d equality constraints, got %d",
                              expected_n_eq, length(dyn$equality_constraints)))

  # Every constraint group should be a 2-element character vector
  expect_true(all(vapply(dyn$equality_constraints,
                          function(g) is.character(g) && length(g) == 2L, logical(1))))

  # thresh_1 should NEVER appear in any tie group
  tied_names <- unlist(dyn$equality_constraints)
  expect_false(any(grepl("_thresh_1$", tied_names)),
               info = "thresh_1 must be period-specific (not tied)")

  # thresh_2 and thresh_3 SHOULD appear for every item
  for (i in seq_len(n_items)) {
    for (k in 2:(n_categories - 1L)) {
      expect_true(any(grepl(sprintf("_m%d_thresh_%d$", i, k), tied_names)),
                   info = sprintf("expected tie for m%d thresh_%d", i, k))
    }
  }

  # ---- Parameter layout through the optimiser's metadata ----
  md <- factorana:::build_parameter_metadata(dyn$model_system)
  # Every component should have a thresh_1 slot (period-specific).
  for (p in c("Y_t1_", "Y_t2_")) {
    for (i in seq_len(n_items)) {
      nm <- paste0(p, "m", i, "_thresh_1")
      expect_true(nm %in% md$names,
                   info = sprintf("%s must exist as a free parameter", nm))
    }
  }

  # ---- Bridge handles oprobit (no _intercept slot expected) ----
  # Build a fake Stage 1 result with plausible parameter values at the
  # layout build_parameter_metadata emits, so that
  # build_dynamic_previous_stage can run without estimation.
  fake_est <- rep(0.5, length(md$names))
  names(fake_est) <- md$names
  # Give plausible increasing thresholds so the sum is strictly monotone.
  for (p in c("Y_t1_", "Y_t2_")) {
    for (i in seq_len(n_items)) {
      fake_est[paste0(p, "m", i, "_thresh_1")] <- -1.0
      fake_est[paste0(p, "m", i, "_thresh_2")] <-  1.0
      fake_est[paste0(p, "m", i, "_thresh_3")] <-  1.0
    }
  }
  fake_result <- list(
    estimates = fake_est,
    std_errors = setNames(rep(0.1, length(fake_est)), names(fake_est)),
    convergence = 0L,
    loglik = 0.0
  )

  dummy <- build_dynamic_previous_stage(dyn, fake_result, dat)
  expect_true(is.list(dummy) && "estimates" %in% names(dummy))
  # Dummy must NOT contain any _intercept parameters (oprobit has none).
  expect_false(any(grepl("_intercept$", names(dummy$estimates))),
               info = "dummy previous_stage for oprobit must not contain _intercept params")
  # Dummy SHOULD contain thresh_1 slots for both periods, carrying the
  # anchor period's value into every slot.
  anchor_thresh1 <- fake_est["Y_t1_m1_thresh_1"]
  for (p in c("Y_t1_", "Y_t2_")) {
    for (i in seq_len(n_items)) {
      nm <- paste0(p, "m", i, "_thresh_1")
      expect_true(nm %in% names(dummy$estimates),
                   info = sprintf("%s must appear in dummy estimates", nm))
    }
  }
  # All m1 thresh_1 slots should equal the anchor value.
  m1_slots <- grep("m1_thresh_1$", names(dummy$estimates), value = TRUE)
  expect_true(length(m1_slots) == 2L)
  expect_equal(unname(dummy$estimates[m1_slots[1]]),
               unname(dummy$estimates[m1_slots[2]]),
               info = "anchor-period thresh_1 must be carried into every period slot")
})


# =============================================================================
# TEST 5: oprobit Stage 1 estimation converges cleanly with no recycling
# warnings.
#
# This exercises the actual Stage 1 fit for the oprobit dynamic wrapper.
# Historically (v1.1.6 and earlier) the C++ initialiser looked up
# threshold-parameter names by the wrong field name ("n_categories"
# instead of "num_choices"), so threshold equality constraints were
# silently dropped. The optimizer then ran over a larger free set than R
# thought; pmax/pmin in the saddle-escape path produced "fractionally
# recycled" warnings, and Stage 1 often landed at conv = 1 with factor
# variances contorted to absorb the period-mean drift.
#
# With the C++ fix, the R and C++ sides agree on the free-parameter set
# and Stage 1 converges strictly with no warnings. We also check that
# the threshold increments (thresh_k for k >= 2) come out exactly tied
# across periods, and that thresh_1 is period-specific with a
# consistent-sign shift reflecting the DGP mean drift.
# =============================================================================
test_that("oprobit Stage 1 fits cleanly with threshold equality constraints enforced", {
  skip_on_cran()

  set.seed(71); n <- 500
  item_load <- c(1.0, 0.9, 1.1)
  cuts      <- c(-1.0, 0.0, 1.0)
  f1 <- rnorm(n, 0, 1)
  f2 <- 0.4 + 0.6 * f1 + rnorm(n, 0, sqrt(0.4))   # positive period-mean drift
  gen_Y <- function(f, i) {
    as.integer(findInterval(item_load[i] * f + rnorm(length(f), 0, 1), cuts) + 1L)
  }
  dat <- data.frame(
    intercept = 1, eval = 1L,
    Y_t1_m1 = gen_Y(f1, 1), Y_t1_m2 = gen_Y(f1, 2), Y_t1_m3 = gen_Y(f1, 3),
    Y_t2_m1 = gen_Y(f2, 1), Y_t2_m2 = gen_Y(f2, 2), Y_t2_m3 = gen_Y(f2, 3)
  )

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "oprobit",
    n_categories         = 4L,
    evaluation_indicator = "eval"
  )

  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)

  # No recycling / fractional-recycled warnings during the fit.
  ws <- capture_warnings(
    s1 <- estimate_model_rcpp(
      dyn$model_system, dat, control = ctrl,
      optimizer = "nlminb", parallel = FALSE, verbose = FALSE
    )
  )
  expect_equal(s1$convergence, 0,
               info = "Stage 1 oprobit must converge strictly (conv = 0)")
  bad <- grep("recycled|recycling", ws, value = TRUE)
  expect_equal(length(bad), 0,
               info = paste("Unexpected recycling warnings:",
                            paste(bad, collapse = "; ")))

  # Threshold increments (k >= 2) tied exactly across periods.
  for (i in 1:3) {
    for (k in 2:3) {
      d1 <- unname(s1$estimates[paste0("Y_t1_m", i, "_thresh_", k)])
      d2 <- unname(s1$estimates[paste0("Y_t2_m", i, "_thresh_", k)])
      expect_equal(d1, d2, tolerance = 1e-10,
                   info = sprintf("m%d thresh_%d must be tied across periods", i, k))
    }
  }

  # thresh_1 period-specific; with a POSITIVE DGP mean drift in f_2, the
  # model's period-2 thresh_1 should be LOWER than period-1 (because
  # higher observed Y's push the cutpoints down when the factor is
  # forced to mean 0 by convention). Every item should show the same
  # sign of shift.
  shifts <- vapply(1:3, function(i) {
    unname(s1$estimates[paste0("Y_t2_m", i, "_thresh_1")] -
           s1$estimates[paste0("Y_t1_m", i, "_thresh_1")])
  }, numeric(1))
  expect_true(all(shifts < 0),
              info = sprintf("Expected negative thresh_1 shifts across items; got %s",
                             paste(sprintf("%+.3f", shifts), collapse = ", ")))
})


# =============================================================================
# TEST 6: Stage 1 FD check (linear) for the dynamic-measurement wrapper.
#
# Stage 1 is a 2-factor independent measurement model with
# equality_constraints tying loadings and sigmas across periods. The FD
# gradient and Hessian are checked at DGP-consistent parameter values
# with the equality constraints active. This path exposes the
# Hessian-accumulation fix in FactorModel::CalcLkhd: before the fix,
# the C++ Hessian computation iterated only over freeparlist and so
# missed cross-derivatives between a primary and its tied derived
# parameters, producing analytical values of ~0 at positions where FD
# reported magnitudes of 10 to 50. After the fix (iterating gradparlist
# plus symmetrising full_hessL before ExtractFreeHessian) the
# aggregation recovers the effective Hessian at primary positions and
# analytical / FD agree to ~1e-6.
# =============================================================================
test_that("dynamic wrapper Stage 1 (linear): FD gradient and Hessian match", {
  skip_on_cran()

  sim <- .simulate_dynamic_single_factor_dgp(n = 800, seed = 81)
  truth <- sim$true
  dat <- sim$wide

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "linear",
    evaluation_indicator = "eval"
  )

  init <- initialize_parameters(dyn$model_system, dat, verbose = FALSE)
  params <- init$init_params
  params["factor_var_1"] <- truth$var_f1
  params["factor_var_2"] <- truth$beta^2 * truth$var_f1 + truth$sigma_e2
  for (i in seq_along(truth$item_int)) {
    params[paste0("Y_t1_m", i, "_intercept")] <- truth$item_int[i]
    params[paste0("Y_t2_m", i, "_intercept")] <- truth$item_int[i] +
      truth$item_load[i] * truth$alpha
    params[paste0("Y_t1_m", i, "_sigma")]     <- truth$item_sigma[i]
    params[paste0("Y_t2_m", i, "_sigma")]     <- truth$item_sigma[i]
    if (i > 1L) {
      params[paste0("Y_t1_m", i, "_loading_1")] <- truth$item_load[i]
      params[paste0("Y_t2_m", i, "_loading_2")] <- truth$item_load[i]
    }
  }

  md <- factorana:::build_parameter_metadata(dyn$model_system)
  cons <- factorana:::setup_parameter_constraints(
    dyn$model_system, params, md, init$factor_variance_fixed, verbose = FALSE
  )
  pfix <- rep(TRUE, length(params))
  pfix[cons$free_idx] <- FALSE

  grad <- check_gradient_accuracy(dyn$model_system, dat, params,
                                   param_fixed = pfix,
                                   tol = 1e-3, verbose = FALSE, n_quad = 8)
  expect_true(grad$pass,
              info = sprintf("Stage 1 linear gradient FD failed (max err: %.2e)",
                             grad$max_error))

  hess <- check_hessian_accuracy(dyn$model_system, dat, params,
                                  param_fixed = pfix,
                                  tol = 5e-3, verbose = FALSE, n_quad = 8)
  expect_true(hess$pass,
              info = sprintf("Stage 1 linear Hessian FD failed (max err: %.2e)",
                             hess$max_error))
})


# =============================================================================
# TEST 7: Stage 1 FD check (oprobit) for the dynamic-measurement wrapper.
#
# As TEST 6 but with ordered-probit indicators and the threshold-based
# equality constraint structure. Exercises the same Hessian-accumulation
# fix on the oprobit likelihood path (now that C++ correctly maps
# _thresh_k names via the num_choices field in rcpp_interface.cpp and
# the Hessian loops iterate over gradparlist in FactorModel.cpp).
# =============================================================================
test_that("dynamic wrapper Stage 1 (oprobit): FD gradient and Hessian match", {
  skip_on_cran()

  set.seed(82); n <- 800
  item_load <- c(1.0, 0.9, 1.1)
  cuts <- c(-1.0, 0.0, 1.0)
  f1 <- rnorm(n, 0, 1); f2 <- 0.4 + 0.6*f1 + rnorm(n, 0, sqrt(0.5))
  gen <- function(f, i) as.integer(findInterval(item_load[i]*f + rnorm(length(f),0,1), cuts) + 1L)
  dat <- data.frame(
    intercept = 1, eval = 1L,
    Y_t1_m1 = gen(f1, 1), Y_t1_m2 = gen(f1, 2), Y_t1_m3 = gen(f1, 3),
    Y_t2_m1 = gen(f2, 1), Y_t2_m2 = gen(f2, 2), Y_t2_m3 = gen(f2, 3)
  )

  dyn <- define_dynamic_measurement(
    data                 = dat,
    items                = c("m1", "m2", "m3"),
    period_prefixes      = c("Y_t1_", "Y_t2_"),
    model_type           = "oprobit",
    n_categories         = 4L,
    evaluation_indicator = "eval"
  )

  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)
  # Use the initialised parameter point (NOT the MLE): at the MLE the
  # gradient is ~0, so FD central-difference noise of order 1e-4 produces
  # a pathological rel_err of ~1 even when analytical and FD agree in
  # absolute terms. Initialised values give gradients of meaningful
  # magnitude for a clean FD comparison.
  init <- initialize_parameters(dyn$model_system, dat, verbose = FALSE)
  params <- init$init_params

  md <- factorana:::build_parameter_metadata(dyn$model_system)
  cons <- factorana:::setup_parameter_constraints(
    dyn$model_system, params, md, init$factor_variance_fixed, verbose = FALSE
  )
  pfix <- rep(TRUE, length(params))
  pfix[cons$free_idx] <- FALSE

  # Oprobit is more numerically noisy than linear; accept a slightly
  # looser tolerance for the gradient (1.5e-3 vs 1e-3 for linear) to
  # tolerate finite-difference noise around the initialisation point
  # without masking real analytical bugs (which produce rel_err of
  # order 1 on the Hessian-accumulation bug that this test guards).
  grad <- check_gradient_accuracy(dyn$model_system, dat, params,
                                   param_fixed = pfix,
                                   tol = 1.5e-3, verbose = FALSE, n_quad = 6)
  expect_true(grad$pass,
              info = sprintf("Stage 1 oprobit gradient FD failed (max err: %.2e)",
                             grad$max_error))

  hess <- check_hessian_accuracy(dyn$model_system, dat, params,
                                  param_fixed = pfix,
                                  tol = 5e-3, verbose = FALSE, n_quad = 6)
  expect_true(hess$pass,
              info = sprintf("Stage 1 oprobit Hessian FD failed (max err: %.2e)",
                             hess$max_error))
})


# =============================================================================
# TEST 8: Oprobit two-stage parameter recovery (wrapper + SE_linear).
#
# End-to-end oprobit workflow: dynamic-measurement wrapper Stage 1 fits
# a 2-factor ordered-probit measurement system with period-specific
# thresh_1 and tied threshold increments; Stage 2 is SE_linear on top.
# The structural parameters (factor_var_1, se_intercept, se_linear_1,
# se_residual_var) are recovered within tolerance. Became feasible once
# the C++ threshold-name lookup and the Hessian-accumulation bugs were
# fixed (v1.1.6 and v1.1.7 respectively).
# =============================================================================
test_that("dynamic wrapper oprobit two-stage: structural parameter recovery", {
  skip_on_cran()

  set.seed(91); n <- 2500
  f1 <- rnorm(n, 0, 1); eps <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.4 + 0.6 * f1 + eps
  item_load <- c(1.0, 0.9, 1.1); cuts <- c(-1.0, 0.0, 1.0)
  gen <- function(f, i) {
    as.integer(findInterval(item_load[i] * f + rnorm(length(f), 0, 1), cuts) + 1L)
  }
  dat <- data.frame(intercept = 1, eval = 1L,
    Y_t1_m1 = gen(f1, 1), Y_t1_m2 = gen(f1, 2), Y_t1_m3 = gen(f1, 3),
    Y_t2_m1 = gen(f2, 1), Y_t2_m2 = gen(f2, 2), Y_t2_m3 = gen(f2, 3))

  dyn <- define_dynamic_measurement(
    data = dat, items = c("m1", "m2", "m3"),
    period_prefixes = c("Y_t1_", "Y_t2_"),
    model_type = "oprobit", n_categories = 4L,
    evaluation_indicator = "eval"
  )
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  s1 <- estimate_model_rcpp(dyn$model_system, dat, control = ctrl,
                             optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
  expect_equal(s1$convergence, 0)

  dummy <- build_dynamic_previous_stage(dyn, s1, dat, anchor_period = 1L)
  fm_s2 <- define_factor_model(n_factors = 2, n_types = 1, factor_structure = "SE_linear")
  ms_s2 <- define_model_system(components = list(), factor = fm_s2, previous_stage = dummy)
  init2 <- initialize_parameters(ms_s2, dat, verbose = FALSE)
  init2$init_params["factor_var_1"]    <- unname(dummy$estimates["factor_var_1"])
  init2$init_params["se_intercept"]    <- 0
  init2$init_params["se_linear_1"]     <- 0.5
  init2$init_params["se_residual_var"] <- 0.5

  r2 <- estimate_model_rcpp(ms_s2, dat, init_params = init2$init_params,
                             control = ctrl, optimizer = "nlminb",
                             parallel = FALSE, verbose = FALSE)
  expect_equal(r2$convergence, 0)

  est <- r2$estimates
  expect_equal(unname(est["factor_var_1"]),    1.0, tolerance = 0.20)
  expect_equal(unname(est["se_intercept"]),    0.4, tolerance = 0.15)
  expect_equal(unname(est["se_linear_1"]),     0.6, tolerance = 0.15)
  expect_equal(unname(est["se_residual_var"]), 0.5, tolerance = 0.20)
  expect_true(all(r2$std_errors[c("factor_var_1", "se_intercept",
                                    "se_linear_1", "se_residual_var")] > 0))
})
