# Tests for equality constraints (measurement invariance)

test_that("equality constraints are stored in model system", {
  set.seed(123)
  n <- 100
  dat <- data.frame(y1 = rnorm(n), y2 = rnorm(n), intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  # Define with equality constraints
  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  # Check that constraints are stored
  expect_equal(length(ms$equality_constraints), 1)
  expect_equal(ms$equality_constraints[[1]], c("m1_sigma", "m2_sigma"))
})

test_that("equality constraints require at least 2 parameters per group", {
  set.seed(124)
  n <- 100
  dat <- data.frame(y1 = rnorm(n), y2 = rnorm(n), intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  # Single parameter constraint should error
  expect_error(
    define_model_system(components = list(mc1, mc2), factor = fm,
                       equality_constraints = list(c("m1_sigma"))),
    regexp = "at least 2"
  )
})

test_that("tied parameters have identical estimates", {
  skip_on_cran()
  set.seed(125)
  n <- 500

  # Generate data with same true sigma for both components
  f <- rnorm(n, 0, 1)
  true_sigma <- 0.6
  y1 <- 1.0 * f + rnorm(n, 0, true_sigma)
  y2 <- 0.8 * f + rnorm(n, 0, true_sigma)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  # Constrain sigmas to be equal
  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)
  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                 optimizer = "nlminb", verbose = FALSE)

  expect_equal(result$convergence, 0,
               info = sprintf("Model with equality constraints did not converge, code: %d", result$convergence))

  # Check that tied parameters are exactly equal
  sigma1 <- unname(result$estimates["m1_sigma"])
  sigma2 <- unname(result$estimates["m2_sigma"])

  expect_equal(sigma1, sigma2, tolerance = 1e-10,
               info = sprintf("Tied sigmas should be exactly equal: %.6f vs %.6f", sigma1, sigma2))

  # Check that estimated sigma is close to true value
  expect_true(abs(sigma1 - true_sigma) < 0.15,
              info = sprintf("Estimated sigma %.3f differs from true %.3f by more than 0.15",
                            sigma1, true_sigma))
})

test_that("gradient is correct with equality constraints", {
  skip_on_cran()
  set.seed(126)
  n <- 300

  f <- rnorm(n, 0, 1)
  y1 <- 1.0 * f + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  # First param (factor_var) is fixed when loading is fixed
  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 1e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass,
              info = sprintf("Gradient check with equality constraints failed, max error: %.2e", result$max_error))
})

test_that("Hessian is correct with equality constraints", {
  skip_on_cran()
  set.seed(127)
  n <- 300

  f <- rnorm(n, 0, 1)
  y1 <- 1.0 * f + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  ms <- define_model_system(components = list(mc1, mc2), factor = fm,
                            equality_constraints = list(c("m1_sigma", "m2_sigma")))

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 1.0, verbose = FALSE, n_quad = 8)
  expect_true(result$pass,
              info = sprintf("Hessian check with equality constraints failed, max error: %.2e", result$max_error))
})


# =============================================================================
# Regression test for the parameter-name-to-index map in
# initialize_factor_model_cpp (rcpp_interface.cpp). The map was historically
# missing typeprob_*_intercept (it stored the wrong name `type_*_intercept`),
# all type_*_loading_*, every factor_mean_<k>_<cov>, and every se_cov_<cov>.
#
# Beyond making equality constraints on those parameter types unrecognised,
# the missing factor_mean and se_cov entries also shifted every subsequent
# component-level idx down by the number of skipped slots, so equality
# constraints on loadings/sigmas/thresholds were silently mapped to factor-
# distribution slots whenever factor_covariates or se_covariates were used.
#
# This test exercises an SE_linear + se_covariates two-stage workflow with
# the dynamic-measurement pattern (loadings and sigmas tied across waves
# via equality_constraints). With the map bug, the wave-2 loading "tied to"
# wave-1 loading would silently land on a factor-mean / se_cov / typeprob
# slot at Stage 1, the constraint would bind a different parameter, and the
# resulting estimate would differ from a no-cov rerun. With the map fixed,
# the tie binds the correct measurement parameters.
# =============================================================================
test_that("equality_constraints + se_covariates does not corrupt component idx mapping", {
  skip_on_cran()

  set.seed(2027)
  n <- 600

  X <- rnorm(n, 0, 1)
  X_dem <- X - mean(X)

  f <- rnorm(n, 0, 1)
  eps <- rnorm(n, 0, 0.5)
  f2 <- 0.5 * f + 0.4 * X_dem + eps

  item_load <- c(1.0, 0.9, 1.1)
  item_sigma <- c(0.6, 0.65, 0.55)
  gen <- function(f, i) item_load[i] * f + rnorm(length(f), 0, item_sigma[i])
  dat <- data.frame(
    intercept = 1, X = X, eval = 1L,
    Y_t1_m1 = gen(f,  1), Y_t1_m2 = gen(f,  2), Y_t1_m3 = gen(f,  3),
    Y_t2_m1 = gen(f2, 1), Y_t2_m2 = gen(f2, 2), Y_t2_m3 = gen(f2, 3)
  )

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear",
                            se_covariates = c("X"))

  mk <- function(name, outcome, norm_idx) {
    norm <- c(0, 0); norm[norm_idx] <- if (grepl("_m1$", name)) 1 else NA_real_
    mc <- define_model_component(
      name = name, data = dat, outcome = outcome, factor = fm,
      covariates = "intercept", model_type = "linear",
      loading_normalization = norm,
      evaluation_indicator = "eval"
    )
    # Fix wave-2 (outcome-factor) measurement intercepts to 0 so se_intercept
    # is properly identified — otherwise se_intercept walks along a flat
    # ridge with the wave-2 intercepts. Same trick as test-se-models.R:574.
    if (norm_idx == 2L) mc <- fix_coefficient(mc, "intercept", 0)
    mc
  }
  comps <- list(
    mk("Y_t1_m1", "Y_t1_m1", 1),
    mk("Y_t1_m2", "Y_t1_m2", 1),
    mk("Y_t1_m3", "Y_t1_m3", 1),
    mk("Y_t2_m1", "Y_t2_m1", 2),
    mk("Y_t2_m2", "Y_t2_m2", 2),
    mk("Y_t2_m3", "Y_t2_m3", 2)
  )

  # Tie wave-1 / wave-2 loadings and sigmas. With the map bug + se_covariates,
  # the component-level lookup would point at wrong indices.
  eq <- list(
    c("Y_t1_m2_loading_1", "Y_t2_m2_loading_2"),
    c("Y_t1_m3_loading_1", "Y_t2_m3_loading_2"),
    c("Y_t1_m1_sigma",     "Y_t2_m1_sigma"),
    c("Y_t1_m2_sigma",     "Y_t2_m2_sigma"),
    c("Y_t1_m3_sigma",     "Y_t2_m3_sigma")
  )

  ms <- define_model_system(components = comps, factor = fm,
                            equality_constraints = eq)

  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result <- estimate_model_rcpp(ms, dat, control = ctrl,
                                 optimizer = "nlminb", parallel = FALSE,
                                 verbose = FALSE)
  expect_equal(result$convergence, 0,
               info = "SE_linear + se_covariates + equality_constraints must converge")

  est <- result$estimates

  # The equality constraints must bind exactly: tied params share their
  # primary's value (no slack from a misaligned lookup).
  expect_equal(unname(est["Y_t1_m2_loading_1"]),
               unname(est["Y_t2_m2_loading_2"]), tolerance = 1e-10,
               info = "loading equality (m2) should bind tied -> primary exactly")
  expect_equal(unname(est["Y_t1_m3_loading_1"]),
               unname(est["Y_t2_m3_loading_2"]), tolerance = 1e-10,
               info = "loading equality (m3) should bind tied -> primary exactly")
  expect_equal(unname(est["Y_t1_m1_sigma"]),
               unname(est["Y_t2_m1_sigma"]),     tolerance = 1e-10,
               info = "sigma equality (m1) should bind tied -> primary exactly")
  expect_equal(unname(est["Y_t1_m2_sigma"]),
               unname(est["Y_t2_m2_sigma"]),     tolerance = 1e-10,
               info = "sigma equality (m2) should bind tied -> primary exactly")
  expect_equal(unname(est["Y_t1_m3_sigma"]),
               unname(est["Y_t2_m3_sigma"]),     tolerance = 1e-10,
               info = "sigma equality (m3) should bind tied -> primary exactly")

  # se_cov_X must recover its DGP value (within sampling noise + finite-n
  # bias). The test guards that the se_cov_X slot holds the SE-covariate
  # coefficient, not some shifted slot's value.
  expect_lt(abs(unname(est["se_cov_X"]) - 0.4), 0.15)

  # Wave-1 measurement loadings and sigmas should recover roughly to their
  # DGP values; this catches gross slot-misalignment that would leave them
  # pinned to factor-distribution slots.
  expect_lt(abs(unname(est["Y_t1_m2_loading_1"]) - 0.9),  0.15)
  expect_lt(abs(unname(est["Y_t1_m3_loading_1"]) - 1.1),  0.15)
  expect_lt(abs(unname(est["Y_t1_m1_sigma"])     - 0.6),  0.15)
  expect_lt(abs(unname(est["Y_t1_m2_sigma"])     - 0.65), 0.15)
  expect_lt(abs(unname(est["Y_t1_m3_sigma"])     - 0.55), 0.15)
})
