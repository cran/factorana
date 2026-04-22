# Tests for structural equation models (SE_linear, SE_quadratic)

test_that("SE_linear gradient is accurate vs finite differences", {
  skip_on_cran()
  set.seed(123)
  n <- 300

  # Generate data with SE_linear structure: f2 = alpha + alpha1*f1 + epsilon
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))  # Residual variance = 0.5
  f2 <- 0.2 + 0.8 * f1 + epsilon  # SE structure

  # Measurement equations
  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)  # f1 measure 1 (loading = 1, fixed)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)  # f1 measure 2
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)  # f2 measure 1 (loading = 1, fixed)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)  # f2 measure 2

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  # Define SE_linear model
  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  # Determine which parameters are fixed
  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE  # factor_var_1 is fixed when loading is fixed

  # Check gradient accuracy - use looser tolerance for SE models due to numerical precision
  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_linear gradient check failed, max error: %.2e", result$max_error))
})

test_that("SE_linear Hessian is accurate vs finite differences", {
  skip_on_cran()
  set.seed(124)
  n <- 300

  # Generate data with SE_linear structure
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.2 + 0.8 * f1 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  # Check Hessian accuracy - use looser tolerance for SE models
  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_linear Hessian check failed, max error: %.2e", result$max_error))
})

test_that("SE_quadratic gradient is accurate vs finite differences", {
  skip_on_cran()
  set.seed(125)
  n <- 300

  # Generate data with SE_quadratic structure: f2 = alpha + alpha1*f1 + alpha2*f1^2 + epsilon
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.2 + 0.6 * f1 + 0.1 * f1^2 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_gradient_accuracy(ms, dat, params, param_fixed = param_fixed,
                                     tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_quadratic gradient check failed, max error: %.2e", result$max_error))
})

test_that("SE_quadratic Hessian is accurate vs finite differences", {
  skip_on_cran()
  set.seed(126)
  n <- 300

  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.2 + 0.6 * f1 + 0.1 * f1^2 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_quadratic")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  mc3 <- define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1))
  mc4 <- define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA))

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params

  param_fixed <- rep(FALSE, length(params))
  param_fixed[1] <- TRUE

  result <- check_hessian_accuracy(ms, dat, params, param_fixed = param_fixed,
                                    tol = 5e-2, verbose = FALSE, n_quad = 8)
  expect_true(result$pass, info = sprintf("SE_quadratic Hessian check failed, max error: %.2e", result$max_error))
})

test_that("SE_linear model converges", {
  skip_on_cran()
  set.seed(128)
  n <- 500

  # Generate data with SE_linear structure
  f1 <- rnorm(n, 0, 1)
  epsilon <- rnorm(n, 0, sqrt(0.5))
  f2 <- 0.7 * f1 + epsilon

  y1 <- 1.0 * f1 + rnorm(n, 0, 0.5)
  y2 <- 0.8 * f1 + rnorm(n, 0, 0.5)
  y3 <- 1.0 * f2 + rnorm(n, 0, 0.5)
  y4 <- 0.9 * f2 + rnorm(n, 0, 0.5)

  dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4, intercept = 1)

  fm <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")

  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(1, 0))
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(NA, 0))
  # Fix the outcome-factor measurement intercepts to 0 for identification.
  # Without this, `se_intercept` and the m3/m4 intercepts form a flat ridge
  # (any shift in se_intercept can be absorbed by shifts in the outcome-
  # measurement intercepts), and nlminb walks along the ridge indefinitely.
  mc3 <- fix_coefficient(define_model_component("m3", dat, "y3", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, 1)), "intercept", 0)
  mc4 <- fix_coefficient(define_model_component("m4", dat, "y4", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = c(0, NA)), "intercept", 0)

  ms <- define_model_system(components = list(mc1, mc2, mc3, mc4), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                 optimizer = "nlminb", verbose = FALSE)

  # Check strict convergence (only code 0 is accepted; see CLAUDE.md).
  expect_equal(result$convergence, 0,
               info = sprintf("SE_linear model did not converge, code: %d",
                              result$convergence))

  # Check SE parameters exist and are finite
  expect_true(is.finite(result$estimates["se_intercept"]))
  expect_true(is.finite(result$estimates["se_linear_1"]))
  expect_true(is.finite(result$estimates["se_residual_var"]))
  expect_true(result$estimates["se_residual_var"] > 0,
              info = "SE residual variance should be positive")
})

# =============================================================================
# Type-specific SE intercept tests
# =============================================================================
#
# These tests exercise the se_intercept_type_{t} parameters introduced to allow
# the outcome-factor intercept in SE_linear / SE_quadratic models to vary by
# latent type. The feature is needed for regime-switching structural models
# (e.g., MH trap project) where the slopes are shared across types but the
# intercept on the outcome factor differs.
#
# Semantic: type probabilities are a function of the INPUT factors only
# (type is drawn before the residual epsilon). Type loadings on the outcome
# factor are fixed to zero at the R level.

.se_type_param_fixed <- function(model_system, data, params) {
  # Compute param_fixed using the same logic the optimizer uses, so the FD
  # helpers know which parameters are internally fixed (e.g., auto-fixed type
  # loadings on the outcome factor for SE models).
  init <- initialize_parameters(model_system, data, verbose = FALSE)
  pm <- factorana:::build_parameter_metadata(model_system)
  pc <- factorana:::setup_parameter_constraints(
    model_system, params, pm, init$factor_variance_fixed, verbose = FALSE)
  pc$param_fixed
}

.build_se_type_model <- function(n_factors = 2, n_types = 2,
                                 factor_structure = "SE_linear",
                                 n_meas_per_factor = 3,
                                 n = 150, seed = 1,
                                 indicator_type = "linear",
                                 n_categories    = 4L,
                                 oprobit_cuts    = c(-1.0, 0.0, 1.0)) {
  stopifnot(indicator_type %in% c("linear", "oprobit"))
  if (indicator_type == "oprobit") {
    stopifnot(length(oprobit_cuts) == n_categories - 1L)
  }

  set.seed(seed)
  # Draw input factors and residual
  n_input <- n_factors - 1L
  f_input <- matrix(rnorm(n * n_input), nrow = n, ncol = n_input)

  # Simple type assignment (deterministic mix of halves) — does not depend on
  # the outcome factor. For the FD tests we do not need realistic type draws.
  type_id <- sample(seq_len(n_types), n, replace = TRUE)

  # Truth for SE equation
  se_int <- 0.1
  se_lin <- rep(0.8, n_input)
  se_quad <- rep(0.2, n_input)
  se_type_offset <- c(0.0, seq_len(n_types - 1L) * 0.4)  # 0, 0.4, 0.8, ...
  se_res_sd <- sqrt(0.5)

  eps <- rnorm(n, 0, se_res_sd)
  if (factor_structure == "SE_quadratic") {
    f_out <- se_int + as.numeric(f_input %*% se_lin) +
      as.numeric((f_input^2) %*% se_quad) + se_type_offset[type_id] + eps
  } else {
    f_out <- se_int + as.numeric(f_input %*% se_lin) + se_type_offset[type_id] + eps
  }

  # Generate measurements: n_meas_per_factor per factor
  dat <- data.frame(intercept = rep(1, n))
  mc_list <- list()
  fm <- define_factor_model(n_factors = n_factors, n_types = n_types,
                            factor_structure = factor_structure)

  for (k in seq_len(n_factors)) {
    fk <- if (k <= n_input) f_input[, k] else f_out
    for (jj in seq_len(n_meas_per_factor)) {
      load_true <- if (jj == 1) 1.0 else 0.8 + 0.1 * jj
      col <- sprintf("f%d_m%d", k, jj)

      if (indicator_type == "linear") {
        dat[[col]] <- load_true * fk + rnorm(n, 0, 0.4)
        loading_norm <- rep(0, n_factors)
        if (jj == 1) loading_norm[k] <- 1 else loading_norm[k] <- NA_real_
        mc_list[[length(mc_list) + 1L]] <- define_model_component(
          col, dat, col, fm,
          covariates = "intercept", model_type = "linear",
          loading_normalization = loading_norm
        )
      } else {
        # Oprobit: discretize the latent Y* into n_categories
        ystar <- load_true * fk + rnorm(n, 0, 1)  # probit scale (sigma = 1)
        dat[[col]] <- as.integer(findInterval(ystar, oprobit_cuts) + 1L)
        loading_norm <- rep(0, n_factors)
        if (jj == 1) loading_norm[k] <- 1 else loading_norm[k] <- NA_real_
        mc_list[[length(mc_list) + 1L]] <- define_model_component(
          col, dat, col, fm,
          covariates = NULL, model_type = "oprobit",
          num_choices = n_categories,
          loading_normalization = loading_norm
        )
      }
    }
  }

  ms <- define_model_system(components = mc_list, factor = fm)
  list(model_system = ms, data = dat, truth = list(
    se_intercept = se_int, se_linear = se_lin, se_quadratic = se_quad,
    se_type_offset = se_type_offset, se_residual_var = se_res_sd^2
  ))
}

test_that("SE_linear with n_types=2 gradient matches finite differences", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_linear",
                              n_meas_per_factor = 3, n = 120, seed = 11)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params

  # Nudge the type intercept and slope to non-zero values so the gradient
  # is tested at a point where the new parameter is active.
  params["se_intercept_type_2"] <- 0.3
  params["se_intercept"] <- 0.1
  params["se_linear_1"] <- 0.7
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  res <- check_gradient_accuracy(mod$model_system, mod$data, params,
                                 param_fixed = pfix,
                                 tol = 1e-4, verbose = FALSE, n_quad = 8)
  expect_true(res$pass,
              info = sprintf("SE_linear ntyp=2 gradient check failed, max error: %.2e",
                             res$max_error))
})

test_that("SE_linear with n_types=2 Hessian matches finite differences", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_linear",
                              n_meas_per_factor = 3, n = 120, seed = 12)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params

  params["se_intercept_type_2"] <- 0.3
  params["se_intercept"] <- 0.1
  params["se_linear_1"] <- 0.7
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  res <- check_hessian_accuracy(mod$model_system, mod$data, params,
                                param_fixed = pfix,
                                tol = 1e-3, verbose = FALSE, n_quad = 8)
  expect_true(res$pass,
              info = sprintf("SE_linear ntyp=2 Hessian check failed, max error: %.2e",
                             res$max_error))
})

# ---- Oprobit indicator variant (no equality constraints) --------------------
# Exercises SE_linear + n_types=2 with ORDERED PROBIT indicators on every
# measurement component. Covers the same likelihood path as the linear
# version above but with threshold-based measurements, verifying the
# analytical gradient and Hessian against finite differences.
test_that("SE_linear with n_types=2 and oprobit indicators: gradient matches FD", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_linear",
                              n_meas_per_factor = 3, n = 150, seed = 13,
                              indicator_type = "oprobit",
                              n_categories = 4L)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params
  params["se_intercept_type_2"] <- 0.3
  params["se_intercept"]        <- 0.1
  params["se_linear_1"]         <- 0.7
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  res <- check_gradient_accuracy(mod$model_system, mod$data, params,
                                  param_fixed = pfix,
                                  tol = 1e-4, verbose = FALSE, n_quad = 8)
  expect_true(res$pass,
              info = sprintf("SE_linear ntyp=2 + oprobit gradient FD failed, max err: %.2e",
                             res$max_error))
})

test_that("SE_linear with n_types=2 and oprobit indicators: Hessian matches FD", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_linear",
                              n_meas_per_factor = 3, n = 150, seed = 14,
                              indicator_type = "oprobit",
                              n_categories = 4L)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params
  params["se_intercept_type_2"] <- 0.3
  params["se_intercept"]        <- 0.1
  params["se_linear_1"]         <- 0.7
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  res <- check_hessian_accuracy(mod$model_system, mod$data, params,
                                 param_fixed = pfix,
                                 tol = 1e-3, verbose = FALSE, n_quad = 8)
  expect_true(res$pass,
              info = sprintf("SE_linear ntyp=2 + oprobit Hessian FD failed, max err: %.2e",
                             res$max_error))
})

test_that("SE_quadratic with n_types=2 gradient and Hessian match finite differences", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 2, n_types = 2,
                              factor_structure = "SE_quadratic",
                              n_meas_per_factor = 3, n = 120, seed = 21)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params

  params["se_intercept_type_2"] <- 0.25
  params["se_intercept"] <- 0.1
  params["se_linear_1"] <- 0.7
  params["se_quadratic_1"] <- 0.15
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  g <- check_gradient_accuracy(mod$model_system, mod$data, params,
                               param_fixed = pfix,
                               tol = 1e-4, verbose = FALSE, n_quad = 8)
  expect_true(g$pass,
              info = sprintf("SE_quadratic ntyp=2 gradient check failed, max error: %.2e",
                             g$max_error))

  h <- check_hessian_accuracy(mod$model_system, mod$data, params,
                              param_fixed = pfix,
                              tol = 1e-3, verbose = FALSE, n_quad = 8)
  expect_true(h$pass,
              info = sprintf("SE_quadratic ntyp=2 Hessian check failed, max error: %.2e",
                             h$max_error))
})

test_that("SE_linear with n_types=3 gradient and Hessian match finite differences", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 2, n_types = 3,
                              factor_structure = "SE_linear",
                              n_meas_per_factor = 3, n = 120, seed = 31)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params

  params["se_intercept_type_2"] <- 0.25
  params["se_intercept_type_3"] <- -0.35
  params["se_intercept"] <- 0.1
  params["se_linear_1"] <- 0.7
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  # Slightly looser tolerance than ntyp=2 tests: more integration points + more
  # type evaluations amplify numerical noise in the gradient / FD comparison.
  g <- check_gradient_accuracy(mod$model_system, mod$data, params,
                               param_fixed = pfix,
                               tol = 5e-4, verbose = FALSE, n_quad = 10)
  expect_true(g$pass,
              info = sprintf("SE_linear ntyp=3 gradient check failed, max error: %.2e",
                             g$max_error))

  h <- check_hessian_accuracy(mod$model_system, mod$data, params,
                              param_fixed = pfix,
                              tol = 5e-3, verbose = FALSE, n_quad = 10)
  expect_true(h$pass,
              info = sprintf("SE_linear ntyp=3 Hessian check failed, max error: %.2e",
                             h$max_error))
})

test_that("SE_linear with 3 factors and n_types=2 gradient and Hessian match finite differences", {
  skip_on_cran()
  mod <- .build_se_type_model(n_factors = 3, n_types = 2,
                              factor_structure = "SE_linear",
                              n_meas_per_factor = 3, n = 150, seed = 41)
  init <- initialize_parameters(mod$model_system, mod$data, verbose = FALSE)
  params <- init$init_params

  params["se_intercept_type_2"] <- 0.3
  params["se_intercept"] <- 0.1
  params["se_linear_1"] <- 0.6
  params["se_linear_2"] <- 0.5
  pfix <- .se_type_param_fixed(mod$model_system, mod$data, params)

  # 3 factors + 2 types → 6^3 = 216 integration points; noise amplifies in FD.
  g <- check_gradient_accuracy(mod$model_system, mod$data, params,
                               param_fixed = pfix,
                               tol = 5e-4, verbose = FALSE, n_quad = 8)
  expect_true(g$pass,
              info = sprintf("SE_linear 3-factor ntyp=2 gradient check failed, max error: %.2e",
                             g$max_error))

  h <- check_hessian_accuracy(mod$model_system, mod$data, params,
                              param_fixed = pfix,
                              tol = 3e-2, verbose = FALSE, n_quad = 8)
  expect_true(h$pass,
              info = sprintf("SE_linear 3-factor ntyp=2 Hessian check failed, max error: %.2e",
                             h$max_error))
})

test_that("SE_linear recovers type-specific intercept from simulated data", {
  skip_on_cran()

  set.seed(2024)
  n <- 1500
  n_meas <- 3

  # Truth
  se_intercept_true <- 0.1
  se_linear_true <- 0.7
  se_intercept_type_2_true <- 0.8  # Large enough to recover with n=1500
  se_residual_sd_true <- sqrt(0.4)

  # Draw input factor and type probabilities (depend on f1 only)
  f1 <- rnorm(n, 0, 1)
  log_odds_type2 <- -0.5 + 1.0 * f1  # typeprob intercept and loading on f1
  p_type2 <- plogis(log_odds_type2)
  type_id <- 1L + (runif(n) < p_type2)  # 1 or 2

  eps <- rnorm(n, 0, se_residual_sd_true)
  se_offset <- c(0, se_intercept_type_2_true)
  f2 <- se_intercept_true + se_linear_true * f1 + se_offset[type_id] + eps

  # Measurements (loadings: 1.0 for first meas (fixed), 0.9 / 0.8 for others)
  loads1 <- c(1.0, 0.9, 0.8)
  loads2 <- c(1.0, 0.9, 0.8)
  dat <- data.frame(intercept = rep(1, n))
  for (j in 1:n_meas) {
    dat[[sprintf("f1_m%d", j)]] <- loads1[j] * f1 + rnorm(n, 0, 0.4)
    dat[[sprintf("f2_m%d", j)]] <- loads2[j] * f2 + rnorm(n, 0, 0.4)
  }

  fm <- define_factor_model(n_factors = 2, n_types = 2, factor_structure = "SE_linear")

  comps <- list()
  for (j in 1:n_meas) {
    ln <- if (j == 1) c(1, 0) else c(NA_real_, 0)
    comps[[length(comps) + 1L]] <- define_model_component(
      sprintf("f1_m%d", j), dat, sprintf("f1_m%d", j), fm,
      covariates = "intercept", model_type = "linear", loading_normalization = ln
    )
  }
  # Outcome-factor (f2) measurements: fix the intercept to 0 so that se_intercept
  # is properly identified. Otherwise se_intercept trades off 1-for-1 with the
  # f2-measurement intercepts and nlminb walks along a flat ridge indefinitely.
  for (j in 1:n_meas) {
    ln <- if (j == 1) c(0, 1) else c(0, NA_real_)
    mc <- define_model_component(
      sprintf("f2_m%d", j), dat, sprintf("f2_m%d", j), fm,
      covariates = "intercept", model_type = "linear", loading_normalization = ln
    )
    comps[[length(comps) + 1L]] <- fix_coefficient(mc, "intercept", 0)
  }

  ms <- define_model_system(components = comps, factor = fm)
  # n_quad = 16 — lower resolution (e.g. 8) can trap the optimizer in a
  # distant local mode for this multi-type SE model.
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                optimizer = "nlminb", verbose = FALSE)

  # Strict convergence: only code 0 is accepted (see CLAUDE.md).
  expect_equal(result$convergence, 0,
               info = sprintf("SE_linear ntyp=2 recovery did not converge (code %d)",
                              result$convergence))

  # Parameter checks.
  # NOTE: Type models are identified only up to a permutation of type labels
  # (label switching). Swapping types 1 and 2 negates se_intercept_type_2 and
  # shifts se_intercept by the same amount. So we check the MAGNITUDE of
  # se_intercept_type_2 (which is the "gap" between the two types' intercepts)
  # rather than its signed value.
  est_type2 <- result$estimates["se_intercept_type_2"]
  se_type2 <- result$std_errors["se_intercept_type_2"]

  expect_true(is.finite(est_type2))
  expect_true(is.finite(se_type2) && se_type2 < 1.0,
              info = sprintf("SE for se_intercept_type_2 unreasonable: %.3f", se_type2))

  # Magnitude within 3 SEs of the true magnitude.
  expect_lt(abs(abs(est_type2) - se_intercept_type_2_true), 3 * se_type2)

  # Cross-check: se_linear_1 recovered (invariant to label swap).
  est_linear <- result$estimates["se_linear_1"]
  se_linear <- result$std_errors["se_linear_1"]
  expect_lt(abs(est_linear - se_linear_true), 3 * se_linear)
})

test_that("SE_quadratic recovers type-specific intercept from simulated data", {
  skip_on_cran()

  set.seed(2025)
  n <- 1500
  n_meas <- 3

  se_intercept_true <- 0.1
  se_linear_true <- 0.6
  se_quadratic_true <- 0.25
  se_intercept_type_2_true <- 0.7
  se_residual_sd_true <- sqrt(0.4)

  f1 <- rnorm(n, 0, 1)
  log_odds_type2 <- -0.4 + 0.8 * f1
  p_type2 <- plogis(log_odds_type2)
  type_id <- 1L + (runif(n) < p_type2)

  eps <- rnorm(n, 0, se_residual_sd_true)
  se_offset <- c(0, se_intercept_type_2_true)
  f2 <- se_intercept_true + se_linear_true * f1 + se_quadratic_true * f1^2 +
    se_offset[type_id] + eps

  loads1 <- c(1.0, 0.9, 0.8)
  loads2 <- c(1.0, 0.9, 0.8)
  dat <- data.frame(intercept = rep(1, n))
  for (j in 1:n_meas) {
    dat[[sprintf("f1_m%d", j)]] <- loads1[j] * f1 + rnorm(n, 0, 0.4)
    dat[[sprintf("f2_m%d", j)]] <- loads2[j] * f2 + rnorm(n, 0, 0.4)
  }

  fm <- define_factor_model(n_factors = 2, n_types = 2, factor_structure = "SE_quadratic")
  comps <- list()
  for (j in 1:n_meas) {
    ln <- if (j == 1) c(1, 0) else c(NA_real_, 0)
    comps[[length(comps) + 1L]] <- define_model_component(
      sprintf("f1_m%d", j), dat, sprintf("f1_m%d", j), fm,
      covariates = "intercept", model_type = "linear", loading_normalization = ln
    )
  }
  # Fix outcome-factor measurement intercepts for identification (see SE_linear
  # recovery test above for explanation).
  for (j in 1:n_meas) {
    ln <- if (j == 1) c(0, 1) else c(0, NA_real_)
    mc <- define_model_component(
      sprintf("f2_m%d", j), dat, sprintf("f2_m%d", j), fm,
      covariates = "intercept", model_type = "linear", loading_normalization = ln
    )
    comps[[length(comps) + 1L]] <- fix_coefficient(mc, "intercept", 0)
  }

  ms <- define_model_system(components = comps, factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 16, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                optimizer = "nlminb", verbose = FALSE)

  # Strict convergence: only code 0 is accepted (see CLAUDE.md).
  expect_equal(result$convergence, 0,
               info = sprintf("SE_quadratic ntyp=2 recovery did not converge (code %d)",
                              result$convergence))

  # Check magnitude to be invariant to type-label switching (see SE_linear test).
  est_type2 <- result$estimates["se_intercept_type_2"]
  se_type2 <- result$std_errors["se_intercept_type_2"]
  expect_true(is.finite(est_type2))
  expect_true(is.finite(se_type2) && se_type2 < 1.0,
              info = sprintf("SE for se_intercept_type_2 unreasonable: %.3f", se_type2))
  expect_lt(abs(abs(est_type2) - se_intercept_type_2_true), 3 * se_type2)
})

test_that("SE models error when user tries to set a non-zero type loading on the outcome factor", {
  skip_on_cran()
  set.seed(99)
  n <- 60
  dat <- data.frame(intercept = 1, y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n),
                    y4 = rnorm(n), y5 = rnorm(n), y6 = rnorm(n))
  fm <- define_factor_model(n_factors = 2, n_types = 2, factor_structure = "SE_linear")
  mcs <- list(
    define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(1, 0)),
    define_model_component("m2", dat, "y2", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(NA_real_, 0)),
    define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(NA_real_, 0)),
    define_model_component("o1", dat, "y4", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, 1)),
    define_model_component("o2", dat, "y5", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, NA_real_)),
    define_model_component("o3", dat, "y6", fm, covariates = "intercept",
                           model_type = "linear", loading_normalization = c(0, NA_real_))
  )
  ms <- define_model_system(components = mcs, factor = fm)

  init <- initialize_parameters(ms, dat, verbose = FALSE)
  params <- init$init_params
  # Explicitly override the type-loading on the outcome factor (index 2) to non-zero.
  params["type_2_loading_2"] <- 0.5
  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)

  expect_error(
    estimate_model_rcpp(ms, dat, init_params = params, control = ctrl,
                        parallel = FALSE, optimizer = "nlminb", verbose = FALSE),
    "outcome factor is not allowed"
  )
})
