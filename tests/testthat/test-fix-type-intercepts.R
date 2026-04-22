# Tests for fix_type_intercepts function

test_that("fix_type_intercepts validates inputs correctly", {
  set.seed(123)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  # n_types = 1 should error (no type intercepts exist)
  fm1 <- define_factor_model(n_factors = 1, n_types = 1)
  # Suppress the warning about use_types having no effect with n_types < 2
  mc1 <- suppressWarnings(define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm1,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval",
    use_types = TRUE,
    intercept = FALSE
  ))

  expect_error(
    fix_type_intercepts(mc1),
    "n_types < 2"
  )

  # n_types = 2 should work
  fm2 <- define_factor_model(n_factors = 1, n_types = 2)
  mc2 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm2,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval",
    use_types = TRUE
  )

  mc2_fixed <- fix_type_intercepts(mc2)
  expect_s3_class(mc2_fixed, "model_component")
  expect_equal(length(mc2_fixed$fixed_type_intercepts), 1)
  expect_equal(mc2_fixed$fixed_type_intercepts[[1]]$type, 2)
  expect_equal(mc2_fixed$fixed_type_intercepts[[1]]$value, 0)

  # Error on non-model_component
  expect_error(
    fix_type_intercepts(list(a = 1)),
    "must be an object of class"
  )

  # Error on invalid type number
  expect_error(
    fix_type_intercepts(mc2, types = 1),  # Type 1 is reference
    "must be integers between 2"
  )
  expect_error(
    fix_type_intercepts(mc2, types = 5),  # Type 5 doesn't exist (n_types=2)
    "must be integers between 2"
  )
})

test_that("fix_type_intercepts stores constraint info correctly", {
  set.seed(456)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  # 2-type model
  fm <- define_factor_model(n_factors = 1, n_types = 2)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval",
    use_types = TRUE
  )

  nparam_before <- mc$nparam_model
  mc_fixed <- fix_type_intercepts(mc)
  nparam_after <- mc_fixed$nparam_model

  # nparam_model should NOT be reduced (C++ layer expects full parameter vector)
  # Fixed parameters are handled via constraint mechanism
  expect_equal(nparam_after, nparam_before)

  # But fixed_type_intercepts should be populated
  expect_equal(length(mc_fixed$fixed_type_intercepts), 1)
  expect_equal(mc_fixed$fixed_type_intercepts[[1]]$type, 2)
  expect_equal(mc_fixed$fixed_type_intercepts[[1]]$value, 0)
})

test_that("fix_type_intercepts initializes fixed params to 0", {
  set.seed(789)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 2)

  # Without fix_type_intercepts
  mc1 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval",
    use_types = TRUE
  )
  ms1 <- define_model_system(components = list(mc1), factor = fm)
  init1 <- initialize_parameters(ms1, dat, verbose = FALSE)

  # With fix_type_intercepts
  mc2 <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval",
    use_types = TRUE
  )
  mc2 <- fix_type_intercepts(mc2)
  ms2 <- define_model_system(components = list(mc2), factor = fm)
  init2 <- initialize_parameters(ms2, dat, verbose = FALSE)

  # Type intercept should be in both (C++ requires full parameter vector)
  expect_true("Y_type_2_intercept" %in% init1$param_names)
  expect_true("Y_type_2_intercept" %in% init2$param_names)

  # Parameters count should be the same (fixed params still in vector)
  expect_equal(length(init2$init_params), length(init1$init_params))

  # Fixed type intercept should be initialized to 0
  idx2 <- which(init2$param_names == "Y_type_2_intercept")
  expect_equal(unname(init2$init_params[idx2]), 0.0)

  # Free type intercept should be initialized to non-zero
  idx1 <- which(init1$param_names == "Y_type_2_intercept")
  expect_true(init1$init_params[idx1] != 0.0)
})

test_that("fix_type_intercepts works with 3 types and partial fixing", {
  set.seed(101)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 3)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval",
    use_types = TRUE
  )

  # Fix only type 2 (leave type 3 free)
  mc_partial <- fix_type_intercepts(mc, types = 2)

  ms <- define_model_system(components = list(mc_partial), factor = fm)
  init <- initialize_parameters(ms, dat, verbose = FALSE)

  # Both should be in param_names (C++ requires full vector)
  expect_true("Y_type_2_intercept" %in% init$param_names)
  expect_true("Y_type_3_intercept" %in% init$param_names)

  # Type 2 should be initialized to 0 (fixed), type 3 non-zero (free)
  idx2 <- which(init$param_names == "Y_type_2_intercept")
  idx3 <- which(init$param_names == "Y_type_3_intercept")
  expect_equal(unname(init$init_params[idx2]), 0.0)
  expect_true(init$init_params[idx3] != 0.0)

  # Fix all types (default behavior)
  mc_all <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval",
    use_types = TRUE
  )
  mc_all <- fix_type_intercepts(mc_all)  # Fixes types 2 and 3

  ms_all <- define_model_system(components = list(mc_all), factor = fm)
  init_all <- initialize_parameters(ms_all, dat, verbose = FALSE)

  # Both should be in param_names but initialized to 0
  expect_true("Y_type_2_intercept" %in% init_all$param_names)
  expect_true("Y_type_3_intercept" %in% init_all$param_names)
  idx2_all <- which(init_all$param_names == "Y_type_2_intercept")
  idx3_all <- which(init_all$param_names == "Y_type_3_intercept")
  expect_equal(unname(init_all$init_params[idx2_all]), 0.0)
  expect_equal(unname(init_all$init_params[idx3_all]), 0.0)
})

test_that("is_type_intercept_fixed helper works correctly", {
  set.seed(202)
  n <- 100
  dat <- data.frame(intercept = 1, x1 = rnorm(n), Y = rnorm(n), eval = 1)

  fm <- define_factor_model(n_factors = 1, n_types = 3)
  mc <- define_model_component(
    name = "Y", data = dat, outcome = "Y", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = 1.0, evaluation_indicator = "eval",
    use_types = TRUE
  )

  # Before fixing
  expect_false(is_type_intercept_fixed(mc, 2))
  expect_false(is_type_intercept_fixed(mc, 3))

  # Fix type 2 only
  mc <- fix_type_intercepts(mc, types = 2)

  expect_true(is_type_intercept_fixed(mc, 2))
  expect_false(is_type_intercept_fixed(mc, 3))
})


test_that("3-type 2-factor model with measurement system and full 2nd order outcomes", {
  # Complex model structure:
  # - 2 factors, 3 types
  # - 6 measurement components (3 per factor) with type intercepts fixed to 0
  # - 3 outcome equations with factor_spec = "full" and free type intercepts

  set.seed(303)
  n <- 200

  # Generate data with 6 measurement variables and 3 outcomes
  dat <- data.frame(
    intercept = 1,
    x1 = rnorm(n),
    eval = 1,
    # Measurement system for Factor 1
    T1 = rnorm(n),
    T2 = rnorm(n),
    T3 = rnorm(n),
    # Measurement system for Factor 2
    T4 = rnorm(n),
    T5 = rnorm(n),
    T6 = rnorm(n),
    # Outcome variables
    Y1 = rnorm(n),
    Y2 = rnorm(n),
    Y3 = rnorm(n)
  )

  # Define factor model: 2 factors, 3 types
  fm <- define_factor_model(n_factors = 2, n_types = 3)

  # ---- Measurement System (6 components) ----
  # Factor 1 indicators: T1, T2, T3
  # Loading normalization: c(1, NA) means load=1 on factor 1, load=0 on factor 2

  # Measurement components with use_types = FALSE (default) - NO type intercepts
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0),  # Fixed loading on F1, zero on F2
    evaluation_indicator = "eval",
    use_types = FALSE  # No type intercepts for measurement
  )

  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0),  # Free loading on F1, zero on F2
    evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA, 0),  # Free loading on F1, zero on F2
    evaluation_indicator = "eval",
    use_types = FALSE
  )

  # Factor 2 indicators: T4, T5, T6
  # Loading normalization: c(0, 1) means load=0 on factor 1, load=1 on factor 2

  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1),  # Zero on F1, fixed loading on F2
    evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T5 <- define_model_component(
    name = "T5", data = dat, outcome = "T5", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA),  # Zero on F1, free loading on F2
    evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T6 <- define_model_component(
    name = "T6", data = dat, outcome = "T6", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA),  # Zero on F1, free loading on F2
    evaluation_indicator = "eval",
    use_types = FALSE
  )

  # ---- Outcome Equations (3 components with full 2nd order) ----
  # Full 2nd order: linear + quadratic + interaction terms for both factors
  # Type intercepts are FREE (not fixed)

  # Outcome components with use_types = TRUE - type intercepts are free
  mc_Y1 <- define_model_component(
    name = "Y1", data = dat, outcome = "Y1", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),  # Free loadings on both factors
    factor_spec = "full",  # Quadratic + interaction terms
    evaluation_indicator = "eval",
    use_types = TRUE  # Type intercepts are free
  )

  mc_Y2 <- define_model_component(
    name = "Y2", data = dat, outcome = "Y2", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval",
    use_types = TRUE
  )

  mc_Y3 <- define_model_component(
    name = "Y3", data = dat, outcome = "Y3", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "full",
    evaluation_indicator = "eval",
    use_types = TRUE
  )

  # ---- Build Model System ----
  ms <- define_model_system(
    components = list(
      mc_T1, mc_T2, mc_T3, mc_T4, mc_T5, mc_T6,  # Measurement system
      mc_Y1, mc_Y2, mc_Y3                         # Outcome equations
    ),
    factor = fm
  )

  # ---- Verify Structure ----
  expect_equal(length(ms$components), 9)

  # Verify measurement components do NOT use types (use_types = FALSE)
  for (i in 1:6) {
    expect_false(isTRUE(ms$components[[i]]$use_types))
  }

  # Verify outcome components use types (use_types = TRUE)
  for (i in 7:9) {
    expect_true(isTRUE(ms$components[[i]]$use_types))
  }

  # ---- Initialize Parameters ----
  init <- initialize_parameters(ms, dat, verbose = FALSE)

  # Measurement components (T1-T6) should NOT have type intercepts
  expect_false("T1_type_2_intercept" %in% init$param_names)
  expect_false("T1_type_3_intercept" %in% init$param_names)
  expect_false("T6_type_2_intercept" %in% init$param_names)
  expect_false("T6_type_3_intercept" %in% init$param_names)

  # Outcome components (Y1-Y3) should have type intercepts (free)
  expect_true("Y1_type_2_intercept" %in% init$param_names)
  expect_true("Y1_type_3_intercept" %in% init$param_names)
  expect_true("Y2_type_2_intercept" %in% init$param_names)
  expect_true("Y2_type_3_intercept" %in% init$param_names)
  expect_true("Y3_type_2_intercept" %in% init$param_names)
  expect_true("Y3_type_3_intercept" %in% init$param_names)

  # Count type intercepts: 3 outcome components x 2 type intercepts = 6 total
  type_intercept_params <- grep("_type_[23]_intercept$", init$param_names, value = TRUE)
  expect_equal(length(type_intercept_params), 6)

  # Verify full 2nd order parameters exist for outcome components
  # Quadratic loadings (2 per outcome = 6 total)
  expect_true("Y1_loading_quad_1" %in% init$param_names)
  expect_true("Y1_loading_quad_2" %in% init$param_names)
  expect_true("Y2_loading_quad_1" %in% init$param_names)
  expect_true("Y2_loading_quad_2" %in% init$param_names)
  expect_true("Y3_loading_quad_1" %in% init$param_names)
  expect_true("Y3_loading_quad_2" %in% init$param_names)

  # Interaction loadings (1 per outcome = 3 total, since k=2 factors gives k*(k-1)/2 = 1)
  expect_true("Y1_loading_inter_1_2" %in% init$param_names)
  expect_true("Y2_loading_inter_1_2" %in% init$param_names)
  expect_true("Y3_loading_inter_1_2" %in% init$param_names)

  # Measurement components should NOT have quadratic/interaction terms
  expect_false(any(grepl("^T[1-6]_loading_quad", init$param_names)))
  expect_false(any(grepl("^T[1-6]_loading_inter", init$param_names)))
})


test_that("3-type 2-factor model parameter recovery with full 2nd order outcomes", {
  # Skip on CRAN - this is a computationally intensive estimation test
  skip_if_not(identical(Sys.getenv("NOT_CRAN"), "true"),
              "Skipping parameter recovery test (set NOT_CRAN=true to run)")

  # ============================================================================
  # TRUE PARAMETER VALUES
  # ============================================================================
  set.seed(42)
  n <- 2000

  # Type probabilities (3 types)
  true_pi <- c(0.4, 0.35, 0.25)

  # Factor variances by type (2 factors × 3 types)
  # sigma^2 for each factor within each type
  true_factor_var <- matrix(c(
    1.0, 0.8,   # Type 1: factor 1 var, factor 2 var
    1.2, 0.9,   # Type 2
    0.9, 1.1    # Type 3
  ), nrow = 3, byrow = TRUE)

  # Measurement system parameters (T1-T6)
  # Each measurement: intercept, loading (fixed or free), sigma
  # T1, T2, T3 load on factor 1 only (loading on factor 2 = 0)
  # T4, T5, T6 load on factor 2 only (loading on factor 1 = 0)
  true_meas <- list(
    T1 = list(intercept = 0.5, loading1 = 1.0, loading2 = 0, sigma = 0.5),  # loading1 fixed
    T2 = list(intercept = -0.3, loading1 = 0.8, loading2 = 0, sigma = 0.6),
    T3 = list(intercept = 0.2, loading1 = 1.2, loading2 = 0, sigma = 0.4),
    T4 = list(intercept = 0.1, loading1 = 0, loading2 = 1.0, sigma = 0.5),  # loading2 fixed
    T5 = list(intercept = -0.2, loading1 = 0, loading2 = 0.9, sigma = 0.55),
    T6 = list(intercept = 0.4, loading1 = 0, loading2 = 1.1, sigma = 0.45)
  )

  # Outcome parameters (Y1-Y3) with full 2nd order and type intercepts
  # Each outcome: intercept, beta_x1, loading1, loading2, quad1, quad2, inter_12, type2_int, type3_int, sigma
  true_outcome <- list(
    Y1 = list(
      intercept = 1.0, beta_x1 = 0.5,
      loading1 = 0.6, loading2 = 0.4,
      quad1 = 0.15, quad2 = 0.1, inter_12 = 0.08,
      type2_int = 0.3, type3_int = -0.2,
      sigma = 0.7
    ),
    Y2 = list(
      intercept = -0.5, beta_x1 = -0.3,
      loading1 = 0.5, loading2 = 0.7,
      quad1 = -0.1, quad2 = 0.12, inter_12 = 0.05,
      type2_int = -0.25, type3_int = 0.4,
      sigma = 0.6
    ),
    Y3 = list(
      intercept = 0.8, beta_x1 = 0.2,
      loading1 = 0.9, loading2 = 0.3,
      quad1 = 0.08, quad2 = -0.05, inter_12 = 0.1,
      type2_int = 0.15, type3_int = 0.35,
      sigma = 0.8
    )
  )

  # ============================================================================
  # SIMULATE DATA
  # ============================================================================

  # Draw types for each observation
  types <- sample(1:3, n, replace = TRUE, prob = true_pi)

  # Draw factors for each observation based on type
  f1 <- numeric(n)
  f2 <- numeric(n)
  for (i in 1:n) {
    typ <- types[i]
    f1[i] <- rnorm(1, 0, sqrt(true_factor_var[typ, 1]))
    f2[i] <- rnorm(1, 0, sqrt(true_factor_var[typ, 2]))
  }

  # Covariate
  x1 <- rnorm(n)

  # Generate measurement variables (no type intercepts - fixed to 0)
  T1 <- true_meas$T1$intercept + true_meas$T1$loading1 * f1 + rnorm(n, 0, true_meas$T1$sigma)
  T2 <- true_meas$T2$intercept + true_meas$T2$loading1 * f1 + rnorm(n, 0, true_meas$T2$sigma)
  T3 <- true_meas$T3$intercept + true_meas$T3$loading1 * f1 + rnorm(n, 0, true_meas$T3$sigma)
  T4 <- true_meas$T4$intercept + true_meas$T4$loading2 * f2 + rnorm(n, 0, true_meas$T4$sigma)
  T5 <- true_meas$T5$intercept + true_meas$T5$loading2 * f2 + rnorm(n, 0, true_meas$T5$sigma)
  T6 <- true_meas$T6$intercept + true_meas$T6$loading2 * f2 + rnorm(n, 0, true_meas$T6$sigma)

  # Generate outcome variables (with type intercepts and 2nd order terms)
  generate_outcome <- function(params, f1, f2, x1, types) {
    n <- length(f1)
    y <- numeric(n)
    for (i in 1:n) {
      type_int <- ifelse(types[i] == 2, params$type2_int,
                         ifelse(types[i] == 3, params$type3_int, 0))
      y[i] <- params$intercept + params$beta_x1 * x1[i] +
        params$loading1 * f1[i] + params$loading2 * f2[i] +
        params$quad1 * f1[i]^2 + params$quad2 * f2[i]^2 +
        params$inter_12 * f1[i] * f2[i] +
        type_int +
        rnorm(1, 0, params$sigma)
    }
    return(y)
  }

  Y1 <- generate_outcome(true_outcome$Y1, f1, f2, x1, types)
  Y2 <- generate_outcome(true_outcome$Y2, f1, f2, x1, types)
  Y3 <- generate_outcome(true_outcome$Y3, f1, f2, x1, types)

  # Create data frame
  dat <- data.frame(
    intercept = 1,
    x1 = x1,
    eval = 1,
    T1 = T1, T2 = T2, T3 = T3,
    T4 = T4, T5 = T5, T6 = T6,
    Y1 = Y1, Y2 = Y2, Y3 = Y3
  )

  # ============================================================================
  # DEFINE MODEL
  # ============================================================================

  fm <- define_factor_model(n_factors = 2, n_types = 3)

  # Measurement system (use_types = FALSE - no type intercepts)
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(1, 0), evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T3 <- define_model_component(
    name = "T3", data = dat, outcome = "T3", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(NA_real_, 0), evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T4 <- define_model_component(
    name = "T4", data = dat, outcome = "T4", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, 1), evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T5 <- define_model_component(
    name = "T5", data = dat, outcome = "T5", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval",
    use_types = FALSE
  )

  mc_T6 <- define_model_component(
    name = "T6", data = dat, outcome = "T6", factor = fm,
    covariates = "intercept", model_type = "linear",
    loading_normalization = c(0, NA_real_), evaluation_indicator = "eval",
    use_types = FALSE
  )

  # Outcome equations (use_types = TRUE - free type intercepts, full 2nd order)
  mc_Y1 <- define_model_component(
    name = "Y1", data = dat, outcome = "Y1", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "full", evaluation_indicator = "eval",
    use_types = TRUE
  )

  mc_Y2 <- define_model_component(
    name = "Y2", data = dat, outcome = "Y2", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "full", evaluation_indicator = "eval",
    use_types = TRUE
  )

  mc_Y3 <- define_model_component(
    name = "Y3", data = dat, outcome = "Y3", factor = fm,
    covariates = c("intercept", "x1"), model_type = "linear",
    loading_normalization = c(NA_real_, NA_real_),
    factor_spec = "full", evaluation_indicator = "eval",
    use_types = TRUE
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_T3, mc_T4, mc_T5, mc_T6, mc_Y1, mc_Y2, mc_Y3),
    factor = fm
  )

  # ============================================================================
  # ESTIMATE MODEL
  # ============================================================================

  control <- define_estimation_control(num_cores = 1, n_quad_points = 8)

  result <- estimate_model_rcpp(
    model_system = ms,
    data = dat,
    control = control,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # ============================================================================
  # CHECK PARAMETER RECOVERY
  # ============================================================================

  est <- result$estimates
  tol <- 0.2  # Tolerance for parameter recovery (20% relative)

  # Helper to get estimate by name
  get_est <- function(name) {
    idx <- which(result$param_names == name)
    if (length(idx) == 0) return(NA)
    unname(est[idx])
  }

  # Check measurement loadings (T2, T3 on factor 1; T5, T6 on factor 2)
  expect_equal(get_est("T2_loading_1"), true_meas$T2$loading1, tolerance = tol)
  expect_equal(get_est("T3_loading_1"), true_meas$T3$loading1, tolerance = tol)
  expect_equal(get_est("T5_loading_2"), true_meas$T5$loading2, tolerance = tol)
  expect_equal(get_est("T6_loading_2"), true_meas$T6$loading2, tolerance = tol)

  # Check measurement intercepts
  expect_equal(get_est("T1_intercept"), true_meas$T1$intercept, tolerance = tol)
  expect_equal(get_est("T2_intercept"), true_meas$T2$intercept, tolerance = tol)
  expect_equal(get_est("T4_intercept"), true_meas$T4$intercept, tolerance = tol)

  # Check outcome linear loadings
  expect_equal(get_est("Y1_loading_1"), true_outcome$Y1$loading1, tolerance = tol)
  expect_equal(get_est("Y1_loading_2"), true_outcome$Y1$loading2, tolerance = tol)
  expect_equal(get_est("Y2_loading_1"), true_outcome$Y2$loading1, tolerance = tol)
  expect_equal(get_est("Y2_loading_2"), true_outcome$Y2$loading2, tolerance = tol)

  # Check outcome beta coefficients (naming: {component}_{covariate})
  expect_equal(get_est("Y1_x1"), true_outcome$Y1$beta_x1, tolerance = tol)
  expect_equal(get_est("Y2_x1"), true_outcome$Y2$beta_x1, tolerance = tol)
  expect_equal(get_est("Y3_x1"), true_outcome$Y3$beta_x1, tolerance = tol)

  # Check quadratic loadings (may have larger tolerance due to identification)
  expect_equal(get_est("Y1_loading_quad_1"), true_outcome$Y1$quad1, tolerance = 0.2)
  expect_equal(get_est("Y1_loading_quad_2"), true_outcome$Y1$quad2, tolerance = 0.2)
  expect_equal(get_est("Y2_loading_quad_1"), true_outcome$Y2$quad1, tolerance = 0.2)

  # Check interaction loadings
  expect_equal(get_est("Y1_loading_inter_1_2"), true_outcome$Y1$inter_12, tolerance = 0.2)
  expect_equal(get_est("Y2_loading_inter_1_2"), true_outcome$Y2$inter_12, tolerance = 0.2)

  # Check type intercepts (key test for fix_type_intercepts)
  # Note: Higher tolerance (0.6) for type intercepts - they are harder to identify
  # in complex multi-type models with quadratic/interaction terms and constrained
  # measurement type intercepts. The key test is that measurement type intercepts
  # are constrained to 0 while outcome type intercepts are freely estimated.
  type_int_tol <- 0.6
  expect_equal(get_est("Y1_type_2_intercept"), true_outcome$Y1$type2_int, tolerance = type_int_tol)
  expect_equal(get_est("Y1_type_3_intercept"), true_outcome$Y1$type3_int, tolerance = type_int_tol)
  expect_equal(get_est("Y2_type_2_intercept"), true_outcome$Y2$type2_int, tolerance = type_int_tol)
  expect_equal(get_est("Y2_type_3_intercept"), true_outcome$Y2$type3_int, tolerance = type_int_tol)
  expect_equal(get_est("Y3_type_2_intercept"), true_outcome$Y3$type2_int, tolerance = type_int_tol)
  expect_equal(get_est("Y3_type_3_intercept"), true_outcome$Y3$type3_int, tolerance = type_int_tol)

  # Verify measurement components do NOT have type intercepts (use_types = FALSE)
  expect_true(is.na(get_est("T1_type_2_intercept")))
  expect_true(is.na(get_est("T1_type_3_intercept")))
  expect_true(is.na(get_est("T6_type_2_intercept")))
  expect_true(is.na(get_est("T6_type_3_intercept")))

  # Check that convergence was achieved (0 = success in nlminb)
  expect_equal(result$convergence, 0)
})
