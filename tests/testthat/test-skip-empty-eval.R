# Tests for warn-and-skip behavior when a component's evaluation_indicator
# selects zero rows. See NEWS 1.3.1.

test_that("define_model_component returns NULL with a warning when eval is all zero", {
  set.seed(1)
  n <- 50
  dat <- data.frame(intercept = 1, y = rnorm(n),
                    eval_off = rep(0L, n), eval_on = rep(1L, n))
  fm <- define_factor_model(n_factors = 1)

  expect_warning(
    mc_null <- define_model_component(
      "m_empty", dat, "y", fm, covariates = "intercept",
      model_type = "linear", loading_normalization = 1,
      evaluation_indicator = "eval_off"
    ),
    "Component 'm_empty' is skipped.*evaluation_indicator 'eval_off'"
  )
  expect_null(mc_null)
})

test_that("define_model_component returns a component when eval has any TRUE row", {
  set.seed(2)
  n <- 50
  ei <- c(rep(0L, n - 1), 1L)  # only one row passes
  dat <- data.frame(intercept = 1, y = rnorm(n), eval = ei)
  fm <- define_factor_model(n_factors = 1)
  expect_silent(
    mc <- define_model_component(
      "m_one", dat, "y", fm, covariates = "intercept",
      model_type = "linear", loading_normalization = 1,
      evaluation_indicator = "eval"
    )
  )
  expect_s3_class(mc, "model_component")
  expect_equal(mc$n_obs, 1L)
})

test_that("define_model_component behavior is unchanged when no evaluation_indicator", {
  set.seed(3)
  n <- 30
  dat <- data.frame(intercept = 1, y = rnorm(n))
  fm <- define_factor_model(n_factors = 1)
  expect_silent(
    mc <- define_model_component(
      "m_no_eval", dat, "y", fm, covariates = "intercept",
      model_type = "linear", loading_normalization = 1
    )
  )
  expect_s3_class(mc, "model_component")
  expect_equal(mc$n_obs, n)
})

test_that("define_model_component still errors on truly empty input data", {
  dat0 <- data.frame(intercept = double(0), y = double(0))
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_component("m", dat0, "y", fm, covariates = "intercept",
                            model_type = "linear", loading_normalization = 1),
    "Evaluation subset has zero rows"
  )
})

test_that("define_model_system silently drops NULL components", {
  set.seed(4)
  n <- 60
  f <- rnorm(n, 0, 1)
  dat <- data.frame(intercept = 1,
                    y1 = f + rnorm(n, 0, 0.5),
                    y2 = 0.9 * f + rnorm(n, 0, 0.5),
                    y3 = 1.1 * f + rnorm(n, 0, 0.5),
                    eval_off = 0L, eval_on = 1L)
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = 1,
                                 evaluation_indicator = "eval_on")
  mc2 <- suppressWarnings(define_model_component(
    "m2_empty", dat, "y2", fm, covariates = "intercept",
    model_type = "linear", loading_normalization = NA_real_,
    evaluation_indicator = "eval_off"))
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = NA_real_,
                                 evaluation_indicator = "eval_on")
  expect_null(mc2)

  ms <- define_model_system(components = list(mc1, mc2, mc3), factor = fm)
  expect_equal(length(ms$components), 2L)
  expect_true(all(vapply(ms$components, inherits, logical(1),
                          "model_component")))
})

test_that("define_model_system errors when all components are NULL and no previous_stage", {
  fm <- define_factor_model(n_factors = 1)
  expect_error(
    define_model_system(components = list(NULL, NULL), factor = fm),
    "All components were skipped"
  )
})

test_that("fix_coefficient passes NULL through silently", {
  expect_silent(out <- fix_coefficient(NULL, "x", 0))
  expect_null(out)
})

test_that("End-to-end: dropped components do not appear in estimates", {
  skip_on_cran()

  set.seed(7)
  n <- 400
  f <- rnorm(n, 0, 1)
  dat <- data.frame(intercept = 1,
                    y1 = f + rnorm(n, 0, 0.5),
                    y2 = 0.9 * f + rnorm(n, 0, 0.5),
                    y3 = 1.1 * f + rnorm(n, 0, 0.5),
                    eval_off = 0L, eval_on = 1L)
  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = 1,
                                 evaluation_indicator = "eval_on")
  mc_dropped <- suppressWarnings(define_model_component(
    "m_dropped", dat, "y2", fm, covariates = "intercept",
    model_type = "linear", loading_normalization = NA_real_,
    evaluation_indicator = "eval_off"))
  mc3 <- define_model_component("m3", dat, "y3", fm, covariates = "intercept",
                                 model_type = "linear", loading_normalization = NA_real_,
                                 evaluation_indicator = "eval_on")
  ms <- define_model_system(components = list(mc1, mc_dropped, mc3),
                             factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)
  r <- estimate_model_rcpp(ms, dat, control = ctrl, optimizer = "nlminb",
                            parallel = FALSE, verbose = FALSE)
  expect_equal(r$convergence, 0)
  expect_false(any(grepl("^m_dropped_", names(r$estimates))),
               info = paste("Estimates names:",
                            paste(names(r$estimates), collapse = ", ")))
})
