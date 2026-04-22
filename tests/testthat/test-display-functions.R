# Smoke tests for display functions (print, summary, tables)
# These tests verify that display functions run without error

# Helper to create a simple model result for testing display functions
create_test_result <- function() {
  set.seed(123)
  n <- 200

  f <- rnorm(n, 0, 1)
  y1 <- 2.0 + 1.0 * f + rnorm(n, 0, 0.5)
  y2 <- 1.5 + 0.8 * f + rnorm(n, 0, 0.6)

  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)

  ms <- define_model_system(components = list(mc1, mc2), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

  result <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                 optimizer = "nlminb", verbose = FALSE)
  return(result)
}

test_that("print.factorana_result works without error", {
  result <- create_test_result()

  # Capture output to verify it doesn't error
  output <- capture.output(print(result))

  expect_true(length(output) > 0)
  expect_true(any(grepl("Factor Model Estimation Results", output)))
  expect_true(any(grepl("Convergence", output)))
  expect_true(any(grepl("Log-lik", output)))
})

test_that("summary.factorana_result works without error", {
  result <- create_test_result()

  # Get summary
  s <- summary(result)

  expect_s3_class(s, "summary.factorana_result")
  expect_true(!is.null(s$coefficients))
  expect_true(!is.null(s$loglik))
  expect_true(!is.null(s$convergence))

  # Print summary
  output <- capture.output(print(s))
  expect_true(length(output) > 0)
  expect_true(any(grepl("Factor Model Estimation Results", output)))
})

test_that("components_table returns data.frame", {
  result <- create_test_result()

  tbl <- components_table(result)

  expect_s3_class(tbl, "data.frame")
  expect_s3_class(tbl, "components_table")
  expect_true(ncol(tbl) >= 2)  # At least Parameter column + 1 component
  expect_true("Parameter" %in% names(tbl))

  # Check printing works
  output <- capture.output(print(tbl))
  expect_true(length(output) > 0)
})

test_that("components_to_latex produces valid LaTeX", {
  result <- create_test_result()

  latex <- components_to_latex(result)

  expect_type(latex, "character")
  expect_true(grepl("\\\\begin\\{table\\}", latex))
  expect_true(grepl("\\\\end\\{table\\}", latex))
  expect_true(grepl("\\\\begin\\{tabular\\}", latex))
  expect_true(grepl("\\\\end\\{tabular\\}", latex))
})

test_that("results_table works with multiple models", {
  result1 <- create_test_result()

  # Create slightly different model for comparison
  set.seed(456)
  n <- 200
  f <- rnorm(n, 0, 1)
  y1 <- 2.0 + 1.0 * f + rnorm(n, 0, 0.5)
  y2 <- 1.5 + 0.8 * f + rnorm(n, 0, 0.6)
  dat <- data.frame(y1 = y1, y2 = y2, intercept = 1)

  fm <- define_factor_model(n_factors = 1)
  mc1 <- define_model_component("m1", dat, "y1", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = 1)
  mc2 <- define_model_component("m2", dat, "y2", fm,
                                 covariates = "intercept", model_type = "linear",
                                 loading_normalization = NA_real_)
  ms <- define_model_system(components = list(mc1, mc2), factor = fm)
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  result2 <- estimate_model_rcpp(ms, dat, control = ctrl, parallel = FALSE,
                                  optimizer = "nlminb", verbose = FALSE)

  tbl <- results_table(result1, result2, model_names = c("Model 1", "Model 2"))

  expect_s3_class(tbl, "data.frame")
  expect_s3_class(tbl, "factorana_table")
  expect_true("Model 1" %in% names(tbl))
  expect_true("Model 2" %in% names(tbl))

  # Check printing works
  output <- capture.output(print(tbl))
  expect_true(length(output) > 0)
  expect_true(any(grepl("Factor Model Comparison", output)))
})

test_that("results_to_latex produces valid LaTeX", {
  result1 <- create_test_result()

  latex <- results_to_latex(result1, model_names = c("Baseline"))

  expect_type(latex, "character")
  expect_true(grepl("\\\\begin\\{table\\}", latex))
  expect_true(grepl("\\\\end\\{table\\}", latex))
  expect_true(grepl("\\\\begin\\{tabular\\}", latex))
  expect_true(grepl("\\\\end\\{tabular\\}", latex))
  expect_true(grepl("Baseline", latex))
})

test_that("results_to_latex with custom caption and label", {
  result <- create_test_result()

  latex <- results_to_latex(result,
                            caption = "My Table Caption",
                            label = "tab:my_table")

  expect_true(grepl("\\\\caption\\{My Table Caption\\}", latex))
  expect_true(grepl("\\\\label\\{tab:my_table\\}", latex))
})

test_that("components_table handles different digit settings", {
  result <- create_test_result()

  tbl2 <- components_table(result, digits = 2)
  tbl5 <- components_table(result, digits = 5)

  # Both should work
  expect_s3_class(tbl2, "components_table")
  expect_s3_class(tbl5, "components_table")
})

test_that("results_table handles stars = FALSE", {
  result <- create_test_result()

  tbl <- results_table(result, stars = FALSE)

  expect_s3_class(tbl, "factorana_table")
  # Should still work without significance stars
  output <- capture.output(print(tbl))
  expect_true(length(output) > 0)
})

test_that("results_table handles se_format options", {
  result <- create_test_result()

  tbl_parens <- results_table(result, se_format = "parentheses")
  tbl_brackets <- results_table(result, se_format = "brackets")
  tbl_below <- results_table(result, se_format = "below")

  expect_s3_class(tbl_parens, "factorana_table")
  expect_s3_class(tbl_brackets, "factorana_table")
  expect_s3_class(tbl_below, "factorana_table")
})
