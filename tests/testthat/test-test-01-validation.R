#Validates the inputs: Valid outcome, covariate, and eval_indicator type

test_that("define_model_component validates inputs", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)

  # good path : calls functino with valid arguments
  mc <- define_model_component("Y1", dat, "Y", fm,
                               evaluation_indicator = "eval_y1",
                               covariates = c("X1"),
                               model_type = "linear",
                               intercept = FALSE)
  expect_s3_class(mc, "model_component")

  # missing outcome ("NOPE" outcome name, doesn't exist)
  expect_error(
    define_model_component("bad", dat, "NOPE", fm,
                           covariates = "X1", model_type = "linear",
                           intercept = FALSE),
    regexp = "Outcome variable.*not found"
  )

  # missing covariate ("NOPE covariate doesn't exist)
  expect_error(
    define_model_component("bad", dat, "Y", fm,
                           covariates = "NOPE", model_type = "linear"),
    regexp = "Covariates not found"
  )

  # bad eval indicator type (eval_bad is a string instead of bool)
  dat2 <- dat; dat2$eval_bad <- "oops"
  expect_error(
    define_model_component("bad", dat2, "Y", fm,
                           evaluation_indicator = "eval_bad",
                           covariates = "X1", model_type = "linear",
                           intercept = FALSE),
    regexp = "evaluation_indicator"
  )
})

test_that("define_model_component detects multicollinearity", {
  set.seed(123)
  n <- 100

  # Create data with perfectly collinear columns
  dat <- data.frame(
    Y = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n),
    eval = 1
  )
  # X3 is a linear combination of X1 and X2
  dat$X3 <- 2 * dat$X1 + 3 * dat$X2

  fm <- define_factor_model(n_factors = 1, n_types = 1)

  # Should error on perfect collinearity
  expect_error(
    define_model_component("Y", dat, "Y", fm,
                           covariates = c("X1", "X2", "X3"),
                           model_type = "linear",
                           evaluation_indicator = "eval",
                           intercept = FALSE),
    regexp = "rank deficient|multicollinearity"
  )

  # Should work without the collinear column
  mc <- define_model_component("Y", dat, "Y", fm,
                               covariates = c("X1", "X2"),
                               model_type = "linear",
                               evaluation_indicator = "eval",
                               intercept = FALSE)
  expect_s3_class(mc, "model_component")
})
