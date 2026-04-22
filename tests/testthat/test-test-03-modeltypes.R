#make_toycreates a Y that's continuous, not binary. the test should throw an error
#saying that it's supposed to be either 0 or 1. this enforces the model type constraint

test_that("probit requires 0/1 outcome", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)
  expect_error(
    define_model_component("sel", dat, "Y", fm,
                           covariates = "Z1", model_type = "probit",
                           intercept = FALSE),
    regexp = "0/1"
  )
})

test_that("ordered probit returns J-1 ordered cutpoints", {
  skip_if_not_installed("MASS")
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)

  # make an ordinal outcome with 4 categories
  dat$Yord <- cut(dat$Y,
                  breaks = quantile(dat$Y, probs = seq(0, 1, by = 0.25)),
                  include.lowest = TRUE, ordered_result = TRUE)

  mc <- define_model_component("Y_ord", dat, "Yord", fm,
                               evaluation_indicator = "eval_y1",
                               covariates = "X1",
                               model_type = "oprobit",
                               num_choices = 4,
                               intercept = FALSE)

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  expect_true(is.numeric(ini$init_params))
  # Should have covariates (1) + loadings (1) + thresholds (3) = 5 params
  # Plus factor variance (1) = 6 total
  expect_length(ini$init_params, 6)
})

test_that("oprobit works with multi-factor loading normalization", {
  set.seed(42)
  dat <- make_toy()

  # create an ordered outcome with 4 categories
  z <- dat$Y
  dat$Yord <- cut(z,
                  breaks = quantile(z, seq(0, 1, 0.25), na.rm = TRUE),
                  include.lowest = TRUE, ordered_result = TRUE)

  # two factors
  fm <- define_factor_model(2, 1)

  # Fix the 2nd loading to 1, leave the 1st free (NA) at component level
  mc <- define_model_component(
    "Y_ord2", dat, "Yord", fm,
    evaluation_indicator = "eval_y1",
    covariates = "X1",
    model_type = "oprobit",
    num_choices = 4,
    loading_normalization = c(NA, 1),
    intercept = FALSE
  )

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  # Verify model type is oprobit (outcome conversion to ordered happens internally)
  expect_equal(mc$model_type, "oprobit")

  # Check loading normalization in component
  expect_length(mc$loading_normalization, 2)
  expect_true(is.na(mc$loading_normalization[1]))
  expect_equal(mc$loading_normalization[2], 1)
})
