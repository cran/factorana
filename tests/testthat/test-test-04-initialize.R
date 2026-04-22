#assert: initialize_parameters returns proper init_params vector
#for linear model with 1 covariate and 1 factor

test_that("initialize_parameters returns expected shapes", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)

  mc <- define_model_component("Y1", dat, "Y", fm,
                               evaluation_indicator = "eval_y1",
                               covariates = c("X1"),
                               model_type = "linear",
                               intercept = FALSE)

  # Create model system for initialization
  ms <- define_model_system(components = list(mc), factor = fm)
  ini <- initialize_parameters(ms, dat)

  # Should return init_params vector
  expect_true(is.numeric(ini$init_params))
  # Should have: factor variance (1) + covariate (1) + loading (1) + sigma (1) = 4 params
  expect_length(ini$init_params, 4)
  expect_true(is.logical(ini$factor_variance_fixed))
})
