#takes a few rows for the Y1 equation and set Y== NA there, expect function to error with message
# "missing values in outcome"
#enforces that each regression must only use rows where it's supposed to
#aka eval_indicator says so, AND with an observed outcome

test_that("only eval==1 and non-missing Y are used", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)

  # inject NA in Y inside eval_y1 subset -> should fail validation
  idx <- which(dat$eval_y1 == 1L)[1:3]
  dat$Y[idx] <- NA_real_

  expect_error(
    define_model_component("Y1", dat, "Y", fm,
                           evaluation_indicator = "eval_y1",
                           covariates = "X1", model_type = "linear",
                           intercept = FALSE),
    regexp = "Missing values in outcome"
  )
})



#creates a scenario where there are no rows for eval_indicator to be 1
#so basically no rows to evaluate.

test_that("zero rows after conditioning throws a clear error", {
  dat <- make_toy()
  fm  <- define_factor_model(1, 1)

  # make eval_y1 all zeros -> no rows
  dat$eval_y1 <- 0L

  expect_error(
    define_model_component(
      name = "Y1",
      data = dat,
      outcome = "Y",
      factor = fm,
      evaluation_indicator = "eval_y1",
      covariates = "X1",
      model_type = "linear",
      intercept = FALSE
    ),
    regexp = "Evaluation subset has zero rows", fixed = TRUE
  )
})
