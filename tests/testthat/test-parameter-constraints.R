test_that("Parameter constraints work correctly", {
  skip_on_cran()

  set.seed(456)
  n <- 200

  # Generate simple data
  f <- rnorm(n)
  T1 <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)
  T2 <- 0.5 + 1.2*f + rnorm(n, 0, 0.6)

  dat <- data.frame(intercept=1, T1=T1, T2=T2, eval=1)

  # Define model with factor variance fixed
  fm <- define_factor_model(n_factors=1, n_types=1)
  mc_T1 <- define_model_component(name="T1", data=dat, outcome="T1", factor=fm,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")
  mc_T2 <- define_model_component(name="T2", data=dat, outcome="T2", factor=fm,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")

  ms <- define_model_system(components=list(mc_T1, mc_T2), factor=fm)

  # Estimate with parameter constraints
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  init_params <- c(1.0, 1.0, 1.0, 0.5, 0.5, 1.2, 0.6)
  result <- estimate_model_rcpp(
    ms, dat,
    init_params = init_params,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Test that estimation converged successfully
  expect_equal(result$convergence, 0)
  expect_true(!is.null(result))
  expect_true(!is.null(result$loglik))
  expect_true(!is.null(result$estimates))

  # Test parameter vector length
  expect_equal(length(result$estimates), 7)

  # Test that factor variance stayed fixed at 1.0
  expect_equal(result$estimates[1], 1.0, tolerance = 1e-10)

  # Test that sigma parameters are positive
  expect_true(result$estimates[4] > 0)  # T1 sigma
  expect_true(result$estimates[7] > 0)  # T2 sigma

  # Test that standard errors were computed
  expect_true(!is.null(result$std_errors))
  expect_equal(length(result$std_errors), 7)

  # Standard error for fixed parameter should be 0 or NA
  expect_true(is.na(result$std_errors[1]) || result$std_errors[1] == 0)
})

test_that("Free parameter optimization works", {
  skip_on_cran()

  set.seed(789)
  n <- 100

  # Simple one-factor model
  f <- rnorm(n)
  Y <- 2.0 + 1.5*f + rnorm(n, 0, 0.3)

  dat <- data.frame(intercept=1, Y=Y, eval=1)

  fm <- define_factor_model(n_factors=1, n_types=1)
  mc <- define_model_component(name="Y", data=dat, outcome="Y", factor=fm,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")

  ms <- define_model_system(components=list(mc), factor=fm)

  # Estimate
  ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  init_params <- c(1.0, 2.0, 1.5, 0.5)
  result <- estimate_model_rcpp(
    ms, dat,
    init_params = init_params,
    control = ctrl,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check convergence
  expect_equal(result$convergence, 0)
  expect_true(!is.null(result$loglik))

  # Check parameters are reasonable
  expect_true(abs(result$estimates[2] - 2.0) < 0.5)  # intercept
  expect_true(abs(result$estimates[3] - 1.5) < 0.5)  # loading
  expect_true(result$estimates[4] > 0)  # sigma
})
