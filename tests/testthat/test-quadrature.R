test_that("Gauss-Hermite quadrature accuracy", {
  skip_on_cran()

  set.seed(888)
  n <- 500

  # Generate measurement system with three tests (like Roy model)
  f <- rnorm(n, 0, 1)
  T1 <- 2.0 + 1.0*f + rnorm(n, 0, 0.5)  # Loading fixed to 1.0
  T2 <- 1.5 + 1.2*f + rnorm(n, 0, 0.6)  # Loading ~1.2
  T3 <- 1.0 + 0.8*f + rnorm(n, 0, 0.4)  # Loading ~0.8

  dat <- data.frame(
    intercept = 1,
    T1 = T1, T2 = T2, T3 = T3,
    eval = 1
  )

  # Test with different quadrature points
  for (n_quad in c(8, 16)) {
    fm <- define_factor_model(n_factors=1, n_types=1)

    # Define three measurement components
    mc1 <- define_model_component(name="T1", data=dat, outcome="T1", factor=fm,
      covariates="intercept", model_type="linear",
      loading_normalization=1.0, evaluation_indicator="eval")

    mc2 <- define_model_component(name="T2", data=dat, outcome="T2", factor=fm,
      covariates="intercept", model_type="linear",
      loading_normalization=NA_real_, evaluation_indicator="eval")

    mc3 <- define_model_component(name="T3", data=dat, outcome="T3", factor=fm,
      covariates="intercept", model_type="linear",
      loading_normalization=NA_real_, evaluation_indicator="eval")

    ms <- define_model_system(components=list(mc1, mc2, mc3), factor=fm)

    ctrl <- define_estimation_control(n_quad_points = n_quad, num_cores = 1)
    result <- estimate_model_rcpp(
      ms, dat,
      init_params = NULL,
      control = ctrl,
      optimizer = "nlminb",
      parallel = FALSE,
      verbose = FALSE
    )

    # Check that estimation completed
    expect_true(!is.null(result$loglik))
    expect_equal(result$convergence, 0)

    # Check that loadings are reasonably recovered (within 20%)
    # Parameters: factor_var, T1_intercept, T1_sigma,
    #             T2_intercept, T2_loading, T2_sigma,
    #             T3_intercept, T3_loading, T3_sigma
    T2_loading_idx <- 5
    T3_loading_idx <- 8

    expect_true(abs(result$estimates[T2_loading_idx] - 1.2) < 0.24,
                info=sprintf("T2 loading: %.3f vs true 1.2 (n_quad=%d)",
                           result$estimates[T2_loading_idx], n_quad))
    expect_true(abs(result$estimates[T3_loading_idx] - 0.8) < 0.16,
                info=sprintf("T3 loading: %.3f vs true 0.8 (n_quad=%d)",
                           result$estimates[T3_loading_idx], n_quad))
  }
})

test_that("n_quad parameter is properly passed to C++", {
  skip_on_cran()

  set.seed(456)
  n <- 200

  f <- rnorm(n)
  Y <- 1.0 + 1.0*f + rnorm(n, 0, 0.5)

  dat <- data.frame(intercept=1, Y=Y, eval=1)

  # Test with n_quad=8
  fm8 <- define_factor_model(n_factors=1, n_types=1)
  mc8 <- define_model_component(name="Y", data=dat, outcome="Y", factor=fm8,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")
  ms8 <- define_model_system(components=list(mc8), factor=fm8)

  ctrl8 <- define_estimation_control(n_quad_points = 8, num_cores = 1)
  init_params <- c(1.0, 1.0, 1.0, 0.5)
  result8 <- estimate_model_rcpp(ms8, dat, init_params=init_params,
    control=ctrl8, optimizer="nlminb", parallel=FALSE, verbose=FALSE)

  # Test with n_quad=16
  fm16 <- define_factor_model(n_factors=1, n_types=1)
  mc16 <- define_model_component(name="Y", data=dat, outcome="Y", factor=fm16,
    covariates="intercept", model_type="linear",
    loading_normalization=NA_real_, evaluation_indicator="eval")
  ms16 <- define_model_system(components=list(mc16), factor=fm16)

  ctrl16 <- define_estimation_control(n_quad_points = 16, num_cores = 1)
  result16 <- estimate_model_rcpp(ms16, dat, init_params=init_params,
    control=ctrl16, optimizer="nlminb", parallel=FALSE, verbose=FALSE)

  # Results should be similar but may differ slightly
  # The key is that both converge successfully
  expect_true(!is.null(result8$loglik))
  expect_true(!is.null(result16$loglik))

  # More quadrature points should generally give better approximation
  # (though for simple cases the difference is small)
  expect_true(abs(result8$loglik - result16$loglik) < 5)
})
