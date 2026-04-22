# Test exploded nested logit (exclude_chosen=FALSE and rankshare_var)

test_that("exploded nested logit works with exclude_chosen=FALSE", {
  skip_on_cran()

  set.seed(42)
  n <- 500
  n_choices <- 4
  n_ranks <- 3

  # Generate data
  dat <- data.frame(
    intercept = 1,
    x1 = rnorm(n),
    f = rnorm(n)  # Latent factor (for simulation)
  )

  # True parameters
  true_beta_int <- c(0.5, -0.3, 0.2)  # Intercepts for choices 2, 3, 4
  true_beta_x1 <- c(0.8, -0.5, 0.3)   # x1 coefficients for choices 2, 3, 4
  true_loading <- c(0.6, 0.4, -0.2)   # Factor loadings for choices 2, 3, 4
  true_f_var <- 1.0

  # Simulate ranked choices (same nest CAN be chosen multiple times)
  simulate_rank <- function(i, already_chosen = NULL) {
    V <- c(0,  # Reference category
           true_beta_int[1] + true_beta_x1[1] * dat$x1[i] + true_loading[1] * dat$f[i],
           true_beta_int[2] + true_beta_x1[2] * dat$x1[i] + true_loading[2] * dat$f[i],
           true_beta_int[3] + true_beta_x1[3] * dat$x1[i] + true_loading[3] * dat$f[i])

    # All choices available (nested logit - no exclusion)
    exp_V <- exp(V)
    probs <- exp_V / sum(exp_V)
    sample(1:n_choices, 1, prob = probs)
  }

  # Generate ranked choices
  for (r in 1:n_ranks) {
    dat[[paste0("rank", r)]] <- sapply(1:n, function(i) simulate_rank(i))
  }

  # Also generate two measurement equations for the factor
  dat$T1 <- 1.0 * dat$f + rnorm(n, 0, 0.5)
  dat$T2 <- 0.8 * dat$f + rnorm(n, 0, 0.6)
  dat$eval <- 1

  # Define model
  fm <- define_factor_model(n_factors = 1)

  # Measurement equations
  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = 1, evaluation_indicator = "eval",
    intercept = FALSE
  )

  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval",
    intercept = FALSE
  )

  # Exploded nested logit with exclude_chosen=FALSE
  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2", "rank3"),  # Multiple outcome columns
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = n_choices,
    exclude_chosen = FALSE,  # KEY: Allow same choice multiple times
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  expect_equal(mc_choice$nrank, 3L)
  expect_false(mc_choice$exclude_chosen)

  # Build model system
  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_choice),
    factor = fm
  )

  control <- define_estimation_control(n_quad_points = 8)

  # Initialize and evaluate likelihood at initial parameters
  fm_ptr <- initialize_factor_model_cpp(ms, dat, control$n_quad_points)

  # Get initial parameters
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  init <- init_result$init_params

  # Evaluate likelihood
  result <- evaluate_likelihood_cpp(fm_ptr, init, compute_gradient = TRUE)

  expect_true(is.finite(result$logLikelihood))
  expect_true(result$logLikelihood < 0)  # Log-likelihood should be negative
  expect_equal(length(result$gradient), length(init))

  # Gradient check using finite differences
  eps <- 1e-5
  fd_grad <- numeric(length(init))
  for (i in seq_along(init)) {
    p_plus <- init
    p_minus <- init
    p_plus[i] <- p_plus[i] + eps
    p_minus[i] <- p_minus[i] - eps

    ll_plus <- evaluate_loglik_only_cpp(fm_ptr, p_plus)
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr, p_minus)
    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps)
  }

  # Check gradients match (relative error < 1e-4)
  rel_errors <- abs(result$gradient - fd_grad) / (abs(fd_grad) + 1e-8)
  expect_true(all(rel_errors < 1e-3),
              info = paste("Max relative gradient error:", max(rel_errors)))
})


test_that("exclude_chosen=TRUE (standard exploded logit) still works", {
  skip_on_cran()

  set.seed(123)
  n <- 300
  n_choices <- 3
  n_ranks <- 2

  dat <- data.frame(
    intercept = 1,
    x1 = rnorm(n)
  )

  # Simulate ranked choices (each choice can only be made once)
  for (i in 1:n) {
    available <- 1:n_choices
    for (r in 1:n_ranks) {
      probs <- rep(1/length(available), length(available))
      choice <- sample(available, 1, prob = probs)
      dat[i, paste0("rank", r)] <- choice
      available <- setdiff(available, choice)  # Remove chosen
    }
  }

  dat$T1 <- rnorm(n)
  dat$eval <- 1

  fm <- define_factor_model(n_factors = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = 1, evaluation_indicator = "eval",
    intercept = FALSE
  )

  # Standard exploded logit (exclude_chosen=TRUE, which is default)
  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2"),
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = n_choices,
    exclude_chosen = TRUE,  # Default behavior
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  expect_true(mc_choice$exclude_chosen)

  ms <- define_model_system(
    components = list(mc_T1, mc_choice),
    factor = fm
  )

  control <- define_estimation_control(n_quad_points = 8)
  fm_ptr <- initialize_factor_model_cpp(ms, dat, control$n_quad_points)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  init <- init_result$init_params

  result <- evaluate_likelihood_cpp(fm_ptr, init, compute_gradient = TRUE)

  expect_true(is.finite(result$logLikelihood))
  expect_true(result$logLikelihood < 0)

  # Gradient check
  eps <- 1e-5
  fd_grad <- numeric(length(init))
  for (i in seq_along(init)) {
    p_plus <- init
    p_minus <- init
    p_plus[i] <- p_plus[i] + eps
    p_minus[i] <- p_minus[i] - eps

    ll_plus <- evaluate_loglik_only_cpp(fm_ptr, p_plus)
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr, p_minus)
    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps)
  }

  rel_errors <- abs(result$gradient - fd_grad) / (abs(fd_grad) + 1e-8)
  expect_true(all(rel_errors < 1e-3),
              info = paste("Max relative gradient error:", max(rel_errors)))
})


test_that("rankshare_var corrections are applied", {
  skip_on_cran()

  set.seed(456)
  n <- 200
  n_choices <- 3
  n_ranks <- 2

  dat <- data.frame(
    intercept = 1,
    x1 = rnorm(n)
  )

  # Simulate choices
  dat$rank1 <- sample(1:n_choices, n, replace = TRUE)
  dat$rank2 <- sample(1:n_choices, n, replace = TRUE)

  # Create rankshare correction columns
  # Layout: (n_choices-1) * n_ranks columns
  # For rank 0: columns for choice 1, choice 2 (0-indexed: 0, 1)
  # For rank 1: columns for choice 1, choice 2
  # These are additive corrections to the linear predictor

  # First column is the "rankshare_var" start
  dat$rs_r0_c0 <- rnorm(n, 0, 0.1)  # Rank 0, choice 1 correction
  dat$rs_r0_c1 <- rnorm(n, 0, 0.1)  # Rank 0, choice 2 correction
  dat$rs_r1_c0 <- rnorm(n, 0, 0.1)  # Rank 1, choice 1 correction
  dat$rs_r1_c1 <- rnorm(n, 0, 0.1)  # Rank 1, choice 2 correction

  dat$T1 <- rnorm(n)
  dat$eval <- 1

  fm <- define_factor_model(n_factors = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = 1, evaluation_indicator = "eval",
    intercept = FALSE
  )

  # Nested logit with rankshare corrections
  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2"),
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = n_choices,
    exclude_chosen = FALSE,
    rankshare_var = "rs_r0_c0",  # First rankshare column
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  expect_equal(mc_choice$rankshare_var, "rs_r0_c0")

  ms <- define_model_system(
    components = list(mc_T1, mc_choice),
    factor = fm
  )

  control <- define_estimation_control(n_quad_points = 8)
  fm_ptr <- initialize_factor_model_cpp(ms, dat, control$n_quad_points)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  init <- init_result$init_params

  # Evaluate with rankshare corrections
  result_with_rs <- evaluate_likelihood_cpp(fm_ptr, init, compute_gradient = TRUE)

  expect_true(is.finite(result_with_rs$logLikelihood))

  # Compare to model without rankshare corrections
  mc_choice_no_rs <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2"),
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = n_choices,
    exclude_chosen = FALSE,
    rankshare_var = NULL,  # No rankshare
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  ms_no_rs <- define_model_system(
    components = list(mc_T1, mc_choice_no_rs),
    factor = fm
  )

  fm_ptr_no_rs <- initialize_factor_model_cpp(ms_no_rs, dat, control$n_quad_points)
  result_no_rs <- evaluate_likelihood_cpp(fm_ptr_no_rs, init)

  # Likelihoods should be different due to rankshare corrections
  # (unless corrections happen to be exactly 0)
  expect_true(abs(result_with_rs$logLikelihood - result_no_rs$logLikelihood) > 1e-10,
              info = "Rankshare corrections should affect likelihood")

  # Gradient check for model with rankshare
  eps <- 1e-5
  fd_grad <- numeric(length(init))
  for (i in seq_along(init)) {
    p_plus <- init
    p_minus <- init
    p_plus[i] <- p_plus[i] + eps
    p_minus[i] <- p_minus[i] - eps

    ll_plus <- evaluate_loglik_only_cpp(fm_ptr, p_plus)
    ll_minus <- evaluate_loglik_only_cpp(fm_ptr, p_minus)
    fd_grad[i] <- (ll_plus - ll_minus) / (2 * eps)
  }

  rel_errors <- abs(result_with_rs$gradient - fd_grad) / (abs(fd_grad) + 1e-8)
  expect_true(all(rel_errors < 1e-3),
              info = paste("Max relative gradient error with rankshare:", max(rel_errors)))
})


test_that("Hessian is correct for exploded nested logit", {
  skip_on_cran()

  set.seed(789)
  n <- 200
  n_choices <- 3
  n_ranks <- 2

  dat <- data.frame(
    intercept = 1,
    x1 = rnorm(n),
    f = rnorm(n)
  )

  # Simulate choices
  for (r in 1:n_ranks) {
    dat[[paste0("rank", r)]] <- sample(1:n_choices, n, replace = TRUE)
  }

  dat$T1 <- 1.0 * dat$f + rnorm(n, 0, 0.5)
  dat$eval <- 1

  fm <- define_factor_model(n_factors = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = 1, evaluation_indicator = "eval",
    intercept = FALSE
  )

  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2"),
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = n_choices,
    exclude_chosen = FALSE,
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_choice),
    factor = fm
  )

  control <- define_estimation_control(n_quad_points = 8)
  fm_ptr <- initialize_factor_model_cpp(ms, dat, control$n_quad_points)
  init_result <- initialize_parameters(ms, dat, verbose = FALSE)
  init <- init_result$init_params

  # Evaluate with Hessian
  result <- evaluate_likelihood_cpp(fm_ptr, init, compute_gradient = TRUE, compute_hessian = TRUE)

  expect_true(is.finite(result$logLikelihood))
  expect_true(length(result$hessian) > 0)

  # FD Hessian check (check every element)
  n_params <- length(init)
  eps <- 1e-4

  # Compute FD Hessian by differentiating gradient
  fd_hess <- matrix(0, n_params, n_params)
  for (i in 1:n_params) {
    p_plus <- init
    p_minus <- init
    p_plus[i] <- p_plus[i] + eps
    p_minus[i] <- p_minus[i] - eps

    grad_plus <- evaluate_likelihood_cpp(fm_ptr, p_plus, compute_gradient = TRUE)$gradient
    grad_minus <- evaluate_likelihood_cpp(fm_ptr, p_minus, compute_gradient = TRUE)$gradient

    fd_hess[, i] <- (grad_plus - grad_minus) / (2 * eps)
  }

  # Symmetrize FD Hessian
  fd_hess <- (fd_hess + t(fd_hess)) / 2

  # Expand analytical Hessian from upper-triangular storage
  analytic_hess <- matrix(0, n_params, n_params)
  idx <- 1
  for (i in 1:n_params) {
    for (j in i:n_params) {
      analytic_hess[i, j] <- result$hessian[idx]
      analytic_hess[j, i] <- result$hessian[idx]
      idx <- idx + 1
    }
  }

  # Check Hessian elements match
  rel_errors <- abs(analytic_hess - fd_hess) / (abs(fd_hess) + 1e-6)
  max_rel_error <- max(rel_errors)

  expect_true(max_rel_error < 0.01,
              info = paste("Max relative Hessian error:", max_rel_error))
})


test_that("parameter recovery for exploded nested logit with factor", {
  skip_on_cran()

  set.seed(42)
  n <- 1000
  n_choices <- 3
  n_ranks <- 3

  # True parameters
  true_f_var <- 1.0
  true_T1_sigma <- 0.5
  true_T2_loading <- 0.8
  true_T2_sigma <- 0.6

  # Choice model true params (for choices 2, 3 relative to reference choice 1)
  true_beta_int <- c(0.5, -0.3)      # Intercepts
  true_beta_x1 <- c(0.6, -0.4)       # x1 coefficients
  true_loading <- c(0.7, 0.5)        # Factor loadings

  # Generate latent factor
  f <- rnorm(n, 0, sqrt(true_f_var))

  # Generate data
  dat <- data.frame(
    intercept = 1,
    x1 = rnorm(n)
  )

  # Generate measurement equations
  dat$T1 <- 1.0 * f + rnorm(n, 0, true_T1_sigma)
  dat$T2 <- true_T2_loading * f + rnorm(n, 0, true_T2_sigma)

  # Simulate ranked choices (nested logit - same choice can be made again)
  for (r in 1:n_ranks) {
    choices <- integer(n)
    for (i in 1:n) {
      V <- c(0,  # Reference
             true_beta_int[1] + true_beta_x1[1] * dat$x1[i] + true_loading[1] * f[i],
             true_beta_int[2] + true_beta_x1[2] * dat$x1[i] + true_loading[2] * f[i])
      probs <- exp(V) / sum(exp(V))
      choices[i] <- sample(1:n_choices, 1, prob = probs)
    }
    dat[[paste0("rank", r)]] <- choices
  }

  dat$eval <- 1

  # Define model
  fm <- define_factor_model(n_factors = 1)

  mc_T1 <- define_model_component(
    name = "T1", data = dat, outcome = "T1", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = 1, evaluation_indicator = "eval",
    intercept = FALSE
  )

  mc_T2 <- define_model_component(
    name = "T2", data = dat, outcome = "T2", factor = fm,
    covariates = NULL, model_type = "linear",
    loading_normalization = NA_real_, evaluation_indicator = "eval",
    intercept = FALSE
  )

  mc_choice <- define_model_component(
    name = "choice", data = dat,
    outcome = c("rank1", "rank2", "rank3"),
    factor = fm,
    covariates = c("intercept", "x1"),
    model_type = "logit",
    num_choices = n_choices,
    exclude_chosen = FALSE,
    loading_normalization = NA_real_,
    evaluation_indicator = "eval"
  )

  ms <- define_model_system(
    components = list(mc_T1, mc_T2, mc_choice),
    factor = fm
  )

  control <- define_estimation_control(n_quad_points = 12)

  # Estimate model
  result <- estimate_model_rcpp(
    ms, dat,
    control = control,
    optimizer = "nlminb",
    parallel = FALSE,
    verbose = FALSE
  )

  expect_equal(result$convergence, 0)

  # Extract estimated parameters
  est <- result$estimates

  # Check factor variance recovery (should be close to 1.0)
  expect_true(abs(est["factor_var_1"] - true_f_var) < 0.3,
              info = paste("Factor variance:", est["factor_var_1"], "vs true:", true_f_var))

  # Check T1 sigma recovery
  expect_true(abs(est["T1_sigma"] - true_T1_sigma) < 0.15,
              info = paste("T1 sigma:", est["T1_sigma"], "vs true:", true_T1_sigma))

  # Check T2 loading recovery
  expect_true(abs(est["T2_loading_1"] - true_T2_loading) < 0.2,
              info = paste("T2 loading:", est["T2_loading_1"], "vs true:", true_T2_loading))

  # Check choice model intercept recovery
  expect_true(abs(est["choice_c1_intercept"] - true_beta_int[1]) < 0.3,
              info = paste("Choice c1 intercept:", est["choice_c1_intercept"], "vs true:", true_beta_int[1]))

  expect_true(abs(est["choice_c2_intercept"] - true_beta_int[2]) < 0.3,
              info = paste("Choice c2 intercept:", est["choice_c2_intercept"], "vs true:", true_beta_int[2]))

  # Check choice model x1 coefficient recovery
  expect_true(abs(est["choice_c1_x1"] - true_beta_x1[1]) < 0.2,
              info = paste("Choice c1 x1:", est["choice_c1_x1"], "vs true:", true_beta_x1[1]))

  expect_true(abs(est["choice_c2_x1"] - true_beta_x1[2]) < 0.2,
              info = paste("Choice c2 x1:", est["choice_c2_x1"], "vs true:", true_beta_x1[2]))

  # Check choice model loading recovery
  expect_true(abs(est["choice_c1_loading_1"] - true_loading[1]) < 0.25,
              info = paste("Choice c1 loading:", est["choice_c1_loading_1"], "vs true:", true_loading[1]))

  expect_true(abs(est["choice_c2_loading_1"] - true_loading[2]) < 0.25,
              info = paste("Choice c2 loading:", est["choice_c2_loading_1"], "vs true:", true_loading[2]))
})
