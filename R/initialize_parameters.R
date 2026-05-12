#' Initialize parameters for factor model estimation
#'
#' Estimates each model component separately in R (ignoring factors) to obtain
#' good starting values. Also checks factor identification.
#'
#' @param model_system A model_system object from define_model_system()
#' @param data Data frame containing all variables
#' @param factor_scores Optional matrix of factor scores (nobs x n_factors). When provided,
#'   factor loadings are estimated by including factor scores as regressors in the
#'   initialization regressions. This is useful for two-stage estimation where factor
#'   scores from Stage 1 can be used to initialize loadings in Stage 2 outcomes.
#' @param verbose Whether to print progress (default TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item \code{init_params} - Initial parameter values
#'     \item \code{factor_variance_fixed} - Logical vector indicating which factor variances must be fixed
#'   }
#'
#' @details
#' Factor identification: For each factor, if NO component has a non-zero fixed loading,
#' then the factor variance is not identified and must be fixed to 1.0.
#'
#' When \code{factor_scores} is provided, the initialization will estimate loadings by
#' treating factor scores as additional regressors. For example, for a linear model:
#' \code{Y ~ X + factor_scores} and the coefficients on factor_scores become the loading
#' estimates. This typically provides better starting values than the default (0.5).
#'
#' @export
initialize_parameters <- function(model_system, data, factor_scores = NULL, verbose = TRUE) {

  if (verbose) {
    message("Initializing parameters...")
  }

  # Get n_types from factor model (needed in all code paths)
  n_types <- model_system$factor$n_types
  if (is.null(n_types)) n_types <- 1L

  # Validate factor_scores if provided
  n_factors <- model_system$factor$n_factors
  use_factor_scores <- !is.null(factor_scores)
  if (use_factor_scores) {
    if (!is.matrix(factor_scores)) {
      factor_scores <- as.matrix(factor_scores)
    }
    if (nrow(factor_scores) != nrow(data)) {
      stop(sprintf("factor_scores has %d rows but data has %d rows",
                   nrow(factor_scores), nrow(data)))
    }
    if (ncol(factor_scores) != n_factors) {
      stop(sprintf("factor_scores has %d columns but model has %d factors",
                   ncol(factor_scores), n_factors))
    }
    if (verbose) {
      message("  Using factor scores to estimate loadings")
    }
  }

  # ---- 1. Handle previous_stage if present ----
  if (!is.null(model_system$previous_stage_info)) {
    n_fixed_comps <- model_system$previous_stage_info$n_components
    allow_diff_struct <- isTRUE(model_system$previous_stage_info$allow_different_structure)

    if (allow_diff_struct) {
      # For SE_linear/SE_quadratic Stage 2: construct new parameter vector
      # - Use Stage 2's factor structure (factor_var_1, factor_var_2, SE params)
      # - Keep measurement parameters from Stage 1 (loadings, thresholds)

      n_factors <- model_system$factor$n_factors
      factor_structure <- model_system$factor$factor_structure

      # Start with factor structure parameters for Stage 2
      if (factor_structure == "SE_linear") {
        # Factor variances for input factors only (n_factors - 1)
        init_params <- rep(1.0, n_factors - 1)
        param_names <- paste0("factor_var_", seq_len(n_factors - 1))

        # SE parameters
        n_input_factors <- n_factors - 1
        init_params <- c(init_params, 0.0)  # se_intercept
        param_names <- c(param_names, "se_intercept")

        for (j in seq_len(n_input_factors)) {
          init_params <- c(init_params, 0.5)  # se_linear_j
          param_names <- c(param_names, paste0("se_linear_", j))
        }

        # Type-specific SE intercepts (only when n_types > 1)
        # Order must match C++: between linear coefs and se_residual_var
        if (n_types > 1L) {
          for (t in 2:n_types) {
            init_params <- c(init_params, 0.0)
            param_names <- c(param_names, paste0("se_intercept_type_", t))
          }
        }

        init_params <- c(init_params, 1.0)  # se_residual_var
        param_names <- c(param_names, "se_residual_var")

        if (verbose) {
          message(sprintf("SE_linear structure: f_%d = intercept + linear_coefs * f_1..%d + epsilon",
                          n_factors, n_input_factors))
        }
      } else if (factor_structure == "SE_quadratic") {
        # Similar handling for SE_quadratic
        init_params <- rep(1.0, n_factors - 1)
        param_names <- paste0("factor_var_", seq_len(n_factors - 1))

        n_input_factors <- n_factors - 1
        init_params <- c(init_params, 0.0)  # se_intercept
        param_names <- c(param_names, "se_intercept")

        for (j in seq_len(n_input_factors)) {
          init_params <- c(init_params, 0.5)  # se_linear_j
          param_names <- c(param_names, paste0("se_linear_", j))
        }

        for (j in seq_len(n_input_factors)) {
          init_params <- c(init_params, 0.0)  # se_quadratic_j
          param_names <- c(param_names, paste0("se_quadratic_", j))
        }

        # Type-specific SE intercepts (only when n_types > 1)
        # Order must match C++: between quadratic coefs and se_residual_var
        if (n_types > 1L) {
          for (t in 2:n_types) {
            init_params <- c(init_params, 0.0)
            param_names <- c(param_names, paste0("se_intercept_type_", t))
          }
        }

        init_params <- c(init_params, 1.0)  # se_residual_var
        param_names <- c(param_names, "se_residual_var")
      }

      # Parameter-ordering invariant (see build_parameter_metadata() comment
      # in optimize_model.R): type-model params (typeprob/type_loading) come
      # immediately after the SE block, BEFORE factor_mean and SE covariate
      # params. Putting covariate slots before the type slots desyncs every
      # gradient/Hessian element involving covariate or type params (the C++
      # FactorModel constructor places typeprob/type_loading at type_param_start
      # = nparam, before SetSECovariates / SetFactorMeanCovariates extend
      # nparam further).

      # Add type model parameters when n_types > 1.
      if (n_types > 1L) {
        n_factors_se <- model_system$factor$n_factors
        # Type probability intercepts: n_types - 1 of them
        for (t in 2:n_types) {
          init_params <- c(init_params, 0.0)
          param_names <- c(param_names, paste0("typeprob_", t, "_intercept"))
        }
        # Type probability loadings on each factor: (n_types - 1) * n_factors.
        # The loading on the OUTCOME factor (k = n_factors_se) is fixed to 0
        # by setup_parameter_constraints; we still need a slot for it.
        for (t in 2:n_types) {
          for (k in seq_len(n_factors_se)) {
            init_params <- c(init_params, 0.0)
            param_names <- c(param_names, paste0("type_", t, "_loading_", k))
          }
        }
      }

      # Add factor mean covariate parameters if specified (for two-stage)
      factor_covariates <- model_system$factor$factor_covariates
      if (!is.null(factor_covariates) && length(factor_covariates) > 0) {
        n_input_factors <- n_factors - 1
        for (k in seq_len(n_input_factors)) {
          for (cov_name in factor_covariates) {
            init_params <- c(init_params, 0.0)
            param_names <- c(param_names, paste0("factor_mean_", k, "_", cov_name))
          }
        }
        if (verbose) {
          message(sprintf("Factor mean covariates: %d covariates x %d factors = %d parameters",
                          length(factor_covariates), n_input_factors,
                          length(factor_covariates) * n_input_factors))
        }
      }

      # Add SE covariate parameters if specified (for two-stage SE_linear/SE_quadratic)
      se_covariates <- model_system$factor$se_covariates
      if (!is.null(se_covariates) && length(se_covariates) > 0) {
        for (cov_name in se_covariates) {
          init_params <- c(init_params, 0.0)
          param_names <- c(param_names, paste0("se_cov_", cov_name))
        }
        if (verbose) {
          message(sprintf("SE covariates (two-stage): %d parameters", length(se_covariates)))
          message(sprintf("  Covariates: %s", paste(se_covariates, collapse = ", ")))
        }
      }

      # Add measurement parameters from Stage 1 (only fixed ones)
      # These are loadings, thresholds, intercepts, sigmas, etc. — not
      # factor_var, se_, factor_mean_, chol_, typeprob_, or type_loading_
      # (all factor-level params, already laid out above in the canonical
      # Stage 2 position — keeping them here too would duplicate them).
      prev_params <- model_system$previous_stage_info$all_param_values
      prev_names <- names(prev_params)
      meas_idx <- !grepl("^(factor_var|se_|chol_|factor_mean_|typeprob_|type_\\d+_loading_)",
                         prev_names)
      meas_params <- prev_params[meas_idx]
      meas_names <- prev_names[meas_idx]

      init_params <- c(init_params, meas_params)
      param_names <- c(param_names, meas_names)

      if (verbose) {
        message(sprintf("Using %d measurement parameters from Stage 1", length(meas_params)))
        message(sprintf("Stage 2 has %d factor distribution parameters", sum(!meas_idx) + length(grep("^(factor_var|se_)", param_names[1:10]))))
      }

    } else {
      # Standard previous_stage behavior: use ALL previous-stage parameter values
      init_params <- model_system$previous_stage_info$all_param_values
      param_names <- if (!is.null(model_system$previous_stage_info$param_names)) {
        model_system$previous_stage_info$param_names
      } else {
        names(init_params)
      }
      if (is.null(param_names)) param_names <- character(0)
    }

    if (verbose) {
      n_fixed <- length(model_system$previous_stage_info$fixed_param_names)
      n_free <- length(model_system$previous_stage_info$free_param_names)
      message(sprintf("Using %d parameters from previous stage (%d fixed, %d free)",
                      length(model_system$previous_stage_info$all_param_values), n_fixed, n_free))
      message(sprintf("  Previous stage had %d components", n_fixed_comps))
      if (n_free > 0) {
        message(sprintf("  Free parameters: %s",
                        paste(model_system$previous_stage_info$free_param_names, collapse = ", ")))
      }
    }

    # Initialize only the new (second-stage) components
    start_comp_idx <- n_fixed_comps + 1
  } else {
    # Standard initialization: start from scratch
    n_factors <- model_system$factor$n_factors

    # Check factor identification
    factor_variance_fixed <- rep(FALSE, n_factors)

    for (comp in model_system$components) {
      if (!is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          # Check if this loading is fixed to 1.0 (identification via unit loading)
          if (!is.na(comp$loading_normalization[k]) &&
              abs(comp$loading_normalization[k] - 1.0) < 1e-6) {
            factor_variance_fixed[k] <- TRUE
          }
        }
      }
    }

    if (verbose) {
      message("Factor identification:")
      for (k in seq_len(n_factors)) {
        status <- if (factor_variance_fixed[k]) "identified by fixed loading (variance will be estimated)" else "NOT identified (variance fixed to 1.0)"
        message(sprintf("  Factor %d: %s", k, status))
      }
    }

    # Get mixture count
    n_mixtures <- model_system$factor$n_mixtures
    if (is.null(n_mixtures)) n_mixtures <- 1L

    # Initialize factor variances to 1.0 (one per factor per mixture)
    # For nmix > 1, we have variances for each mixture component
    init_params <- rep(1.0, n_factors * n_mixtures)
    param_names <- character(0)
    for (imix in seq_len(n_mixtures)) {
      for (ifac in seq_len(n_factors)) {
        param_names <- c(param_names,
                         if (n_mixtures == 1) paste0("factor_var_", ifac)
                         else paste0("mix", imix, "_factor_var_", ifac))
      }
    }

    # Add mixture means for nmix > 1 (for non-reference mixtures only)
    # Mean constraint: Σ w_m * μ_m = 0, so last mixture mean is derived
    if (n_mixtures > 1) {
      for (imix in seq_len(n_mixtures - 1)) {
        for (ifac in seq_len(n_factors)) {
          # Initialize means to spread around 0: mix1 positive, mix2 negative, etc.
          # This helps separation of mixture components
          init_mean <- if (imix == 1) 0.5 else if (imix == 2) -0.5 else 0.0
          init_params <- c(init_params, init_mean)
          param_names <- c(param_names, paste0("mix", imix, "_factor_mean_", ifac))
        }
      }

      # Add mixture log-weights for non-reference mixtures
      # Using softmax: w_m = exp(fw_m) / (1 + Σ exp(fw_j))
      # Initialize to 0 for equal weights (all w_m = 1/nmix)
      for (imix in seq_len(n_mixtures - 1)) {
        init_params <- c(init_params, 0.0)  # log-weight = 0 -> equal weights
        param_names <- c(param_names, paste0("mix", imix, "_logweight"))
      }

      if (verbose) {
        message(sprintf("Mixture of normals: %d components", n_mixtures))
        message(sprintf("  %d variance params + %d mean params + %d log-weight params",
                        n_factors * n_mixtures,
                        n_factors * (n_mixtures - 1),
                        n_mixtures - 1))
      }
    }

    # Add structure-specific parameters based on factor_structure
    factor_structure <- model_system$factor$factor_structure
    if (is.null(factor_structure)) factor_structure <- "independent"

    if (factor_structure == "SE_linear") {
      # SE_linear: f_k = alpha + alpha_1*f_1 + ... + epsilon
      # For SE_linear, only input factors (n_factors - 1) have variance parameters
      # Mixtures apply to input factors; epsilon is always single normal
      n_input_factors <- n_factors - 1

      # For SE models, we need to rebuild init_params/param_names to only include input factors
      # Current state: nmix * n_factors variances + (nmix-1)*n_factors means + (nmix-1) log-weights
      # Target state: nmix * n_input_factors variances + (nmix-1)*n_input_factors means + (nmix-1) log-weights + SE params
      init_params_new <- numeric(0)
      param_names_new <- character(0)

      # Input factor variances (one per input factor per mixture)
      for (imix in seq_len(n_mixtures)) {
        for (ifac in seq_len(n_input_factors)) {
          init_params_new <- c(init_params_new, 1.0)
          param_names_new <- c(param_names_new,
                               if (n_mixtures == 1) paste0("factor_var_", ifac)
                               else paste0("mix", imix, "_factor_var_", ifac))
        }
      }

      # Input factor means (for non-reference mixtures)
      if (n_mixtures > 1) {
        for (imix in seq_len(n_mixtures - 1)) {
          for (ifac in seq_len(n_input_factors)) {
            init_mean <- if (imix == 1) 0.5 else if (imix == 2) -0.5 else 0.0
            init_params_new <- c(init_params_new, init_mean)
            param_names_new <- c(param_names_new, paste0("mix", imix, "_factor_mean_", ifac))
          }
        }
        # Mixture log-weights
        for (imix in seq_len(n_mixtures - 1)) {
          init_params_new <- c(init_params_new, 0.0)
          param_names_new <- c(param_names_new, paste0("mix", imix, "_logweight"))
        }
      }

      init_params <- init_params_new
      param_names <- param_names_new

      # Update factor_variance_fixed (only applies to input factors)
      factor_variance_fixed <- factor_variance_fixed[1:n_input_factors]

      # SE parameters: intercept, linear coefficients, residual variance
      # SE intercept (initialize to 0)
      init_params <- c(init_params, 0.0)
      param_names <- c(param_names, "se_intercept")

      # SE linear coefficients (initialize to reasonable values)
      for (j in seq_len(n_input_factors)) {
        init_params <- c(init_params, 0.5)
        param_names <- c(param_names, paste0("se_linear_", j))
      }

      # Type-specific SE intercepts (only when n_types > 1)
      # Order must match C++: between linear coefs and se_residual_var
      if (n_types > 1L) {
        for (t in 2:n_types) {
          init_params <- c(init_params, 0.0)
          param_names <- c(param_names, paste0("se_intercept_type_", t))
        }
      }

      # SE residual variance (initialize to 1.0)
      init_params <- c(init_params, 1.0)
      param_names <- c(param_names, "se_residual_var")

      if (verbose) {
        message(sprintf("SE_linear structure: f_%d = intercept + linear_coefs * f_1..%d + epsilon",
                        n_factors, n_input_factors))
      }

    } else if (factor_structure == "SE_quadratic") {
      # SE_quadratic: f_k = alpha + alpha_1*f_1 + alpha_q1*f_1^2 + ... + epsilon
      # Similar to SE_linear but with quadratic terms added
      n_input_factors <- n_factors - 1

      # Rebuild init_params/param_names for SE_quadratic (same as SE_linear)
      init_params_new <- numeric(0)
      param_names_new <- character(0)

      # Input factor variances (one per input factor per mixture)
      for (imix in seq_len(n_mixtures)) {
        for (ifac in seq_len(n_input_factors)) {
          init_params_new <- c(init_params_new, 1.0)
          param_names_new <- c(param_names_new,
                               if (n_mixtures == 1) paste0("factor_var_", ifac)
                               else paste0("mix", imix, "_factor_var_", ifac))
        }
      }

      # Input factor means (for non-reference mixtures)
      if (n_mixtures > 1) {
        for (imix in seq_len(n_mixtures - 1)) {
          for (ifac in seq_len(n_input_factors)) {
            init_mean <- if (imix == 1) 0.5 else if (imix == 2) -0.5 else 0.0
            init_params_new <- c(init_params_new, init_mean)
            param_names_new <- c(param_names_new, paste0("mix", imix, "_factor_mean_", ifac))
          }
        }
        # Mixture log-weights
        for (imix in seq_len(n_mixtures - 1)) {
          init_params_new <- c(init_params_new, 0.0)
          param_names_new <- c(param_names_new, paste0("mix", imix, "_logweight"))
        }
      }

      init_params <- init_params_new
      param_names <- param_names_new

      # Update factor_variance_fixed (only applies to input factors)
      factor_variance_fixed <- factor_variance_fixed[1:n_input_factors]

      # SE parameters: intercept, linear coefficients, quadratic coefficients, residual variance
      # SE intercept (initialize to 0)
      init_params <- c(init_params, 0.0)
      param_names <- c(param_names, "se_intercept")

      # SE linear coefficients (initialize to reasonable values)
      for (j in seq_len(n_input_factors)) {
        init_params <- c(init_params, 0.5)
        param_names <- c(param_names, paste0("se_linear_", j))
      }

      # SE quadratic coefficients (initialize to 0 - neutral starting point)
      for (j in seq_len(n_input_factors)) {
        init_params <- c(init_params, 0.0)
        param_names <- c(param_names, paste0("se_quadratic_", j))
      }

      # Type-specific SE intercepts (only when n_types > 1)
      # Order must match C++: between quadratic coefs and se_residual_var
      if (n_types > 1L) {
        for (t in 2:n_types) {
          init_params <- c(init_params, 0.0)
          param_names <- c(param_names, paste0("se_intercept_type_", t))
        }
      }

      # SE residual variance (initialize to 1.0)
      init_params <- c(init_params, 1.0)
      param_names <- c(param_names, "se_residual_var")

      if (verbose) {
        message(sprintf("SE_quadratic structure: f_%d = intercept + linear_coefs * f_1..%d + quadratic_coefs * f_1^2..%d + epsilon",
                        n_factors, n_input_factors, n_input_factors))
      }

    } else if (factor_structure == "correlation" && n_factors == 2) {
      # For 2-factor correlated model, add one correlation parameter
      # Initialize to 0 (uncorrelated) as a neutral starting point
      init_params <- c(init_params, 0.0)
      param_names <- c(param_names, "factor_corr_1_2")

    } else if (isTRUE(model_system$factor$correlation) && n_factors == 2) {
      # Backward compatibility: correlation = TRUE
      init_params <- c(init_params, 0.0)
      param_names <- c(param_names, "factor_corr_1_2")
    }

    # Parameter-ordering invariant — see build_parameter_metadata() for the
    # full rationale. Type-model params (typeprob/type_loading) MUST appear
    # immediately after the SE block, BEFORE factor_mean and SE covariate
    # params, because the C++ FactorModel constructor places them there.

    # Add type model parameters if n_types > 1 AND at least one component uses types
    # Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    # (n_types - 1) intercepts + (n_types - 1) * n_factors loadings (type 1 is reference)
    any_uses_types <- any(sapply(model_system$components, function(c) isTRUE(c$use_types)))
    # SE_linear / SE_quadratic with ntyp > 1 implies types at the structural level
    # (via se_intercept_type_{t}), so the type probability model is needed even when
    # no measurement component sets use_types = TRUE.
    .fs_for_types <- model_system$factor$factor_structure
    if (!is.null(.fs_for_types) && .fs_for_types %in% c("SE_linear", "SE_quadratic") &&
        n_types > 1L) {
      any_uses_types <- TRUE
    }
    if (n_types > 1L && any_uses_types) {
      # Type probability intercepts (n_types - 1)
      typeprob_intercepts <- rep(0.0, n_types - 1L)
      typeprob_intercept_names <- paste0("typeprob_", 2:n_types, "_intercept")
      init_params <- c(init_params, typeprob_intercepts)
      param_names <- c(param_names, typeprob_intercept_names)

      # Type probability loadings ((n_types - 1) * n_factors)
      type_loadings <- rep(0.0, (n_types - 1L) * n_factors)
      type_loading_names <- character(0)
      for (t in 2:n_types) {
        for (k in seq_len(n_factors)) {
          type_loading_names <- c(type_loading_names, paste0("type_", t, "_loading_", k))
        }
      }
      init_params <- c(init_params, type_loadings)
      param_names <- c(param_names, type_loading_names)
    }

    # Add factor mean covariate parameters if specified
    # These coefficients shift the factor mean: E[f_k | X] = X * gamma_k
    factor_covariates <- model_system$factor$factor_covariates
    if (!is.null(factor_covariates) && length(factor_covariates) > 0) {
      # Determine which factors get mean covariates
      # For SE models, only input factors get covariates (outcome factor mean is from SE)
      if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
        n_factors_with_mean <- n_factors - 1
      } else {
        n_factors_with_mean <- n_factors
      }

      # Add parameters: factor_mean_<factor>_<covariate>
      # Initialize to 0 (no effect)
      for (k in seq_len(n_factors_with_mean)) {
        for (cov_name in factor_covariates) {
          init_params <- c(init_params, 0.0)
          param_names <- c(param_names, paste0("factor_mean_", k, "_", cov_name))
        }
      }

      if (verbose) {
        message(sprintf("Factor mean covariates: %d covariates x %d factors = %d parameters",
                        length(factor_covariates), n_factors_with_mean,
                        length(factor_covariates) * n_factors_with_mean))
        message(sprintf("  Covariates: %s", paste(factor_covariates, collapse = ", ")))
      }
    }

    # Add SE covariate parameters if specified (for SE_linear/SE_quadratic)
    # These coefficients directly affect the outcome factor: f_k = ... + beta * X + epsilon
    se_covariates <- model_system$factor$se_covariates
    if (!is.null(se_covariates) && length(se_covariates) > 0) {
      # Add parameters: se_cov_<covariate>
      # Initialize to 0 (no effect)
      for (cov_name in se_covariates) {
        init_params <- c(init_params, 0.0)
        param_names <- c(param_names, paste0("se_cov_", cov_name))
      }

      if (verbose) {
        message(sprintf("SE covariates: %d parameters", length(se_covariates)))
        message(sprintf("  Covariates: %s", paste(se_covariates, collapse = ", ")))
      }
    }

    start_comp_idx <- 1
  }

  # ---- 2. Estimate each component separately ----
  # If previous_stage, only estimate new components
  # Skip if no new components to estimate
  n_total_comps <- length(model_system$components)
  if (start_comp_idx <= n_total_comps) {
  for (i_comp in start_comp_idx:n_total_comps) {
    comp <- model_system$components[[i_comp]]

    if (verbose) {
      message(sprintf("\nEstimating component %d (%s)...", i_comp, comp$name))
    }

    # Apply evaluation indicator if present
    if (!is.null(comp$evaluation_indicator)) {
      idx <- data[[comp$evaluation_indicator]] == 1 & !is.na(data[[comp$evaluation_indicator]])
      comp_data <- data[idx, , drop = FALSE]
    } else {
      comp_data <- data
    }

    # Get outcome - for dynamic models, use zeros; otherwise from data
    # For exploded logit (multiple outcomes), use first outcome column for initialization
    outcome_col <- if (length(comp$outcome) > 1) comp$outcome[1] else comp$outcome
    if (outcome_col %in% names(comp_data)) {
      outcome <- comp_data[[outcome_col]]
    } else if (isTRUE(comp$is_dynamic)) {
      # For dynamic models, create zero outcome (the dummy outcome column)
      outcome <- rep(0, nrow(comp_data))
    } else {
      stop(sprintf("Outcome '%s' not found in data for component '%s'",
                   outcome_col, comp$name))
    }

    # Get covariates - create intercept column if needed
    if (length(comp$covariates) > 0) {
      X <- matrix(NA_real_, nrow = nrow(comp_data), ncol = length(comp$covariates))
      colnames(X) <- comp$covariates
      for (cov in comp$covariates) {
        if (cov %in% names(comp_data)) {
          X[, cov] <- comp_data[[cov]]
        } else if (cov == "intercept" || cov == "constant") {
          X[, cov] <- 1
        } else {
          stop(sprintf("Covariate '%s' not found in data for component '%s'",
                       cov, comp$name))
        }
      }
    } else {
      X <- matrix(nrow = nrow(comp_data), ncol = 0)
    }

    # Count free factor loadings for this component
    n_free_loadings <- sum(is.na(comp$loading_normalization))

    # Estimate model depending on type
    comp_params <- NULL
    comp_param_names <- NULL

    # Get fixed coefficients list (may be empty)
    fixed_coefs <- if (!is.null(comp$fixed_coefficients)) comp$fixed_coefficients else list()

    if (comp$model_type == "linear") {
      # Linear regression - optionally include factor scores to estimate loadings
      loading_init <- rep(0.5, n_free_loadings)  # Default loading initialization

      if (use_factor_scores && n_free_loadings > 0 && !isTRUE(comp$is_dynamic)) {
        # Get factor scores for FREE loadings only (subset to relevant factors)
        free_factor_idx <- which(is.na(comp$loading_normalization))
        # Use comp_data rows (respects evaluation_indicator subsetting)
        fs_subset <- factor_scores[, free_factor_idx, drop = FALSE]
        if (!is.null(comp$evaluation_indicator)) {
          eval_mask <- data[[comp$evaluation_indicator]] == 1
          fs_subset <- fs_subset[eval_mask, , drop = FALSE]
        }

        # Combine X and factor scores for regression
        if (ncol(X) > 0) {
          X_full <- cbind(X, fs_subset)
        } else {
          X_full <- fs_subset
        }

        fit <- lm(outcome ~ X_full - 1)
        all_coefs <- coef(fit)
        sigma <- summary(fit)$sigma

        # Split coefficients: first ncol(X) are betas, rest are loadings
        n_x <- ncol(X)
        coefs <- if (n_x > 0) all_coefs[1:n_x] else numeric(0)
        loading_init <- all_coefs[(n_x + 1):length(all_coefs)]
        names(loading_init) <- NULL  # Remove names for clean output

      } else if (ncol(X) > 0) {
        fit <- lm(outcome ~ X - 1)  # -1 because intercept is already in covariates
        coefs <- coef(fit)
        sigma <- summary(fit)$sigma
      } else {
        # No covariates - just compute residual variance
        coefs <- numeric(0)
        sigma <- sd(outcome)
      }

      # For dynamic models (outcome = 0), sigma would be 0 or very small
      # Initialize to 1.0 for a reasonable starting point
      if (isTRUE(comp$is_dynamic)) {
        sigma <- 1.0
        # Also reset coefficients to small values for dynamic models
        if (length(coefs) > 0) {
          coefs <- rep(0.1, length(coefs))
        }
      }

      # Apply any fixed coefficient values
      coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

      if (verbose) {
        n_fixed <- length(fixed_coefs)
        msg <- sprintf("  Linear model: %d covariates, sigma = %.4f", length(coefs), sigma)
        if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
        if (use_factor_scores && n_free_loadings > 0) {
          msg <- paste0(msg, sprintf(" (loadings from factor scores: %s)",
                                     paste(sprintf("%.3f", loading_init), collapse = ", ")))
        }
        message(msg)
      }

      # Build parameter vector: coefs, linear loadings, [second-order loadings], sigma
      # Note: For linear models with intercept=TRUE, "intercept" is automatically
      # added to covariates in define_model_component(), so it's handled as a covariate.
      # Get second-order loading initializations
      second_order <- get_second_order_loading_init(comp)

      comp_params <- c(coefs, loading_init, second_order$values, sigma)

      # Build parameter names
      coef_names <- if (length(comp$covariates) > 0) {
        paste0(comp$name, "_", comp$covariates)
      } else {
        character(0)
      }
      loading_names <- character(0)
      if (n_free_loadings > 0) {
        free_factor_idx <- which(is.na(comp$loading_normalization))
        loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
      }
      comp_param_names <- c(coef_names, loading_names, second_order$names, paste0(comp$name, "_sigma"))

      # Add type-specific intercepts if component uses types and n_types > 1
      # Only components with use_types = TRUE get type intercepts
      if (isTRUE(comp$use_types) && n_types > 1L) {
        type_intercepts <- numeric(n_types - 1L)
        type_intercept_names <- character(n_types - 1L)
        for (t in 2:n_types) {
          idx <- t - 1L
          type_intercept_names[idx] <- paste0(comp$name, "_type_", t, "_intercept")
          if (is_type_intercept_fixed(comp, t, choice = NULL)) {
            type_intercepts[idx] <- 0.0  # Fixed value
          } else {
            type_intercepts[idx] <- 0.1 * (t - 1L)  # 0.1, 0.2, 0.3, ...
          }
        }
        comp_params <- c(comp_params, type_intercepts)
        comp_param_names <- c(comp_param_names, type_intercept_names)
      }

    } else if (comp$model_type == "probit") {
      # Binary probit - optionally include factor scores to estimate loadings
      loading_init <- rep(0.5, n_free_loadings)  # Default loading initialization

      if (use_factor_scores && n_free_loadings > 0) {
        # Get factor scores for FREE loadings only
        free_factor_idx <- which(is.na(comp$loading_normalization))
        fs_subset <- factor_scores[, free_factor_idx, drop = FALSE]
        if (!is.null(comp$evaluation_indicator)) {
          eval_mask <- data[[comp$evaluation_indicator]] == 1
          fs_subset <- fs_subset[eval_mask, , drop = FALSE]
        }

        # Combine X and factor scores for regression
        if (ncol(X) > 0) {
          X_full <- cbind(X, fs_subset)
        } else {
          X_full <- fs_subset
        }

        fit <- glm(outcome ~ X_full - 1, family = binomial(link = "probit"))
        all_coefs <- coef(fit)

        # Split coefficients
        n_x <- ncol(X)
        coefs <- if (n_x > 0) all_coefs[1:n_x] else numeric(0)
        loading_init <- all_coefs[(n_x + 1):length(all_coefs)]
        names(loading_init) <- NULL
      } else {
        fit <- glm(outcome ~ X - 1, family = binomial(link = "probit"))
        coefs <- coef(fit)
      }

      # Apply any fixed coefficient values
      coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

      if (verbose) {
        n_fixed <- length(fixed_coefs)
        msg <- sprintf("  Probit model: %d covariates", length(coefs))
        if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
        if (use_factor_scores && n_free_loadings > 0) {
          msg <- paste0(msg, sprintf(" (loadings from factor scores: %s)",
                                     paste(sprintf("%.3f", loading_init), collapse = ", ")))
        }
        message(msg)
      }

      # Get second-order loading initializations
      second_order <- get_second_order_loading_init(comp)

      comp_params <- c(coefs, loading_init, second_order$values)

      # Build parameter names
      coef_names <- if (length(comp$covariates) > 0) {
        paste0(comp$name, "_", comp$covariates)
      } else {
        character(0)
      }
      loading_names <- character(0)
      if (n_free_loadings > 0) {
        free_factor_idx <- which(is.na(comp$loading_normalization))
        loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
      }
      comp_param_names <- c(coef_names, loading_names, second_order$names)

      # Add type-specific intercepts if component uses types and n_types > 1
      if (isTRUE(comp$use_types) && n_types > 1L) {
        type_intercepts <- numeric(n_types - 1L)
        type_intercept_names <- character(n_types - 1L)
        for (t in 2:n_types) {
          idx <- t - 1L
          type_intercept_names[idx] <- paste0(comp$name, "_type_", t, "_intercept")
          if (is_type_intercept_fixed(comp, t, choice = NULL)) {
            type_intercepts[idx] <- 0.0  # Fixed value
          } else {
            type_intercepts[idx] <- 0.1 * (t - 1L)  # 0.1, 0.2, 0.3, ...
          }
        }
        comp_params <- c(comp_params, type_intercepts)
        comp_param_names <- c(comp_param_names, type_intercept_names)
      }

    } else if (comp$model_type == "logit") {
      if (comp$num_choices == 2) {
        # Binary logit - optionally include factor scores to estimate loadings
        loading_init <- rep(0.5, n_free_loadings)  # Default loading initialization

        # Convert 1/2 coded outcomes to 0/1 for glm (which expects 0/1)
        # The model uses 1-indexed outcomes (1, 2, ..., K) but glm needs 0/1
        outcome_01 <- outcome - 1  # Convert 1/2 to 0/1

        if (use_factor_scores && n_free_loadings > 0) {
          # Get factor scores for FREE loadings only
          free_factor_idx <- which(is.na(comp$loading_normalization))
          fs_subset <- factor_scores[, free_factor_idx, drop = FALSE]
          if (!is.null(comp$evaluation_indicator)) {
            eval_mask <- data[[comp$evaluation_indicator]] == 1
            fs_subset <- fs_subset[eval_mask, , drop = FALSE]
          }

          # Combine X and factor scores for regression
          if (ncol(X) > 0) {
            X_full <- cbind(X, fs_subset)
          } else {
            X_full <- fs_subset
          }

          fit <- glm(outcome_01 ~ X_full - 1, family = binomial(link = "logit"))
          all_coefs <- coef(fit)

          # Split coefficients
          n_x <- ncol(X)
          coefs <- if (n_x > 0) all_coefs[1:n_x] else numeric(0)
          loading_init <- all_coefs[(n_x + 1):length(all_coefs)]
          names(loading_init) <- NULL
        } else {
          fit <- glm(outcome_01 ~ X - 1, family = binomial(link = "logit"))
          coefs <- coef(fit)
        }

        # Apply any fixed coefficient values
        coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

        if (verbose) {
          n_fixed <- length(fixed_coefs)
          msg <- sprintf("  Binary logit: %d covariates", length(coefs))
          if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
          if (use_factor_scores && n_free_loadings > 0) {
            msg <- paste0(msg, sprintf(" (loadings from factor scores: %s)",
                                       paste(sprintf("%.3f", loading_init), collapse = ", ")))
          }
          message(msg)
        }

        # Get second-order loading initializations
        second_order <- get_second_order_loading_init(comp)

        comp_params <- c(coefs, loading_init, second_order$values)

        # Build parameter names
        coef_names <- if (length(comp$covariates) > 0) {
          paste0(comp$name, "_", comp$covariates)
        } else {
          character(0)
        }
        loading_names <- character(0)
        if (n_free_loadings > 0) {
          free_factor_idx <- which(is.na(comp$loading_normalization))
          loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
        }
        comp_param_names <- c(coef_names, loading_names, second_order$names)

        # Add type-specific intercepts if component uses types and n_types > 1
        if (isTRUE(comp$use_types) && n_types > 1L) {
          type_intercepts <- numeric(n_types - 1L)
          type_intercept_names <- character(n_types - 1L)
          for (t in 2:n_types) {
            idx <- t - 1L
            type_intercept_names[idx] <- paste0(comp$name, "_type_", t, "_intercept")
            if (is_type_intercept_fixed(comp, t, choice = NULL)) {
              type_intercepts[idx] <- 0.0  # Fixed value
            } else {
              type_intercepts[idx] <- 0.1 * (t - 1L)  # 0.1, 0.2, 0.3, ...
            }
          }
          comp_params <- c(comp_params, type_intercepts)
          comp_param_names <- c(comp_param_names, type_intercept_names)
        }

      } else {
        # Multinomial logit
        # Prepare factor scores if using them for loading estimation
        loading_init_by_choice <- NULL
        if (use_factor_scores && n_free_loadings > 0) {
          # Get factor scores for FREE loadings only
          free_factor_idx <- which(is.na(comp$loading_normalization))
          fs_subset <- factor_scores[, free_factor_idx, drop = FALSE]
          if (!is.null(comp$evaluation_indicator)) {
            eval_mask <- data[[comp$evaluation_indicator]] == 1
            fs_subset <- fs_subset[eval_mask, , drop = FALSE]
          }

          # Estimate loadings for each choice using linear probability model
          # (LPM is simpler and faster than running multinomial logit with factor scores)
          loading_init_by_choice <- list()
          for (choice in seq_len(comp$num_choices - 1)) {
            # Binary indicator: 1 if this choice, 0 otherwise
            choice_indicator <- as.numeric(outcome == (choice + 1))  # +1 because reference is 1

            # Run LPM: choice ~ X + factor_scores
            if (ncol(X) > 0) {
              X_full <- cbind(X, fs_subset)
            } else {
              X_full <- fs_subset
            }

            lpm_fit <- tryCatch(
              lm(choice_indicator ~ X_full - 1),
              error = function(e) NULL
            )

            if (!is.null(lpm_fit)) {
              all_coefs <- coef(lpm_fit)
              n_x <- ncol(X)
              loading_init_by_choice[[choice]] <- all_coefs[(n_x + 1):length(all_coefs)]
              names(loading_init_by_choice[[choice]]) <- NULL
            } else {
              loading_init_by_choice[[choice]] <- rep(0.5, n_free_loadings)
            }
          }
        }

        # Calculate number of weights needed
        n_weights <- (comp$num_choices - 1) * ncol(X)
        max_nwts <- max(1000, n_weights + 100)  # Allow enough for the model

        # Try to fit, fall back to zeros if too large or fails
        fit <- tryCatch(
          nnet::multinom(outcome ~ X - 1, trace = FALSE, MaxNWts = max_nwts),
          error = function(e) {
            if (verbose) message("  Note: multinom initialization failed, using zeros")
            NULL
          }
        )

        if (is.null(fit)) {
          # Fall back to zero initialization
          coefs_mat <- matrix(0, nrow = comp$num_choices - 1, ncol = ncol(X))
          rownames(coefs_mat) <- as.character(2:comp$num_choices)
          colnames(coefs_mat) <- colnames(X)
        } else {
          coefs_mat <- coef(fit)
          # Ensure coefs_mat is a matrix (multinom returns vector for binary case)
          if (!is.matrix(coefs_mat)) {
            coefs_mat <- matrix(coefs_mat, nrow = 1)
            colnames(coefs_mat) <- colnames(X)
          }
          # Check if we got the expected number of rows
          if (nrow(coefs_mat) < comp$num_choices - 1) {
            # Some choices might be missing - pad with zeros
            full_mat <- matrix(0, nrow = comp$num_choices - 1, ncol = ncol(X))
            colnames(full_mat) <- colnames(X)
            full_mat[seq_len(nrow(coefs_mat)), ] <- coefs_mat
            coefs_mat <- full_mat
          }
        }

        if (verbose) {
          n_fixed <- length(fixed_coefs)
          msg <- sprintf("  Multinomial logit: %d choices, %d covariates", comp$num_choices, ncol(X))
          if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
          if (use_factor_scores && n_free_loadings > 0) {
            msg <- paste0(msg, " (loadings estimated from factor scores)")
          }
          message(msg)
        }

        # Flatten parameters: for each choice (except reference), add covariates + loadings
        comp_params <- c()
        comp_param_names <- c()
        for (choice in seq_len(comp$num_choices - 1)) {
          # Get coefficients for this choice
          # Note: multinom names rows by actual choice values (2, 3, ...), not indices (1, 2, ...)
          if (comp$num_choices == 2) {
            choice_coefs <- coefs_mat
          } else if (is.matrix(coefs_mat) && nrow(coefs_mat) >= choice) {
            choice_coefs <- coefs_mat[choice, ]
          } else {
            # Fallback to zeros if coefs_mat doesn't have this choice
            choice_coefs <- rep(0, ncol(X))
            names(choice_coefs) <- colnames(X)
          }

          # Apply any fixed coefficient values for this choice
          choice_coefs <- apply_fixed_coefficients(choice_coefs, comp$covariates, fixed_coefs, choice = choice)

          # Get second-order loading initializations for this choice
          second_order <- get_second_order_loading_init(comp, choice = choice)

          # Get loading initialization for this choice
          if (!is.null(loading_init_by_choice) && choice <= length(loading_init_by_choice)) {
            choice_loading_init <- loading_init_by_choice[[choice]]
          } else {
            choice_loading_init <- rep(0.5, n_free_loadings)
          }

          # Add covariates, linear loadings, and second-order loadings
          comp_params <- c(comp_params, choice_coefs, choice_loading_init, second_order$values)

          # Build parameter names for this choice
          coef_names <- if (length(comp$covariates) > 0) {
            paste0(comp$name, "_c", choice, "_", comp$covariates)
          } else {
            character(0)
          }
          loading_names <- character(0)
          if (n_free_loadings > 0) {
            free_factor_idx <- which(is.na(comp$loading_normalization))
            loading_names <- paste0(comp$name, "_c", choice, "_loading_", free_factor_idx)
          }
          comp_param_names <- c(comp_param_names, coef_names, loading_names, second_order$names)
        }

        # Add type-specific intercepts if component uses types and n_types > 1
        # For multinomial logit, each non-reference choice gets type-specific intercepts
        if (isTRUE(comp$use_types) && n_types > 1L) {
          for (choice in seq_len(comp$num_choices - 1)) {
            type_intercepts <- numeric(0)
            type_intercept_names <- character(0)
            for (t in 2:n_types) {
              if (!is_type_intercept_fixed(comp, t, choice = choice)) {
                type_intercepts <- c(type_intercepts, 0.1 * (t - 1L))  # 0.1, 0.2, 0.3, ...
                type_intercept_names <- c(type_intercept_names, paste0(comp$name, "_c", choice, "_type_", t, "_intercept"))
              }
            }
            comp_params <- c(comp_params, type_intercepts)
            comp_param_names <- c(comp_param_names, type_intercept_names)
          }
        }
      }

    } else if (comp$model_type == "oprobit") {
      # Ordered probit - optionally include factor scores to estimate loadings
      loading_init <- rep(0.5, n_free_loadings)  # Default loading initialization

      # For ordered probit, intercept is not identified (absorbed into thresholds)
      # Remove intercept column if present
      X_no_int <- X
      intercept_col <- which(tolower(comp$covariates) == "intercept")
      if (length(intercept_col) > 0) {
        X_no_int <- X[, -intercept_col, drop = FALSE]
      }

      # Estimate loadings from factor scores using linear approximation
      if (use_factor_scores && n_free_loadings > 0) {
        # Get factor scores for FREE loadings only
        free_factor_idx <- which(is.na(comp$loading_normalization))
        fs_subset <- factor_scores[, free_factor_idx, drop = FALSE]
        if (!is.null(comp$evaluation_indicator)) {
          eval_mask <- data[[comp$evaluation_indicator]] == 1
          fs_subset <- fs_subset[eval_mask, , drop = FALSE]
        }

        # Use linear approximation: treat ordered outcome as continuous
        # This gives reasonable starting values for loadings
        outcome_numeric <- as.numeric(outcome)
        if (ncol(X_no_int) > 0) {
          X_full <- cbind(X_no_int, fs_subset)
        } else {
          X_full <- fs_subset
        }

        lpm_fit <- tryCatch(
          lm(outcome_numeric ~ X_full - 1),
          error = function(e) NULL
        )

        if (!is.null(lpm_fit)) {
          all_coefs <- coef(lpm_fit)
          n_x <- ncol(X_no_int)
          loading_init <- all_coefs[(n_x + 1):length(all_coefs)]
          names(loading_init) <- NULL
        }
      }

      # Fit model without intercept (for betas and thresholds)
      # Suppress MASS::polr warning about intercept (it's expected for ordered models)
      if (ncol(X_no_int) > 0) {
        fit <- suppressWarnings(
          MASS::polr(as.ordered(outcome) ~ X_no_int - 1, method = "probit")
        )
        coefs_no_int <- coef(fit)
      } else {
        # No covariates: fit intercept-only model to get thresholds
        fit <- suppressWarnings(
          MASS::polr(as.ordered(outcome) ~ 1, method = "probit")
        )
        coefs_no_int <- numeric(0)
      }
      thresholds_abs <- fit$zeta

      # Convert absolute thresholds to incremental parameterization
      # thresh1, thresh2, thresh3 -> thresh1, (thresh2-thresh1), (thresh3-thresh2)
      thresholds <- c(thresholds_abs[1])
      if (length(thresholds_abs) > 1) {
        for (i in 2:length(thresholds_abs)) {
          thresholds <- c(thresholds, thresholds_abs[i] - thresholds_abs[i-1])
        }
      }

      # Scale thresholds to match factor model variance at initialization
      # MASS::polr assumes Var(y*) = 1, but the factor model has:
      # Var(y*) = sum(lambda_k^2 * sigma_k^2) + 1
      # where lambda_k are loadings and sigma_k^2 are factor variances (all = 1.0 at init)

      # Calculate expected variance contribution from factors at initialization
      variance_from_factors <- 0.0

      if (!is.null(comp$loading_normalization)) {
        for (k in seq_along(comp$loading_normalization)) {
          loading_value <- comp$loading_normalization[k]

          if (is.na(loading_value)) {
            # Free loading: will be initialized to 0.5
            variance_from_factors <- variance_from_factors + (0.5^2) * 1.0
          } else if (abs(loading_value) > 1e-6) {
            # Fixed non-zero loading (e.g., 1.0)
            variance_from_factors <- variance_from_factors + (loading_value^2) * 1.0
          }
          # Fixed zero loadings contribute 0
        }
      }

      # Total variance at initialization: factor contributions + residual (1.0)
      total_variance <- variance_from_factors + 1.0

      # Scale thresholds if total variance differs from 1.0 (polr's assumption)
      if (abs(total_variance - 1.0) > 0.01) {
        scale_factor <- sqrt(total_variance)
        thresholds <- thresholds * scale_factor

        if (verbose) {
          message(sprintf("  Scaling thresholds by %.3f (factor model variance = %.2f at initialization)",
                          scale_factor, total_variance))
        }
      }

      # Build parameter vector: put intercept (0.0) first if it was in covariates
      # Note: For oprobit, intercept is theoretically absorbed into thresholds,
      # but we keep it in the parameter vector (as 0.0) to match C++ expectations
      if (length(intercept_col) > 0) {
        # Reconstruct in original order: intercept gets 0.0, others get fitted values
        coefs <- rep(0.0, ncol(X))
        coefs[-intercept_col] <- coefs_no_int
      } else {
        coefs <- coefs_no_int
      }

      # Apply any fixed coefficient values
      coefs <- apply_fixed_coefficients(coefs, comp$covariates, fixed_coefs)

      if (verbose) {
        n_fixed <- length(fixed_coefs)
        msg <- sprintf("  Ordered probit: %d covariates, %d thresholds (incremental form)", length(coefs), length(thresholds))
        if (n_fixed > 0) msg <- paste0(msg, sprintf(" (%d fixed)", n_fixed))
        if (use_factor_scores && n_free_loadings > 0) {
          msg <- paste0(msg, sprintf(" (loadings from factor scores: %s)",
                                     paste(sprintf("%.3f", loading_init), collapse = ", ")))
        }
        message(msg)
      }

      # Get second-order loading initializations
      second_order <- get_second_order_loading_init(comp)

      comp_params <- c(coefs, loading_init, second_order$values, thresholds)

      # Build parameter names
      coef_names <- character(0)
      if (length(coefs) > 0 && length(comp$covariates) > 0) {
        coef_names <- paste0(comp$name, "_", comp$covariates)
      }
      loading_names <- character(0)
      if (n_free_loadings > 0) {
        free_factor_idx <- which(is.na(comp$loading_normalization))
        loading_names <- paste0(comp$name, "_loading_", free_factor_idx)
      }
      n_thresholds <- length(thresholds)
      threshold_names <- paste0(comp$name, "_thresh_", seq_len(n_thresholds))
      comp_param_names <- c(coef_names, loading_names, second_order$names, threshold_names)

      # Add type-specific intercepts if component uses types and n_types > 1
      # For oprobit, type-specific intercepts shift all thresholds by a constant
      if (isTRUE(comp$use_types) && n_types > 1L) {
        type_intercepts <- numeric(n_types - 1L)
        type_intercept_names <- character(n_types - 1L)
        for (t in 2:n_types) {
          idx <- t - 1L
          type_intercept_names[idx] <- paste0(comp$name, "_type_", t, "_intercept")
          if (is_type_intercept_fixed(comp, t, choice = NULL)) {
            type_intercepts[idx] <- 0.0  # Fixed value
          } else {
            type_intercepts[idx] <- 0.1 * (t - 1L)  # 0.1, 0.2, 0.3, ...
          }
        }
        comp_params <- c(comp_params, type_intercepts)
        comp_param_names <- c(comp_param_names, type_intercept_names)
      }
    }

    # Add component parameters to overall parameter vector
    init_params <- c(init_params, comp_params)
    param_names <- c(param_names, comp_param_names)
  }
  }  # End of if (start_comp_idx <= n_total_comps)

  if (verbose) {
    message(sprintf("\nInitialized %d parameters total", length(init_params)))
  }

  # Determine factor_variance_fixed status
  # Must be a per-factor logical VECTOR (not scalar!) so that
  # setup_parameter_constraints can correctly fix/free each factor's variance.
  if (!is.null(model_system$previous_stage_info)) {
    # Compute per-factor identification from loadings (same logic as standard path)
    factor_variance_fixed_status <- rep(FALSE, n_factors)
    for (comp in model_system$components) {
      if (!is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          if (!is.na(comp$loading_normalization[k]) &&
              abs(comp$loading_normalization[k] - 1.0) < 1e-6) {
            factor_variance_fixed_status[k] <- TRUE
          }
        }
      }
    }
    # If factor_var_k is explicitly in free_params, mark factor k as identified
    # so its variance stays free regardless of loading normalization.
    .fp <- model_system$previous_stage_info$free_param_names
    if (!is.null(.fp)) {
      for (fpn in .fp) {
        m <- regmatches(fpn, regexec("^(?:mix\\d+_)?factor_var_(\\d+)$", fpn))[[1]]
        if (length(m) >= 2) {
          k <- as.integer(m[2])
          if (k <= n_factors) factor_variance_fixed_status[k] <- TRUE
        }
      }
    }
  } else {
    # Standard case: use the computed status (vector for multifactor models)
    factor_variance_fixed_status <- factor_variance_fixed
  }

  # Assign names to parameter vector
  names(init_params) <- param_names

  # Compute param_fixed using setup_parameter_constraints
  # Apply user-fixed factor-distribution parameters from fix_factor_param()
  # BEFORE setup_parameter_constraints. This ensures init_params carries the
  # user-supplied value at the fixed position so the C++ side
  # (SetParameterConstraintsWithValues) starts there. setup_parameter_constraints
  # will then mark these positions as fixed independently; the values match.
  if (!is.null(model_system$factor$fixed_params) &&
      length(model_system$factor$fixed_params) > 0L) {
    fp <- model_system$factor$fixed_params
    fp_idx <- match(names(fp), param_names)
    valid <- !is.na(fp_idx)
    if (any(valid)) {
      init_params[fp_idx[valid]] <- unname(fp)[valid]
    }
  }

  # This is needed for gradient/Hessian checking in tests
  param_metadata <- build_parameter_metadata(model_system)
  param_constraints <- setup_parameter_constraints(
    model_system, init_params, param_metadata,
    factor_variance_fixed_status, verbose = FALSE
  )

  list(
    init_params = init_params,
    param_names = param_names,
    factor_variance_fixed = factor_variance_fixed_status,
    param_fixed = param_constraints$param_fixed
  )
}


#' Get second-order loading initialization values and names
#'
#' Helper function to generate initial values and parameter names for
#' quadratic and interaction factor loadings.
#'
#' @param comp Model component object
#' @param choice Integer or NULL. For multinomial logit, which choice.
#' @return List with values and names for second-order loadings
#' @keywords internal
get_second_order_loading_init <- function(comp, choice = NULL) {
  values <- numeric(0)
  names_vec <- character(0)

  k <- comp$k
  comp_name <- comp$name
  choice_suffix <- if (!is.null(choice)) paste0("_c", choice) else ""

  # For dynamic models, skip the outcome factor
  is_dynamic <- isTRUE(comp$is_dynamic)
  outcome_factor <- if (is_dynamic) comp$outcome_factor else -1L

  # Quadratic loadings: one per factor (skip outcome factor for dynamic models)
  if (!is.null(comp$n_quadratic_loadings) && comp$n_quadratic_loadings > 0) {
    values <- c(values, rep(0.1, comp$n_quadratic_loadings))
    quad_names <- character(0)
    for (fac_idx in seq_len(k)) {
      if (is_dynamic && fac_idx == outcome_factor) next
      quad_names <- c(quad_names, paste0(comp_name, choice_suffix, "_loading_quad_", fac_idx))
    }
    names_vec <- c(names_vec, quad_names)
  }

  # Interaction loadings: one per unique pair j < k (skip pairs involving outcome factor for dynamic models)
  if (!is.null(comp$n_interaction_loadings) && comp$n_interaction_loadings > 0) {
    values <- c(values, rep(0.1, comp$n_interaction_loadings))
    inter_names <- character(0)
    for (j in seq_len(k - 1)) {
      if (is_dynamic && j == outcome_factor) next
      for (kk in (j + 1):k) {
        if (is_dynamic && kk == outcome_factor) next
        inter_names <- c(inter_names, paste0(comp_name, choice_suffix, "_loading_inter_", j, "_", kk))
      }
    }
    names_vec <- c(names_vec, inter_names)
  }

  list(values = values, names = names_vec)
}


#' Apply fixed coefficient values to estimated coefficients
#'
#' Helper function to replace estimated coefficients with fixed values
#' where specified in the component's fixed_coefficients list.
#'
#' @param coefs Named numeric vector of estimated coefficients
#' @param covariates Character vector of covariate names (in order)
#' @param fixed_coefficients List of fixed coefficient constraints from component
#' @param choice Integer or NULL. For multinomial logit, which choice these coefs belong to.
#' @return Numeric vector with fixed values applied
#' @keywords internal
apply_fixed_coefficients <- function(coefs, covariates, fixed_coefficients, choice = NULL) {
  if (length(fixed_coefficients) == 0) {
    return(coefs)
  }

  for (fc in fixed_coefficients) {
    # Check if this constraint applies (matching choice for mlogit)
    if (!identical(fc$choice, choice)) {
      next
    }

    # Find the position of this covariate
    pos <- match(fc$covariate, covariates)
    if (!is.na(pos)) {
      coefs[pos] <- fc$value
    }
  }

  return(coefs)
}
