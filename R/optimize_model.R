# ---- Parallel cleanup utilities ----

#' Clean up orphaned parallel worker processes
#'
#' When using parallel estimation with Ctrl-C interrupts, worker processes may
#' continue running after the main process exits. This function finds and
#' terminates any orphaned R worker processes started by the parallel package.
#'
#' @param signal Signal to send (default "TERM" for graceful shutdown, use "KILL" for force)
#' @param verbose Whether to print messages (default TRUE)
#' @param list_only If TRUE, only list potential workers without killing them (default FALSE)
#' @return Invisible NULL (or vector of PIDs if list_only = TRUE)
#' @export
#'
#' @examples
#' # Safe to run: list potential orphaned parallel workers without killing them.
#' cleanup_parallel_workers(list_only = TRUE, verbose = FALSE)
#'
#' \donttest{
#' # After interrupting a parallel job with Ctrl-C, terminate orphaned workers:
#' cleanup_parallel_workers()
#'
#' # Force kill if graceful shutdown doesn't work:
#' cleanup_parallel_workers(signal = "KILL")
#' }
cleanup_parallel_workers <- function(signal = "TERM", verbose = TRUE, list_only = FALSE) {
  # Get current R process ID to avoid killing ourselves
  my_pid <- Sys.getpid()

  if (.Platform$OS.type == "unix") {
    # On Unix/Linux/macOS, find R processes that look like workers
    # Try multiple patterns to catch different worker types

    # Pattern 1: PSOCK workers with workRSOCK (most reliable - unique to parallel workers)
    # Pattern 2: R workers with --no-echo --no-restore (common parallel worker flags)
    # Use grep -v grep to exclude grep from matching itself
    patterns <- c(
      "ps aux | grep 'workRSOCK' | grep -v grep",
      "ps aux | grep 'slaveRSOCK' | grep -v grep",
      "ps aux | grep '/R.*--no-echo.*--no-restore.*parallel' | grep -v grep"
    )

    all_pids <- integer(0)
    all_info <- character(0)

    for (pattern in patterns) {
      tryCatch({
        result <- suppressWarnings(system(pattern, intern = TRUE, ignore.stderr = TRUE))
        if (length(result) > 0) {
          # Extract PIDs (second field) and full line for info
          for (line in result) {
            # Skip lines that are clearly not R worker processes
            if (grepl("grep", line, fixed = TRUE)) next
            fields <- strsplit(trimws(line), "\\s+")[[1]]
            if (length(fields) >= 2) {
              pid <- as.integer(fields[2])
              if (!is.na(pid) && !(pid %in% all_pids)) {
                all_pids <- c(all_pids, pid)
                all_info <- c(all_info, line)
              }
            }
          }
        }
      }, error = function(e) NULL)
    }

    # Remove our own PID if somehow included
    keep <- !(all_pids %in% my_pid)
    all_pids <- all_pids[keep]
    all_info <- all_info[keep]

    if (length(all_pids) == 0) {
      if (verbose) message("No orphaned parallel workers found.")
      return(invisible(NULL))
    }

    if (verbose || list_only) {
      message(sprintf("Found %d potential worker process(es):", length(all_pids)))
      for (i in seq_along(all_pids)) {
        # Truncate long lines for display
        info <- if (nchar(all_info[i]) > 100) {
          paste0(substr(all_info[i], 1, 97), "...")
        } else {
          all_info[i]
        }
        message(sprintf("  PID %d: %s", all_pids[i], info))
      }
    }

    if (list_only) {
      return(invisible(all_pids))
    }

    # Send signal to each worker
    for (pid in all_pids) {
      tryCatch({
        tools::pskill(pid, signal = if (signal == "KILL") tools::SIGKILL else tools::SIGTERM)
        if (verbose) message(sprintf("  Sent %s to PID %d", signal, pid))
      }, error = function(e) {
        if (verbose) message(sprintf("  Failed to signal PID %d: %s", pid, e$message))
      })
    }

  } else {
    # On Windows, use taskkill
    # Worker processes have "Rscript" in their name
    if (verbose) {
      message("On Windows, you can manually kill orphaned R workers with:")
      message("  taskkill /F /IM Rscript.exe")
    }
  }

  invisible(NULL)
}

# ---- Internal helper functions ----

# Build parameter metadata from model system
build_parameter_metadata <- function(model_system) {
  n_factors <- model_system$factor$n_factors
  n_mixtures <- model_system$factor$n_mixtures
  if (is.null(n_mixtures)) n_mixtures <- 1L

  # Initialize vectors
  param_names <- character(0)
  param_types <- character(0)  # "factor_var", "factor_corr", "intercept", "beta", "loading", "sigma", "cutpoint"
  component_id <- integer(0)

  # Add factor variances (for SE models, only input factors have variance params)
  factor_structure <- model_system$factor$factor_structure
  if (is.null(factor_structure)) factor_structure <- "independent"

  # Determine number of factors with variance parameters
  if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
    n_factors_for_mixture <- n_factors - 1L  # Only input factors
  } else {
    n_factors_for_mixture <- n_factors
  }

  if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
    # SE models: only first (n_factors - 1) have variance parameters
    n_input_factors <- n_factors - 1
    for (imix in seq_len(n_mixtures)) {
      for (k in seq_len(n_input_factors)) {
        if (n_mixtures == 1) {
          param_names <- c(param_names, sprintf("factor_var_%d", k))
        } else {
          param_names <- c(param_names, sprintf("mix%d_factor_var_%d", imix, k))
        }
        param_types <- c(param_types, "factor_var")
        component_id <- c(component_id, 0)  # 0 = factor model
      }
    }

    # Mixture means (for non-reference mixtures only)
    if (n_mixtures > 1) {
      for (imix in seq_len(n_mixtures - 1)) {
        for (k in seq_len(n_input_factors)) {
          param_names <- c(param_names, sprintf("mix%d_factor_mean_%d", imix, k))
          param_types <- c(param_types, "factor_mean_mix")
          component_id <- c(component_id, 0)
        }
      }
      # Mixture log-weights (for non-reference mixtures)
      for (imix in seq_len(n_mixtures - 1)) {
        param_names <- c(param_names, sprintf("mix%d_logweight", imix))
        param_types <- c(param_types, "mix_logweight")
        component_id <- c(component_id, 0)
      }
    }

    # SE intercept
    param_names <- c(param_names, "se_intercept")
    param_types <- c(param_types, "se_intercept")
    component_id <- c(component_id, 0)

    # SE linear coefficients
    for (j in seq_len(n_input_factors)) {
      param_names <- c(param_names, sprintf("se_linear_%d", j))
      param_types <- c(param_types, "se_linear")
      component_id <- c(component_id, 0)
    }

    # SE quadratic coefficients (for SE_quadratic only)
    if (factor_structure == "SE_quadratic") {
      for (j in seq_len(n_input_factors)) {
        param_names <- c(param_names, sprintf("se_quadratic_%d", j))
        param_types <- c(param_types, "se_quadratic")
        component_id <- c(component_id, 0)
      }
    }

    # Type-specific SE intercepts (only when n_types > 1)
    # Order must match C++ layout in FactorModel.h::GetSETypeInterceptIndex:
    # between (quadratic) coefs and se_residual_var.
    .n_types_se <- model_system$factor$n_types
    if (!is.null(.n_types_se) && .n_types_se > 1L) {
      for (t in 2:.n_types_se) {
        param_names <- c(param_names, sprintf("se_intercept_type_%d", t))
        param_types <- c(param_types, "se_intercept_type")
        component_id <- c(component_id, 0)
      }
    }

    # SE residual variance
    param_names <- c(param_names, "se_residual_var")
    param_types <- c(param_types, "se_residual_var")
    component_id <- c(component_id, 0)

  } else {
    # Standard case: all n_factors have variance parameters
    for (imix in seq_len(n_mixtures)) {
      for (k in seq_len(n_factors)) {
        if (n_mixtures == 1) {
          param_names <- c(param_names, sprintf("factor_var_%d", k))
        } else {
          param_names <- c(param_names, sprintf("mix%d_factor_var_%d", imix, k))
        }
        param_types <- c(param_types, "factor_var")
        component_id <- c(component_id, 0)  # 0 = factor model
      }
    }

    # Mixture means (for non-reference mixtures only)
    if (n_mixtures > 1) {
      for (imix in seq_len(n_mixtures - 1)) {
        for (k in seq_len(n_factors)) {
          param_names <- c(param_names, sprintf("mix%d_factor_mean_%d", imix, k))
          param_types <- c(param_types, "factor_mean_mix")
          component_id <- c(component_id, 0)
        }
      }
      # Mixture log-weights (for non-reference mixtures)
      for (imix in seq_len(n_mixtures - 1)) {
        param_names <- c(param_names, sprintf("mix%d_logweight", imix))
        param_types <- c(param_types, "mix_logweight")
        component_id <- c(component_id, 0)
      }
    }
  }

  # Add factor correlation if correlation = TRUE and n_factors = 2
  if (isTRUE(model_system$factor$correlation) && n_factors == 2) {
    param_names <- c(param_names, "factor_corr_1_2")
    param_types <- c(param_types, "factor_corr")
    component_id <- c(component_id, 0)  # 0 = factor model
  }

  # Add factor mean covariate parameters if specified
  factor_covariates <- model_system$factor$factor_covariates
  if (!is.null(factor_covariates) && length(factor_covariates) > 0) {
    # Determine how many factors get mean covariates
    if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
      n_factors_with_mean <- n_factors - 1L  # Only input factors
    } else {
      n_factors_with_mean <- n_factors
    }
    for (k in seq_len(n_factors_with_mean)) {
      for (cov_name in factor_covariates) {
        param_names <- c(param_names, sprintf("factor_mean_%d_%s", k, cov_name))
        param_types <- c(param_types, "factor_mean")
        component_id <- c(component_id, 0)  # 0 = factor model
      }
    }
  }

  # Add SE covariate parameters if specified (for SE_linear/SE_quadratic)
  se_covariates <- model_system$factor$se_covariates
  if (!is.null(se_covariates) && length(se_covariates) > 0) {
    for (cov_name in se_covariates) {
      param_names <- c(param_names, sprintf("se_cov_%s", cov_name))
      param_types <- c(param_types, "se_covariate")
      component_id <- c(component_id, 0)  # 0 = factor model
    }
  }

  # Add type model parameters if n_types > 1 and at least one component uses types
  # Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
  # (n_types - 1) intercepts + (n_types - 1) * n_factors loadings (type 1 is reference)
  n_types <- model_system$factor$n_types
  any_uses_types <- any(sapply(model_system$components, function(c) isTRUE(c$use_types)))
  # SE_linear / SE_quadratic with ntyp > 1 implies types at the structural level
  # (se_intercept_type_{t}), so the type probability model is needed here too.
  .fs_opt <- model_system$factor$factor_structure
  if (!is.null(.fs_opt) && .fs_opt %in% c("SE_linear", "SE_quadratic") &&
      !is.null(n_types) && n_types > 1L) {
    any_uses_types <- TRUE
  }
  if (!is.null(n_types) && n_types > 1L && any_uses_types) {
    # Type probability intercepts
    for (t in 2:n_types) {
      param_names <- c(param_names, sprintf("typeprob_%d_intercept", t))
      param_types <- c(param_types, "typeprob_intercept")
      component_id <- c(component_id, 0)  # 0 = factor model
    }
    # Type probability loadings
    for (t in 2:n_types) {
      for (k in seq_len(n_factors)) {
        param_names <- c(param_names, sprintf("type_%d_loading_%d", t, k))
        param_types <- c(param_types, "type_loading")
        component_id <- c(component_id, 0)  # 0 = factor model
      }
    }
  }

  # Add component parameters
  for (i in seq_along(model_system$components)) {
    comp <- model_system$components[[i]]
    comp_name <- comp$name

    # Special handling for multinomial logit (mlogit) with >2 choices
    # Parameters are organized per choice (except reference choice 0)
    if (comp$model_type == "logit" && !is.null(comp$num_choices) && comp$num_choices > 2) {
      n_alternatives <- comp$num_choices - 1  # Exclude reference category

      # Check if intercept is already in covariates (as "intercept" or "constant")
      has_intercept_in_covariates <- !is.null(comp$covariates) &&
        any(comp$covariates %in% c("intercept", "constant"))

      for (alt in seq_len(n_alternatives)) {
        # Intercept for this alternative - only add if not already in covariates
        # C++ counts all covariates including "intercept"/"constant" as regressors
        if (comp$intercept && !has_intercept_in_covariates) {
          param_names <- c(param_names, sprintf("%s_intercept_alt%d", comp_name, alt))
          param_types <- c(param_types, "intercept")
          component_id <- c(component_id, i)
        }

        # Covariate coefficients for this alternative
        # Include ALL covariates (even fixed ones) to match init_params vector
        # Note: "intercept" and "constant" are counted as regular betas to match C++
        if (!is.null(comp$covariates) && length(comp$covariates) > 0) {
          for (cov in comp$covariates) {
            param_names <- c(param_names, sprintf("%s_beta_%s_alt%d", comp_name, cov, alt))
            param_types <- c(param_types, "beta")
            component_id <- c(component_id, i)
          }
        }

        # Factor loadings for this alternative
        if (is.null(comp$loading_normalization)) {
          for (k in seq_len(n_factors)) {
            param_names <- c(param_names, sprintf("%s_loading_%d_alt%d", comp_name, k, alt))
            param_types <- c(param_types, "loading")
            component_id <- c(component_id, i)
          }
        } else {
          for (k in seq_len(n_factors)) {
            if (is.na(comp$loading_normalization[k]) ||
                abs(comp$loading_normalization[k]) < 1e-10) {
              param_names <- c(param_names, sprintf("%s_loading_%d_alt%d", comp_name, k, alt))
              param_types <- c(param_types, "loading")
              component_id <- c(component_id, i)
            }
          }
        }

        # Quadratic factor loadings for this alternative (if factor_spec includes quadratic)
        if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("quadratic", "full")) {
          for (k in seq_len(n_factors)) {
            param_names <- c(param_names, sprintf("%s_loading_quad_%d_alt%d", comp_name, k, alt))
            param_types <- c(param_types, "loading_quad")
            component_id <- c(component_id, i)
          }
        }

        # Interaction factor loadings for this alternative (if factor_spec includes interactions)
        if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("interactions", "full") && n_factors >= 2) {
          for (j in seq_len(n_factors - 1)) {
            for (kk in (j + 1):n_factors) {
              param_names <- c(param_names, sprintf("%s_loading_inter_%d_%d_alt%d", comp_name, j, kk, alt))
              param_types <- c(param_types, "loading_inter")
              component_id <- c(component_id, i)
            }
          }
        }
      }
    } else {
      # Standard handling for all other model types
      # Note: For linear models with intercept=TRUE, "intercept" is automatically
      # added to covariates in define_model_component(), so no special handling needed here.

      # Covariate coefficients - naming must match initialize_parameters.R: comp_name_covariate
      # Include ALL covariates (even fixed ones) to match init_params vector
      if (!is.null(comp$covariates) && length(comp$covariates) > 0) {
        for (cov in comp$covariates) {
          param_names <- c(param_names, sprintf("%s_%s", comp_name, cov))
          # Mark intercept type separately for constraint handling
          # Both "intercept" and "constant" are treated as intercept
          if (cov %in% c("intercept", "constant")) {
            param_types <- c(param_types, "intercept")
          } else {
            param_types <- c(param_types, "beta")
          }
          component_id <- c(component_id, i)
        }
      }

      # Factor loadings
      if (is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          param_names <- c(param_names, sprintf("%s_loading_%d", comp_name, k))
          param_types <- c(param_types, "loading")
          component_id <- c(component_id, i)
        }
      } else {
        for (k in seq_len(n_factors)) {
          if (is.na(comp$loading_normalization[k])) {
            param_names <- c(param_names, sprintf("%s_loading_%d", comp_name, k))
            param_types <- c(param_types, "loading")
            component_id <- c(component_id, i)
          }
        }
      }

      # Quadratic factor loadings (if factor_spec includes quadratic)
      if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("quadratic", "full")) {
        for (k in seq_len(n_factors)) {
          param_names <- c(param_names, sprintf("%s_loading_quad_%d", comp_name, k))
          param_types <- c(param_types, "loading_quad")
          component_id <- c(component_id, i)
        }
      }

      # Interaction factor loadings (if factor_spec includes interactions)
      if (!is.null(comp$factor_spec) && comp$factor_spec %in% c("interactions", "full") && n_factors >= 2) {
        for (j in seq_len(n_factors - 1)) {
          for (kk in (j + 1):n_factors) {
            param_names <- c(param_names, sprintf("%s_loading_inter_%d_%d", comp_name, j, kk))
            param_types <- c(param_types, "loading_inter")
            component_id <- c(component_id, i)
          }
        }
      }
    }

    # Residual variance (only for linear models - oprobit has normalized variance)
    if (comp$model_type == "linear") {
      param_names <- c(param_names, sprintf("%s_sigma", comp_name))
      param_types <- c(param_types, "sigma")
      component_id <- c(component_id, i)
    }

    # Thresholds for ordered probit (num_choices - 1 thresholds for identification)
    if (comp$model_type == "oprobit" && !is.null(comp$num_choices) && comp$num_choices > 1) {
      n_thresholds <- comp$num_choices - 1
      for (j in seq_len(n_thresholds)) {
        param_names <- c(param_names, sprintf("%s_thresh_%d", comp_name, j))
        param_types <- c(param_types, "cutpoint")
        component_id <- c(component_id, i)
      }
    }

    # Type-specific intercepts for this component (only if use_types = TRUE and n_types > 1)
    # Added after all other component parameters to match initialize_parameters.R ordering
    if (isTRUE(comp$use_types) && !is.null(n_types) && n_types > 1L) {
      # For multinomial logit, each non-reference choice gets type-specific intercepts
      if (comp$model_type == "logit" && !is.null(comp$num_choices) && comp$num_choices > 2) {
        for (choice in seq_len(comp$num_choices - 1)) {
          for (t in 2:n_types) {
            param_names <- c(param_names, sprintf("%s_c%d_type_%d_intercept", comp_name, choice, t))
            param_types <- c(param_types, "type_intercept")
            component_id <- c(component_id, i)
          }
        }
      } else {
        # Standard case: one type intercept per type
        for (t in 2:n_types) {
          param_names <- c(param_names, sprintf("%s_type_%d_intercept", comp_name, t))
          param_types <- c(param_types, "type_intercept")
          component_id <- c(component_id, i)
        }
      }
    }
  }

  return(list(
    names = param_names,
    types = param_types,
    component_id = component_id,
    n_params = length(param_names)
  ))
}

# Setup parameter bounds and identify fixed parameters
setup_parameter_constraints <- function(model_system, init_params, param_metadata, factor_variance_fixed = NULL, verbose = FALSE) {
  n_params <- length(init_params)
  n_factors <- model_system$factor$n_factors

  lower_bounds <- rep(-Inf, n_params)
  upper_bounds <- rep(Inf, n_params)

  # Identify which factor variances are identified
  # If factor_variance_fixed was passed from initialize_parameters, use it
  # Otherwise, compute it by checking for fixed non-zero loadings
  if (!is.null(factor_variance_fixed)) {
    # Use the passed value: factor_variance_fixed=TRUE means variance is estimated (identified by loading)
    # So factor_variance_identified = factor_variance_fixed
    factor_variance_identified <- factor_variance_fixed
    if (verbose) {
      message(sprintf("Factor variance identification (from initialize_parameters): %s",
                      paste(factor_variance_identified, collapse=", ")))
    }
  } else {
    # Fallback: compute by checking for fixed non-zero loadings
    # A factor is identified if at least one component has a fixed non-zero loading on that factor
    factor_variance_identified <- rep(FALSE, n_factors)
    for (comp in model_system$components) {
      if (!is.null(comp$loading_normalization)) {
        for (k in seq_len(n_factors)) {
          if (!is.na(comp$loading_normalization[k]) &&
              abs(comp$loading_normalization[k]) > 1e-6) {
            factor_variance_identified[k] <- TRUE
          }
        }
      }
    }
    if (verbose) {
      message(sprintf("Factor variance identification (computed fallback): %s",
                      paste(factor_variance_identified, collapse=", ")))
    }
  }

  param_fixed <- rep(FALSE, n_params)

  # Handle previous_stage: mark previous-stage parameters as fixed (by name, not position)
  if (!is.null(model_system$previous_stage_info)) {
    fixed_param_names <- model_system$previous_stage_info$fixed_param_names
    fixed_param_values <- model_system$previous_stage_info$fixed_param_values

    if (length(fixed_param_names) > 0) {
      # Match fixed parameter names to indices in current model
      fixed_indices <- match(fixed_param_names, param_metadata$names)

      # Remove NA matches (parameters not in current model)
      valid_matches <- !is.na(fixed_indices)
      if (any(!valid_matches) && verbose) {
        message(sprintf("Note: %d previous_stage parameters not in current model (expected for new components)",
                        sum(!valid_matches)))
      }

      fixed_indices <- fixed_indices[valid_matches]
      fixed_values <- fixed_param_values[valid_matches]

      if (length(fixed_indices) > 0) {
        param_fixed[fixed_indices] <- TRUE
        lower_bounds[fixed_indices] <- fixed_values
        upper_bounds[fixed_indices] <- fixed_values

        if (verbose) {
          message(sprintf("Fixed %d parameters from previous_stage", length(fixed_indices)))
          if (!is.null(model_system$previous_stage_info$free_param_names)) {
            message(sprintf("Free parameters from previous_stage: %s",
                            paste(model_system$previous_stage_info$free_param_names, collapse = ", ")))
          }
        }
      }
    }
  }

  # Track cutpoint indices per component to identify incremental thresholds
  cutpoint_counter <- list()

  for (i in seq_len(n_params)) {
    param_type <- param_metadata$types[i]
    comp_id <- param_metadata$component_id[i]

    # Fix non-identified factor variances
    if (param_type == "factor_var") {
      # For factor variances, check identification
      # Extract factor index from parameter name (e.g., "factor_var_1" or "mix1_factor_var_1")
      param_name <- param_metadata$names[i]
      factor_idx_match <- regmatches(param_name, regexpr("factor_var_(\\d+)$", param_name))
      if (length(factor_idx_match) > 0) {
        factor_idx <- as.integer(sub("factor_var_", "", factor_idx_match))
        if (factor_idx <= length(factor_variance_identified)) {
          if (verbose) {
            message(sprintf("  Constraint check: param %d (%s) -> factor_idx=%d, identified=%s",
                            i, param_name, factor_idx, factor_variance_identified[factor_idx]))
          }
          if (!factor_variance_identified[factor_idx]) {
            param_fixed[i] <- TRUE
            lower_bounds[i] <- init_params[i]
            upper_bounds[i] <- init_params[i]
            if (verbose) {
              message(sprintf("    -> FIXED at %.6f", init_params[i]))
            }
          } else {
            # Set lower bound for free factor variances to prevent numerical issues
            # (division by sqrt(factor_var) in gradient/Hessian chain rule)
            lower_bounds[i] <- 0.01
          }
        }
      }
    }

    # Mixture log-weights: no special constraints (softmax handles the range)
    # Mixture means: no special constraints

    # Factor mean mixture parameters: set reasonable bounds
    if (param_type == "factor_mean_mix") {
      # Means can be any real number, but large values may cause numerical issues
      # No bounds needed - the constraint E[f]=0 is handled by the last mixture's mean
    }

    # Set bounds for factor correlation parameters
    if (param_type == "factor_corr") {
      # Correlation must be between -1 and 1
      lower_bounds[i] <- -0.99
      upper_bounds[i] <- 0.99
      if (verbose) {
        message(sprintf("  Constraint: param %d (factor_corr) -> bounds = [-0.99, 0.99]", i))
      }
    }

    # Set lower bound for sigma parameters
    if (param_type == "sigma") {
      lower_bounds[i] <- 0.01
    }

    # Set lower bound for SE residual variance. Without a positive lower bound,
    # the optimizer can walk into negative territory; historically this was
    # "patched" in C++ by computing sigma_eps = sqrt(|se_residual_var|), which
    # makes the likelihood at +x and -x identical but breaks the reported
    # estimate (negative value shown to user) and the Hessian-based SE at the
    # boundary. Enforcing a positive lower bound removes both issues.
    if (param_type == "se_residual_var") {
      lower_bounds[i] <- 0.01
    }

    # Set lower bounds for ordered probit cutpoints (incremental parameterization)
    # The first cutpoint is unrestricted (can be any value)
    # Subsequent cutpoints are increments and must be positive to ensure ordering
    if (param_type == "cutpoint") {
      comp_key <- as.character(comp_id)
      if (is.null(cutpoint_counter[[comp_key]])) {
        cutpoint_counter[[comp_key]] <- 1
        # First cutpoint is unrestricted
      } else {
        cutpoint_counter[[comp_key]] <- cutpoint_counter[[comp_key]] + 1
        # Subsequent cutpoints are increments, must be positive
        lower_bounds[i] <- 0.01
      }
    }

    # For SE_linear / SE_quadratic factor structures, type probabilities must be a
    # function of the INPUT factors only — the outcome factor is a deterministic
    # function of the inputs + residual + type, so letting type loadings depend on
    # it would create a circular dependency. We enforce this by fixing the type
    # loading on the outcome factor (factor index = n_factors) to 0 for each
    # non-reference type, and erroring if the user explicitly set a non-zero init.
    if (param_type == "type_loading") {
      fs <- model_system$factor$factor_structure
      if (!is.null(fs) && fs %in% c("SE_linear", "SE_quadratic")) {
        pname <- param_metadata$names[i]
        m <- regmatches(pname, regexec("^type_([0-9]+)_loading_([0-9]+)$", pname))[[1]]
        if (length(m) >= 3) {
          k_idx <- as.integer(m[3])
          if (k_idx == n_factors) {
            # Threshold: tolerate finite-difference perturbations (~1.5e-8) and
            # other numerical noise, but reject user-supplied non-zero values.
            if (abs(init_params[i]) > 1e-6) {
              stop(sprintf(
                paste0("Type loading on the outcome factor is not allowed for ",
                       "SE_linear / SE_quadratic factor structures (parameter %s = %g). ",
                       "Type probabilities must depend only on input factors; the outcome ",
                       "factor is a deterministic function of type, so a type loading on ",
                       "it would create a circular dependency. Leave the initial value at 0 ",
                       "and it will be fixed automatically."),
                pname, init_params[i]))
            }
            param_fixed[i] <- TRUE
            lower_bounds[i] <- 0.0
            upper_bounds[i] <- 0.0
            if (verbose) {
              message(sprintf(
                "  Fixed type loading on outcome factor for SE model: param %d (%s) at 0",
                i, pname))
            }
          }
        }
      }
    }

    # Fix oprobit intercepts at 0 (absorbed into thresholds)
    if (param_type == "intercept") {
      comp <- model_system$components[[comp_id]]
      if (!is.null(comp) && comp$model_type == "oprobit") {
        param_fixed[i] <- TRUE
        lower_bounds[i] <- 0.0
        upper_bounds[i] <- 0.0
        if (verbose) {
          message(sprintf("  Fixed oprobit intercept: param %d (%s) at 0",
                          i, param_metadata$names[i]))
        }
      }
    }

    # Fix type intercepts that were marked as fixed via fix_type_intercepts()
    if (param_type == "type_intercept") {
      comp <- model_system$components[[comp_id]]
      if (!is.null(comp) && !is.null(comp$fixed_type_intercepts) && length(comp$fixed_type_intercepts) > 0) {
        # Extract type number from parameter name (e.g., "Y_type_2_intercept" -> 2)
        param_name <- param_metadata$names[i]
        type_match <- regmatches(param_name, regexec("_type_([0-9]+)_intercept$", param_name))[[1]]
        if (length(type_match) >= 2) {
          type_num <- as.integer(type_match[2])
          if (is_type_intercept_fixed(comp, type_num, choice = NULL)) {
            param_fixed[i] <- TRUE
            lower_bounds[i] <- 0.0
            upper_bounds[i] <- 0.0
            if (verbose) {
              message(sprintf("  Fixed type intercept: param %d (%s) at 0",
                              i, param_name))
            }
          }
        }
      }
    }

    # Fix coefficients that were marked as fixed via fix_coefficient()
    if (param_type == "beta" || param_type == "intercept") {
      comp <- model_system$components[[comp_id]]
      if (!is.null(comp) && !is.null(comp$fixed_coefficients) && length(comp$fixed_coefficients) > 0) {
        param_name <- param_metadata$names[i]

        # Parse covariate name from parameter name
        # For standard models: compname_covariate (e.g., "Y_x1")
        # For multinomial logit: compname_beta_covariate_altN (e.g., "Y_beta_x1_alt1")
        comp_name <- comp$name
        covariate <- NULL
        choice <- NULL

        if (comp$model_type == "logit" && !is.null(comp$num_choices) && comp$num_choices > 2) {
          # Multinomial logit: parse "compname_beta_covariate_altN" or "compname_intercept_altN"
          if (param_type == "intercept") {
            # Format: compname_intercept_altN
            pattern <- sprintf("^%s_intercept_alt([0-9]+)$", comp_name)
            match <- regmatches(param_name, regexec(pattern, param_name))[[1]]
            if (length(match) >= 2) {
              covariate <- "intercept"
              choice <- as.integer(match[2])
            }
          } else {
            # Format: compname_beta_covariate_altN
            pattern <- sprintf("^%s_beta_(.+)_alt([0-9]+)$", comp_name)
            match <- regmatches(param_name, regexec(pattern, param_name))[[1]]
            if (length(match) >= 3) {
              covariate <- match[2]
              choice <- as.integer(match[3])
            }
          }
        } else {
          # Standard model: parse "compname_covariate"
          pattern <- sprintf("^%s_(.+)$", comp_name)
          match <- regmatches(param_name, regexec(pattern, param_name))[[1]]
          if (length(match) >= 2) {
            covariate <- match[2]
          }
        }

        if (!is.null(covariate) && is_coefficient_fixed(comp, covariate, choice)) {
          fixed_value <- get_fixed_coefficient_value(comp, covariate, choice)
          param_fixed[i] <- TRUE
          lower_bounds[i] <- fixed_value
          upper_bounds[i] <- fixed_value
          if (verbose) {
            message(sprintf("  Fixed coefficient: param %d (%s) at %.6f",
                            i, param_name, fixed_value))
          }
        }
      }
    }
  }

  # Handle equality constraints
  # For each constraint group, mark derived params (2nd, 3rd, ...) as fixed
  # and store the mapping to their primary parameter
  equality_mapping <- list()  # Maps derived param index -> primary param index

  if (!is.null(model_system$equality_constraints)) {
    for (constraint in model_system$equality_constraints) {
      # Find parameter indices for each name in the constraint
      param_indices <- match(constraint, param_metadata$names)

      # Check that all parameter names were found
      if (any(is.na(param_indices))) {
        missing <- constraint[is.na(param_indices)]
        stop(sprintf("Equality constraint references unknown parameter(s): %s\nAvailable parameters: %s",
                     paste(missing, collapse = ", "),
                     paste(param_metadata$names, collapse = ", ")))
      }

      # First parameter is the primary, rest are derived
      primary_idx <- param_indices[1]

      # Check that primary is not already fixed or derived
      if (param_fixed[primary_idx]) {
        stop(sprintf("Equality constraint primary parameter '%s' is already fixed", constraint[1]))
      }

      for (i in 2:length(param_indices)) {
        derived_idx <- param_indices[i]

        # Check that derived is not already fixed or derived
        if (param_fixed[derived_idx]) {
          stop(sprintf("Equality constraint parameter '%s' is already fixed", constraint[i]))
        }
        if (as.character(derived_idx) %in% names(equality_mapping)) {
          stop(sprintf("Parameter '%s' appears in multiple equality constraints", constraint[i]))
        }

        # Mark derived as fixed and store mapping
        param_fixed[derived_idx] <- TRUE
        equality_mapping[[as.character(derived_idx)]] <- primary_idx

        if (verbose) {
          message(sprintf("  Equality constraint: param %d (%s) = param %d (%s)",
                          derived_idx, constraint[i], primary_idx, constraint[1]))
        }
      }
    }
  }

  free_idx <- which(!param_fixed)
  n_free <- length(free_idx)

  return(list(
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    param_fixed = param_fixed,
    free_idx = free_idx,
    n_free = n_free,
    lower_bounds_free = lower_bounds[free_idx],
    upper_bounds_free = upper_bounds[free_idx],
    equality_mapping = equality_mapping  # Maps derived -> primary
  ))
}

# ---- Helper function for eigenvector-based saddle point escape ----

#' Escape saddle point using eigenvector-based restarts
#'
#' When optimization converges to a saddle point (detected by positive Hessian
#' eigenvalues), this function perturbs parameters in the direction of positive
#' eigenvectors to escape the saddle point.
#'
#' @param params Current parameter estimates at saddle point
#' @param hessian_fn Function to compute Hessian at given parameters
#' @param objective_fn Function to compute objective at given parameters
#' @param param_constraints List with lower/upper bounds and free_idx
#' @param verbose Whether to print progress
#' @return List with new_params and found_escape (TRUE if escape direction found)
#' @keywords internal
find_saddle_escape_direction <- function(params, hessian_fn, objective_fn,
                                          param_constraints, verbose = FALSE) {
  # Compute Hessian at current point
  hess_mat <- hessian_fn(params)

  # Get eigenvalue decomposition
  eig <- eigen(hess_mat, symmetric = TRUE)

  # For minimization: saddle point has at least one NEGATIVE eigenvalue

  # (positive eigenvalues are good - they indicate local minimum directions)
  neg_idx <- which(eig$values < -1e-6)

  if (length(neg_idx) == 0) {
    if (verbose) message("  No negative eigenvalues found - not a saddle point")
    return(list(new_params = params, found_escape = FALSE))
  }

  if (verbose) {
    message(sprintf("  Found %d negative eigenvalue(s): %s",
                   length(neg_idx),
                   paste(round(eig$values[neg_idx], 2), collapse = ", ")))
  }

  # Current objective value
  current_obj <- objective_fn(params)
  best_obj <- current_obj
  best_params <- params

  # Try moving in direction of negative eigenvector(s)
  for (idx in neg_idx) {
    eigvec <- eig$vectors[, idx]

    # Scale step: use a minimum step size to avoid getting stuck with tiny steps
    # Larger negative eigenvalues suggest steeper escape, but ensure at least 0.1
    base_step <- max(0.1, min(1.0, sqrt(abs(eig$values[idx]))))

    # Try both directions with different step sizes
    for (direction in c(-1, 1)) {
      for (step_mult in c(0.5, 1.0, 2.0, 5.0)) {
        step <- base_step * step_mult

        new_params <- params + direction * step * eigvec

        # Apply bounds
        new_params <- pmax(new_params, param_constraints$lower_bounds_free)
        new_params <- pmin(new_params, param_constraints$upper_bounds_free)

        # Evaluate objective at new point
        new_obj <- tryCatch(
          objective_fn(new_params),
          error = function(e) Inf
        )

        if (is.finite(new_obj) && new_obj < best_obj) {
          best_obj <- new_obj
          best_params <- new_params
          if (verbose) {
            message(sprintf("    Eigenvector %d, dir=%+d, step=%.2f: obj %.4f -> %.4f",
                           idx, direction, step, current_obj, new_obj))
          }
        }
      }
    }
  }

  found_escape <- best_obj < current_obj - 1e-6
  if (verbose && found_escape) {
    message(sprintf("  Found escape direction: obj improved by %.4f", current_obj - best_obj))
  }

  return(list(new_params = best_params, found_escape = found_escape))
}

# ---- Adaptive quadrature summary ----

#' Print adaptive quadrature summary
#'
#' Computes and displays statistics about adaptive integration point allocation.
#' This is used internally to show the distribution of integration points
#' across observations when adaptive quadrature is enabled.
#'
#' @param factor_ses Matrix (n_obs x n_factors) of factor score standard errors
#' @param factor_vars Numeric vector of factor variances
#' @param threshold Threshold for determining quadrature points
#' @param max_quad Maximum quadrature points per factor
#' @return Invisible list with summary statistics
#' @keywords internal
print_adaptive_quadrature_summary <- function(factor_ses, factor_vars, threshold, max_quad) {
  n_obs <- nrow(factor_ses)
  n_fac <- ncol(factor_ses)

  # Compute per-observation total integration points
  nquad_counts <- integer(0)
  total_points <- 0
  total_points_standard <- max_quad^n_fac

  for (i in seq_len(n_obs)) {
    obs_total <- 1
    for (j in seq_len(n_fac)) {
      f_se <- factor_ses[i, j]
      f_var <- factor_vars[j]
      ratio <- f_se / f_var / threshold
      nq <- 1 + 2 * floor(ratio)
      if (nq < 1) nq <- 1
      if (nq > max_quad) nq <- max_quad
      if (f_se > sqrt(f_var)) nq <- max_quad
      obs_total <- obs_total * nq
    }
    key <- as.character(obs_total)
    if (is.na(nquad_counts[key])) nquad_counts[key] <- 0
    nquad_counts[key] <- nquad_counts[key] + 1
    total_points <- total_points + obs_total
  }

  avg_points <- total_points / n_obs
  reduction <- 100 * (1 - avg_points / total_points_standard)

  # Print summary
  message("\nAdaptive Integration Summary")
  message("----------------------------")
  message(sprintf("Threshold: %.2f, Max quad points: %d\n", threshold, max_quad))

  message("Integration points per observation:")
  message("  Points   Observations   Percent")

  # Sort by number of points
  sorted_keys <- sort(as.numeric(names(nquad_counts)))
  for (pts in sorted_keys) {
    key <- as.character(pts)
    count <- nquad_counts[key]
    pct <- 100 * count / n_obs
    message(sprintf("%8d%15d%10.1f%%", pts, count, pct))
  }

  message(sprintf("\nAverage integration points: %.1f (vs %d standard)",
                 avg_points, as.integer(total_points_standard)))
  message(sprintf("Computational reduction: %.1f%%\n", reduction))

  invisible(list(
    avg_points = avg_points,
    total_points_standard = total_points_standard,
    reduction_pct = reduction,
    distribution = nquad_counts
  ))
}

# ---- Main estimation function ----

#' Estimate factor model using R-based optimization
#'
#' This function estimates a factor model by optimizing the likelihood using
#' C++ for fast evaluation and R for optimization and parallelization.
#'
#' @param model_system A model_system object from define_model_system()
#' @param data Data frame containing all variables
#' @param init_params Initial parameter values (optional)
#' @param control Estimation control object from define_estimation_control()
#' @param optimizer Optimizer to use (default "nloptr"):
#'   \itemize{
#'     \item \code{"nloptr"} - L-BFGS via nloptr (gradient only, no Hessian)
#'     \item \code{"optim"} - L-BFGS-B via stats::optim (gradient only, no Hessian)
#'     \item \code{"nlminb"} - Uses analytical gradient AND Hessian (more efficient!)
#'     \item \code{"trust"} - Trust region method with analytical Hessian (requires trustOptim package)
#'   }
#' @param parallel Whether to use parallel computation (default TRUE)
#' @param verbose Whether to print progress (default TRUE)
#' @param max_restarts Maximum number of eigenvector-based restarts for escaping
#'   saddle points (default 5). Set to 0 to disable.
#' @param factor_scores Matrix (n_obs x n_factors) of factor score estimates from
#'   a previous stage. Used for adaptive quadrature. (default NULL)
#' @param factor_ses Matrix (n_obs x n_factors) of factor score standard errors
#'   from a previous stage. Used for adaptive quadrature. (default NULL)
#' @param factor_vars Named numeric vector of factor variances from a previous
#'   stage. Used for adaptive quadrature. (default NULL)
#' @param init_factor_scores Matrix (n_obs x n_factors) of factor scores to use
#'   for initializing factor loadings. When provided, loadings are estimated by
#'   treating factor scores as regressors, giving better starting values than
#'   the default (0.5). Useful for two-stage estimation where Stage 1 factor
#'   scores can improve Stage 2 initialization. (default NULL)
#' @param checkpoint_file Path to file for saving checkpoint parameters during
#'   optimization. When specified, parameters are saved each time the Hessian is
#'   evaluated at a point with improved likelihood. Useful for long-running
#'   estimations that may need to be restarted. The file contains parameter names
#'   and values as CSV with metadata headers. (default NULL = no checkpointing)
#'
#' @details
#' For maximum efficiency, use \code{optimizer = "nlminb"} or \code{optimizer = "trust"}
#' which exploit the analytical Hessian computed in C++. The default L-BFGS methods
#' only use the gradient and approximate the Hessian from gradient history.
#'
#' When optimization fails to converge (possibly at a saddle point), the function
#' will attempt to escape by moving in the direction of negative Hessian eigenvalues.
#' This is controlled by the \code{max_restarts} parameter.
#'
#' @return List with parameter estimates, standard errors, log-likelihood, etc.
#' @examples
#' \donttest{
#' # Simulate a simple one-factor model with two linear indicators
#' set.seed(1); n <- 100
#' f <- rnorm(n)
#' dat <- data.frame(intercept = 1,
#'   y1 = 1.0 * f + rnorm(n, 0, 0.5),
#'   y2 = 0.8 * f + rnorm(n, 0, 0.5))
#' fm <- define_factor_model(n_factors = 1)
#' mc1 <- define_model_component("m1", dat, "y1", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = 1)
#' mc2 <- define_model_component("m2", dat, "y2", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = NA_real_)
#' ms <- define_model_system(components = list(mc1, mc2), factor = fm)
#' ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
#' result <- estimate_model_rcpp(ms, dat, control = ctrl,
#'   optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
#' result$estimates
#' }
#' @export
estimate_model_rcpp <- function(model_system, data, init_params = NULL,
                                control = NULL, optimizer = "nlminb",
                                parallel = TRUE, verbose = TRUE,
                                max_restarts = 5,
                                factor_scores = NULL, factor_ses = NULL,
                                factor_vars = NULL,
                                init_factor_scores = NULL,
                                checkpoint_file = NULL) {

  # WORKAROUND: Deep copy model_system to avoid C++ reuse bug
  # Use serialize/unserialize for true deep copy
  model_system <- unserialize(serialize(model_system, NULL))

  # Validate init_params (catch common mistake of passing control as 3rd positional arg)
  if (!is.null(init_params) && inherits(init_params, "estimation_control")) {
    stop("'init_params' appears to be a control object. Did you mean: control = ...?
  Correct usage: estimate_model_rcpp(model_system, data, control = my_control)")
  }

  # Default control if not provided
  if (is.null(control)) {
    control <- define_estimation_control()
  }

  # Print version info
  if (verbose) {
    pkg_version <- utils::packageVersion("factorana")
    message(sprintf("factorana version %s", pkg_version))
  }

  # Convert data to matrix
  data_mat <- as.matrix(data)

  # Validate observation weights if specified
  if (!is.null(model_system$weights)) {
    weights_var <- model_system$weights
    if (!weights_var %in% colnames(data_mat)) {
      stop(sprintf("Weights variable '%s' not found in data", weights_var))
    }
    weights_vec <- data_mat[, weights_var]
    if (any(is.na(weights_vec))) {
      stop("Observation weights contain NA values")
    }
    if (any(weights_vec <= 0)) {
      warning("Some observation weights are <= 0. This may cause issues.")
    }
  }

  # Setup parallel cluster if requested
  cl <- NULL
  if (parallel && control$num_cores > 1) {
    # Split data across workers first to determine actual number of workers needed
    n_obs <- nrow(data_mat)
    n_per_worker <- ceiling(n_obs / control$num_cores)
    data_splits <- split(1:n_obs, ceiling(1:n_obs / n_per_worker))

    # Use actual number of data splits (may be less than requested cores)
    n_workers <- length(data_splits)

    if (verbose) {
      message(sprintf("Setting up parallel cluster with %d cores (requested %d)...",
                     n_workers, control$num_cores))
    }

    # Determine cluster type based on control setting
    cluster_type <- control$cluster_type
    if (is.null(cluster_type)) cluster_type <- "auto"  # Backwards compatibility

    if (cluster_type == "auto") {
      # Auto: use FORK on Unix, PSOCK on Windows
      if (.Platform$OS.type == "unix") {
        cluster_type <- "FORK"
      } else {
        cluster_type <- "PSOCK"
      }
    }

    # Create the cluster
    if (cluster_type == "FORK") {
      cl <- tryCatch({
        if (verbose) message("Using FORK cluster (shared memory)...")
        parallel::makeForkCluster(n_workers)
      }, error = function(e) {
        if (verbose) message("FORK cluster failed, using PSOCK...")
        parallel::makeCluster(n_workers)
      })
    } else {
      if (verbose) message("Using PSOCK cluster...")
      cl <- parallel::makeCluster(n_workers)
    }
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel estimation. ",
           "Install with: install.packages('doParallel')")
    }
    doParallel::registerDoParallel(cl)

    # Robust cleanup: use on.exit AND store cluster in parent frame for interrupt handling
    # The on.exit handles normal exits and errors
    on.exit({
      if (!is.null(cl)) {
        tryCatch(
          parallel::stopCluster(cl),
          error = function(e) NULL  # Ignore errors during cleanup
        )
      }
    }, add = TRUE)

  } else {
    # Single worker
    data_splits <- list(1:nrow(data_mat))
  }

  # Get initial parameters FIRST (needed for fixed coefficient values in C++)
  factor_variance_fixed <- NULL
  full_init_params <- NULL  # Full parameter vector (including fixed)

  if (is.null(init_params)) {
    # Use smart initialization based on separate component estimation
    # If init_factor_scores provided, use them to estimate loadings
    init_result <- initialize_parameters(model_system, data,
                                          factor_scores = init_factor_scores,
                                          verbose = verbose)
    full_init_params <- init_result$init_params
    factor_variance_fixed <- init_result$factor_variance_fixed
  } else {
    # User provided init_params (assumed to be full vector)
    full_init_params <- init_params
  }

  # Initialize factor models on each worker
  if (verbose) message("Initializing C++ likelihood evaluators...")

  # Get n_quad from control
  n_quad <- control$n_quad_points

  if (!is.null(cl)) {
    # Check if this is a FORK cluster
    is_fork <- inherits(cl[[1]], "forknode")

    if (is_fork) {
      # FORK cluster: load library on workers (they inherit parent's state but need explicit library)
      parallel::clusterEvalQ(cl, library(factorana))
    } else {
      # PSOCK cluster: need to set library paths and load factorana
      current_lib_paths <- .libPaths()
      parallel::clusterExport(cl, "current_lib_paths", envir = environment())
      parallel::clusterEvalQ(cl, {
        .libPaths(current_lib_paths)
        library(factorana)
      })
    }

    # Export necessary objects to workers
    parallel::clusterExport(cl, c("model_system", "data_mat", "n_quad", "data_splits", "full_init_params"),
                           envir = environment())

    # Set worker IDs (1 to n_workers) in each worker's global environment
    for (i in seq_along(cl)) {
      parallel::clusterCall(cl[i], assign, ".self_id", i, envir = .GlobalEnv)
    }

    # Initialize FactorModel on each worker with its data subset
    # IMPORTANT: Store pointer in worker's global env to avoid serialization
    parallel::clusterEvalQ(cl, {
      worker_id <- .self_id
      idx <- data_splits[[worker_id]]
      data_subset <- data_mat[idx, , drop = FALSE]
      .fm_ptr <- initialize_factor_model_cpp(model_system, data_subset, n_quad, full_init_params)
    })

    fm_ptrs <- NULL  # Not used; pointers stored on workers

    # Set observation weights if specified in model_system
    if (!is.null(model_system$weights)) {
      weights_var <- model_system$weights
      # Export weights variable name
      parallel::clusterExport(cl, "weights_var", envir = environment())
      # Each worker extracts and sets weights for its data subset
      parallel::clusterEvalQ(cl, {
        worker_id <- .self_id
        idx <- data_splits[[worker_id]]
        weights_subset <- data_mat[idx, weights_var]
        set_observation_weights_cpp(.fm_ptr, weights_subset)
      })
      if (verbose) {
        weights_vec <- data_mat[, weights_var]
        message(sprintf("Using observation weights from '%s' (range: %.3f to %.3f)",
                       weights_var, min(weights_vec), max(weights_vec)))
      }
    }

    # Set up adaptive quadrature if factor_scores provided
    if (!is.null(factor_scores) && !is.null(factor_ses) && !is.null(factor_vars)) {
      adapt_thresh <- control$adapt_int_thresh
      adapt_max_quad <- n_quad
      # Export adaptive quadrature settings to workers
      parallel::clusterExport(cl, c("factor_scores", "factor_ses", "factor_vars",
                                    "adapt_thresh", "adapt_max_quad"),
                             envir = environment())
      # Each worker sets up adaptive quadrature for its data subset
      parallel::clusterEvalQ(cl, {
        worker_id <- .self_id
        idx <- data_splits[[worker_id]]
        scores_subset <- factor_scores[idx, , drop = FALSE]
        ses_subset <- factor_ses[idx, , drop = FALSE]
        set_adaptive_quadrature_cpp(.fm_ptr, scores_subset, ses_subset, factor_vars,
                                    adapt_thresh, adapt_max_quad, verbose = FALSE)
      })
      # Print summary from R (works with parallel output)
      if (verbose) {
        print_adaptive_quadrature_summary(factor_ses, factor_vars, adapt_thresh, adapt_max_quad)
      }
    }
  } else {
    # Single worker initialization
    fm_ptrs <- list(initialize_factor_model_cpp(model_system, data_mat, n_quad, full_init_params))

    # Set observation weights if specified in model_system
    if (!is.null(model_system$weights)) {
      weights_var <- model_system$weights
      weights_vec <- data_mat[, weights_var]
      set_observation_weights_cpp(fm_ptrs[[1]], weights_vec)
      if (verbose) {
        message(sprintf("Using observation weights from '%s' (range: %.3f to %.3f)",
                       weights_var, min(weights_vec), max(weights_vec)))
      }
    }

    # Set up adaptive quadrature if factor_scores provided
    if (!is.null(factor_scores) && !is.null(factor_ses) && !is.null(factor_vars)) {
      adapt_thresh <- control$adapt_int_thresh
      adapt_max_quad <- n_quad
      set_adaptive_quadrature_cpp(fm_ptrs[[1]], factor_scores, factor_ses, factor_vars,
                                  adapt_thresh, adapt_max_quad, verbose = FALSE)
      # Print summary from R for consistency with parallel mode
      if (verbose) {
        print_adaptive_quadrature_summary(factor_ses, factor_vars, adapt_thresh, adapt_max_quad)
      }
    }
  }

  # Get parameter count and extract free parameters
  if (!is.null(cl)) {
    # For parallel case, get from first worker
    param_info <- parallel::clusterEvalQ(cl, get_parameter_info_cpp(.fm_ptr))[[1]]
  } else {
    param_info <- get_parameter_info_cpp(fm_ptrs[[1]])
  }
  n_params_free <- param_info$n_param_free
  n_params_total <- param_info$n_param

  # Extract free parameters from full init_params for the optimizer
  if (n_params_free == n_params_total) {
    # No fixed params - use full vector
    init_params <- full_init_params
  } else {
    # Some params are fixed - extract only free ones
    if (!is.null(cl)) {
      # For parallel case, use first worker to extract
      parallel::clusterExport(cl, "full_init_params", envir = environment())
      init_params <- parallel::clusterEvalQ(cl, {
        extract_free_params_cpp(.fm_ptr, full_init_params)
      })[[1]]
    } else {
      init_params <- extract_free_params_cpp(fm_ptrs[[1]], full_init_params)
    }

    if (verbose) {
      n_fixed <- n_params_total - n_params_free
      message(sprintf("  %d parameters total, %d fixed, %d free", n_params_total, n_fixed, n_params_free))
    }
  }

  # Validate parameter count
  if (length(init_params) != n_params_free) {
    stop(sprintf("Parameter initialization returned %d free parameters but C++ expects %d",
                 length(init_params), n_params_free))
  }

  # Build parameter metadata (names, types, etc.)
  param_metadata <- build_parameter_metadata(model_system)
  param_names <- param_metadata$names  # For checkpointing

  # Validate parameter counts match between R metadata and C++
  if (param_metadata$n_params != n_params_total) {
    warning(sprintf(
      paste0("Parameter count mismatch: R metadata has %d params but C++ has %d. ",
             "This may cause parameter names to be misaligned. ",
             "Check build_parameter_metadata() for multinomial logit components."),
      param_metadata$n_params, n_params_total
    ))
    # Also warn about which component might be causing the issue
    if (verbose) {
      message("Debugging parameter count mismatch:")
      message(sprintf("  param_metadata$n_params = %d", param_metadata$n_params))
      message(sprintf("  n_params_total (from C++) = %d", n_params_total))
      message(sprintf("  length(full_init_params) = %d", length(full_init_params)))
    }
  }

  # Setup constraints and identify fixed/free parameters
  # Use full_init_params (not init_params) because param_metadata covers ALL parameters
  param_constraints <- setup_parameter_constraints(model_system, full_init_params, param_metadata, factor_variance_fixed, verbose)

  if (verbose) {
    message(sprintf("Total parameters: %d", n_params_total))
    message(sprintf("Free parameters: %d", n_params_free))
    if (any(param_constraints$param_fixed)) {
      fixed_names <- param_metadata$names[param_constraints$param_fixed]
      message(sprintf("Fixed parameters (%d): %s",
                     sum(param_constraints$param_fixed),
                     paste(fixed_names, collapse = ", ")))
    }
    n_sigma <- sum(param_metadata$types == "sigma")
    if (n_sigma > 0) {
      message(sprintf("Sigma parameters (%d) have lower bound = 0.01", n_sigma))
    }
  }

  # Benchmark likelihood/gradient/Hessian computation before optimization
  # NOTE: C++ functions expect FREE params only (not full params)
  if (verbose) {
    message("\nBenchmarking computation times (single evaluation)...")

    # Benchmark log-likelihood only
    t_loglik <- system.time({
      if (!is.null(cl)) {
        parallel::clusterExport(cl, "init_params", envir = environment())
        loglik_parts <- parallel::clusterEvalQ(cl, {
          evaluate_loglik_only_cpp(.fm_ptr, init_params)
        })
        loglik_test <- sum(unlist(loglik_parts))
      } else {
        loglik_test <- evaluate_loglik_only_cpp(fm_ptrs[[1]], init_params)
      }
    })[3]

    # Benchmark gradient
    t_grad <- system.time({
      if (!is.null(cl)) {
        grad_parts <- parallel::clusterEvalQ(cl, {
          result <- evaluate_likelihood_cpp(.fm_ptr, init_params,
                                           compute_gradient = TRUE,
                                           compute_hessian = FALSE)
          result$gradient
        })
      } else {
        result <- evaluate_likelihood_cpp(fm_ptrs[[1]], init_params,
                                         compute_gradient = TRUE,
                                         compute_hessian = FALSE)
      }
    })[3]

    # Benchmark Hessian
    t_hess <- system.time({
      if (!is.null(cl)) {
        hess_parts <- parallel::clusterEvalQ(cl, {
          result <- evaluate_likelihood_cpp(.fm_ptr, init_params,
                                           compute_gradient = FALSE,
                                           compute_hessian = TRUE)
          result$hessian
        })
      } else {
        result <- evaluate_likelihood_cpp(fm_ptrs[[1]], init_params,
                                         compute_gradient = FALSE,
                                         compute_hessian = TRUE)
      }
    })[3]

    message(sprintf("  Log-likelihood:  %.3f sec", t_loglik))
    message(sprintf("  Gradient:        %.3f sec", t_grad))
    message(sprintf("  Hessian:         %.3f sec", t_hess))
    message(sprintf("  Initial loglik:  %.4f", loglik_test))
    message("")
  }

  # Helper function to apply equality constraints to full parameter vector
  apply_equality_constraints <- function(params_full) {
    if (length(param_constraints$equality_mapping) > 0) {
      for (derived_idx_str in names(param_constraints$equality_mapping)) {
        derived_idx <- as.integer(derived_idx_str)
        primary_idx <- param_constraints$equality_mapping[[derived_idx_str]]
        params_full[derived_idx] <- params_full[primary_idx]
      }
    }
    return(params_full)
  }

  # Checkpointing: track best likelihood for smart checkpoint saving
  # Use an environment so closures can modify these values
  checkpoint_env <- new.env(parent = emptyenv())
  checkpoint_env$best_loglik <- -Inf
  checkpoint_env$best_params_free <- NULL
  checkpoint_env$n_checkpoints <- 0

  # Define objective function (operates on free parameters only)
  objective_fn <- function(params_free) {
    # C++ expects only FREE params (it has fixed values stored internally)
    # params_free is already the free params from the optimizer
    # We just need to pass it directly to C++

    if (!is.null(cl)) {
      # Parallel evaluation: aggregate across workers
      parallel::clusterExport(cl, "params_free", envir = environment())
      loglik_parts <- parallel::clusterEvalQ(cl, {
        evaluate_loglik_only_cpp(.fm_ptr, params_free)
      })
      loglik <- sum(unlist(loglik_parts))
    } else {
      # Single worker
      loglik <- evaluate_loglik_only_cpp(fm_ptrs[[1]], params_free)
    }

    # Track best likelihood for checkpointing
    if (loglik > checkpoint_env$best_loglik) {
      checkpoint_env$best_loglik <- loglik
      checkpoint_env$best_params_free <- params_free
    }

    return(-loglik)  # Negative for minimization
  }

  # Define gradient function (returns gradient for free parameters only)
  gradient_fn <- function(params_free) {
    # C++ expects only FREE params (it has fixed values stored internally)
    # C++ returns gradient for free params only

    if (!is.null(cl)) {
      # Parallel evaluation: aggregate gradients
      parallel::clusterExport(cl, "params_free", envir = environment())
      grad_parts <- parallel::clusterEvalQ(cl, {
        result <- evaluate_likelihood_cpp(.fm_ptr, params_free,
                                         compute_gradient = TRUE,
                                         compute_hessian = FALSE)
        result$gradient
      })
      grad_free <- Reduce(`+`, grad_parts)
    } else {
      result <- evaluate_likelihood_cpp(fm_ptrs[[1]], params_free,
                                       compute_gradient = TRUE,
                                       compute_hessian = FALSE)
      grad_free <- result$gradient
    }

    # C++ returns gradient for free parameters only, already extracted
    # Return negative for minimization
    return(-grad_free)
  }

  # Define Hessian function (returns Hessian for free parameters only)
  hessian_fn <- function(params_free) {
    # C++ expects only FREE params (it has fixed values stored internally)
    # C++ returns Hessian for free params only

    if (!is.null(cl)) {
      # Parallel evaluation: aggregate Hessians
      parallel::clusterExport(cl, "params_free", envir = environment())
      hess_parts <- parallel::clusterEvalQ(cl, {
        result <- evaluate_likelihood_cpp(.fm_ptr, params_free,
                                         compute_gradient = FALSE,
                                         compute_hessian = TRUE)
        result$hessian
      })
      # Aggregate Hessian matrices (they're upper triangles, addition is correct)
      hess <- Reduce(`+`, hess_parts)
    } else {
      result <- evaluate_likelihood_cpp(fm_ptrs[[1]], params_free,
                                       compute_gradient = FALSE,
                                       compute_hessian = TRUE)
      hess <- result$hessian
    }

    # C++ returns Hessian as upper triangle vector for FREE params
    # Convert to full symmetric matrix
    n_params_free <- length(params_free)
    hess_vec_len <- length(hess)
    expected_len <- n_params_free * (n_params_free + 1) / 2

    if (verbose && hess_vec_len != expected_len) {
      message(sprintf("WARNING: Hessian vector length mismatch!"))
      message(sprintf("  Expected: %d (for %d free params)", expected_len, n_params_free))
      message(sprintf("  Actual: %d", hess_vec_len))
    }

    hess_mat_free <- matrix(0, n_params_free, n_params_free)
    idx <- 1
    for (i in 1:n_params_free) {
      for (j in i:n_params_free) {
        if (idx > hess_vec_len) {
          if (verbose) {
            message(sprintf("ERROR: Trying to access hess[%d] but length is only %d", idx, hess_vec_len))
          }
          stop("Hessian vector length insufficient for reconstruction")
        }
        hess_mat_free[i, j] <- hess[idx]
        hess_mat_free[j, i] <- hess[idx]  # Symmetrize
        idx <- idx + 1
      }
    }

    # Debug: Check for NA/NaN values
    if (any(is.na(hess_mat_free)) || any(is.infinite(hess_mat_free))) {
      if (verbose) {
        message("WARNING: Hessian contains NA/NaN/Inf values!")
        message(sprintf("  NA count: %d", sum(is.na(hess_mat_free))))
        message(sprintf("  Inf count: %d", sum(is.infinite(hess_mat_free))))
      }
    }

    # Checkpointing: save parameters when Hessian is evaluated at best point
    if (!is.null(checkpoint_file) && !is.null(checkpoint_env$best_params_free)) {
      # Check if current params are at or near the best
      params_match <- isTRUE(all.equal(params_free, checkpoint_env$best_params_free,
                                        tolerance = 1e-10))
      if (params_match) {
        # Reconstruct full parameter vector from free parameters
        params_full_checkpoint <- full_init_params  # Start with initial full params
        params_full_checkpoint[param_constraints$free_idx] <- params_free
        params_full_checkpoint <- apply_equality_constraints(params_full_checkpoint)

        # Save checkpoint with parameter names
        checkpoint_data <- data.frame(
          parameter = param_names,
          value = params_full_checkpoint,
          stringsAsFactors = FALSE
        )

        # Add metadata
        checkpoint_header <- sprintf(
          "# Checkpoint saved: %s\n# Log-likelihood: %.8f\n# Iteration: %d\n",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          checkpoint_env$best_loglik,
          checkpoint_env$n_checkpoints + 1
        )

        # Write to file
        tryCatch({
          con <- file(checkpoint_file, "w")
          writeLines(checkpoint_header, con)
          utils::write.csv(checkpoint_data, con, row.names = FALSE)
          close(con)
          checkpoint_env$n_checkpoints <- checkpoint_env$n_checkpoints + 1
          if (verbose) {
            message(sprintf("  [Checkpoint %d saved: LL=%.4f]",
                           checkpoint_env$n_checkpoints, checkpoint_env$best_loglik))
          }
        }, error = function(e) {
          if (verbose) {
            message(sprintf("  [Warning: Failed to save checkpoint: %s]", e$message))
          }
        })
      }
    }

    return(-hess_mat_free)  # Negative for minimization
  }

  # Run optimization with eigenvector-based restarts
  if (verbose) message(sprintf("Running optimization using %s...", optimizer))

  # Track current starting point for restarts
  # IMPORTANT: Use init_params (extracted via C++) not param_constraints$free_idx
  # because param_metadata may have different parameter count than C++
  current_start <- init_params
  n_restarts_used <- 0

  # Track optimization timing
  optimization_start_time <- Sys.time()

  # Helper function to run one optimization attempt
  run_one_optimization <- function(start_params, verbose_opt) {
    if (optimizer == "nloptr") {
      if (!requireNamespace("nloptr", quietly = TRUE)) {
        stop("nloptr package required but not installed")
      }

      opt_result <- nloptr::nloptr(
        x0 = start_params,
        eval_f = objective_fn,
        eval_grad_f = gradient_fn,
        lb = param_constraints$lower_bounds_free,
        ub = param_constraints$upper_bounds_free,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          maxeval = 1000,
          xtol_rel = 1e-6,
          print_level = if (verbose_opt) 2 else 0
        )
      )
      return(list(
        par = opt_result$solution,
        value = opt_result$objective,
        convergence = opt_result$status
      ))

    } else if (optimizer == "optim") {
      opt_result <- stats::optim(
        par = start_params,
        fn = objective_fn,
        gr = gradient_fn,
        method = "L-BFGS-B",
        lower = param_constraints$lower_bounds_free,
        upper = param_constraints$upper_bounds_free,
        control = list(trace = if (verbose_opt) 1 else 0)
      )
      return(list(
        par = opt_result$par,
        value = opt_result$value,
        convergence = opt_result$convergence
      ))

    } else if (optimizer == "nlminb") {
      if (verbose_opt) message("  Using analytical gradient and Hessian")

      opt_result <- stats::nlminb(
        start = start_params,
        objective = objective_fn,
        gradient = gradient_fn,
        hessian = hessian_fn,
        lower = param_constraints$lower_bounds_free,
        upper = param_constraints$upper_bounds_free,
        control = list(
          trace = if (verbose_opt) 1 else 0,
          eval.max = 5000,
          iter.max = 10000
        )
      )
      return(list(
        par = opt_result$par,
        value = opt_result$objective,
        convergence = opt_result$convergence,
        iterations = opt_result$iterations,
        evaluations = opt_result$evaluations
      ))

    } else if (optimizer == "trust") {
      if (!requireNamespace("trustOptim", quietly = TRUE)) {
        stop("trustOptim package required but not installed. Install with: install.packages('trustOptim')")
      }

      if (any(param_constraints$param_fixed)) {
        warning("trust optimizer does not support fixed parameters. ",
                "Fixed parameters may drift from their initial values. ",
                "Consider using 'nlminb' or 'optim' instead.")
      }

      if (verbose_opt) message("  Using trust region with analytical gradient and Hessian")

      opt_result <- trustOptim::trust.optim(
        x = start_params,
        fn = objective_fn,
        gr = gradient_fn,
        hs = hessian_fn,
        method = "SR1",
        control = list(
          report.level = if (verbose_opt) 2 else 0,
          maxit = 500
        )
      )
      return(list(
        par = opt_result$solution,
        value = opt_result$value,
        convergence = if (opt_result$converged) 0 else 1
      ))

    } else {
      stop("Unknown optimizer: ", optimizer, "\n",
           "Available: 'nloptr' (L-BFGS), 'optim' (L-BFGS-B), 'nlminb' (uses Hessian), 'trust' (uses Hessian)")
    }
  }

  # Determine convergence success based on optimizer
  is_converged <- function(conv_code) {
    if (optimizer %in% c("optim", "nlminb")) {
      return(conv_code == 0)
    } else if (optimizer == "nloptr") {
      return(conv_code >= 1 && conv_code <= 4)  # nloptr success codes
    } else if (optimizer == "trust") {
      return(conv_code == 0)  # Already converted to 0 = success
    }
    return(FALSE)
  }

  # Run initial optimization
  opt_result <- run_one_optimization(current_start, verbose)
  estimates_free <- opt_result$par
  loglik <- -opt_result$value
  convergence <- opt_result$convergence

  # Eigenvector-based restart loop
  if (max_restarts > 0 && !is_converged(convergence)) {
    for (restart in seq_len(max_restarts)) {
      if (verbose) {
        message(sprintf("  Optimization did not converge. Attempting eigenvector-based restart %d/%d...",
                       restart, max_restarts))
      }

      # Try to find escape direction using eigenvector analysis
      escape_result <- find_saddle_escape_direction(
        params = estimates_free,
        hessian_fn = hessian_fn,
        objective_fn = objective_fn,
        param_constraints = param_constraints,
        verbose = verbose
      )

      if (!escape_result$found_escape) {
        if (verbose) message("  No escape direction found. Stopping restarts.")
        break
      }

      n_restarts_used <- restart

      # Restart from new point
      current_start <- escape_result$new_params
      opt_result <- run_one_optimization(current_start, verbose)
      estimates_free <- opt_result$par
      loglik <- -opt_result$value
      convergence <- opt_result$convergence

      if (verbose) {
        message(sprintf("  Restart %d: loglik = %.4f, conv = %d",
                       restart, loglik, convergence))
      }

      if (is_converged(convergence)) {
        if (verbose) message(sprintf("  Converged after %d restart(s)!", restart))
        break
      }
    }
  }

  # Reconstruct full parameter vector
  # Use full_init_params (not init_params) because free_idx indexes into the FULL parameter vector
  estimates <- full_init_params
  estimates[param_constraints$free_idx] <- estimates_free

  # Apply equality constraints to final estimates
  estimates <- apply_equality_constraints(estimates)

  # Compute Hessian for standard errors
  # C++ expects FREE params only
  if (verbose) message("Computing Hessian for standard errors...")

  if (!is.null(cl)) {
    parallel::clusterExport(cl, "estimates_free", envir = environment())
    hess_parts <- parallel::clusterEvalQ(cl, {
      result <- evaluate_likelihood_cpp(.fm_ptr, estimates_free,
                                       compute_gradient = TRUE,
                                       compute_hessian = TRUE)
      result$hessian
    })
    # Aggregate Hessian matrices (they're upper triangles, addition is correct)
    hessian <- Reduce(`+`, hess_parts)
  } else {
    result <- evaluate_likelihood_cpp(fm_ptrs[[1]], estimates_free,
                                     compute_gradient = TRUE,
                                     compute_hessian = TRUE)
    hessian <- result$hessian
  }

  # IMPORTANT: The C++ code returns the Hessian of the log-likelihood.
  # Since we're minimizing the NEGATIVE log-likelihood, we need to negate
  # the Hessian to get the correct second derivatives.
  hessian <- -hessian

  # Compute standard errors from Hessian
  # The Hessian is for FREE params only; fixed params get SE = 0
  # Initialize full vector for all params
  std_errors <- rep(NA, length(estimates))
  names(std_errors) <- names(estimates)

  tryCatch({
    # C++ returns Hessian for FREE params only as upper triangle vector
    # Convert to full symmetric matrix
    n_free <- length(estimates_free)
    hess_free <- matrix(0, n_free, n_free)

    # Fill in the upper triangle
    idx <- 1
    for (i in 1:n_free) {
      for (j in i:n_free) {
        hess_free[i, j] <- hessian[idx]
        if (i != j) {
          hess_free[j, i] <- hessian[idx]  # Symmetric
        }
        idx <- idx + 1
      }
    }

    # Verify the matrix is symmetric
    if (verbose) {
      max_asym <- max(abs(hess_free - t(hess_free)))
      if (max_asym > 1e-10) {
        warning(sprintf("Hessian matrix is not symmetric (max diff: %.2e)", max_asym))
      }
    }

    # Identify fixed vs free parameters for result indexing
    fixed_params <- which(param_constraints$param_fixed)
    free_params <- param_constraints$free_idx

    if (n_free > 0) {
      # Check condition number of Hessian
      if (verbose) {
        eig_vals <- eigen(hess_free, only.values = TRUE)$values
        cond_num <- max(abs(eig_vals)) / min(abs(eig_vals))
        message(sprintf("  Hessian condition number: %.2e", cond_num))
        message(sprintf("  Eigenvalues: min=%.2e, max=%.2e", min(abs(eig_vals)), max(abs(eig_vals))))
        if (cond_num > 1e10) {
          message("  Warning: Hessian is poorly conditioned")
        }
      }

      # Invert the Hessian to get the covariance matrix for free parameters
      # Use SVD-based pseudoinverse for numerical stability
      if (verbose) {
        message("\n  Using SVD-based pseudoinverse for numerical stability")
      }

      svd_result <- svd(hess_free)

      # Set tolerance for singular values
      tol <- max(dim(hess_free)) * max(svd_result$d) * .Machine$double.eps * 100  # More conservative
      pos_idx <- svd_result$d > tol

      if (verbose) {
        message(sprintf("  SVD: %d/%d singular values kept (tol=%.2e)", sum(pos_idx), length(pos_idx), tol))
        message(sprintf("  Singular values range: [%.2e, %.2e]", min(svd_result$d), max(svd_result$d)))
      }

      # Compute pseudoinverse
      d_inv <- rep(0, length(svd_result$d))
      d_inv[pos_idx] <- 1 / svd_result$d[pos_idx]

      cov_free <- svd_result$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(svd_result$u)

      # Check if the covariance matrix is positive definite
      if (verbose) {
        cov_eig <- eigen(cov_free, symmetric = TRUE, only.values = TRUE)$values
        n_neg_eig <- sum(cov_eig < -1e-10)  # Allow small numerical errors
        if (n_neg_eig > 0) {
          message(sprintf("  Warning: Covariance matrix has %d negative eigenvalues", n_neg_eig))
          message(sprintf("  Eigenvalues range: [%.2e, %.2e]", min(cov_eig), max(cov_eig)))
        } else {
          message(sprintf("  Covariance matrix is positive semi-definite"))
        }
      }

      # Standard errors for free parameters
      cov_diag <- diag(cov_free)
      se_free <- sqrt(pmax(0, cov_diag))  # pmax ensures non-negative
      std_errors[free_params] <- se_free

      # Fixed parameters have zero standard error
      # Exception: previous_stage parameters should retain their SEs
      std_errors[fixed_params] <- 0.0

      # If we have previous_stage, use those standard errors (match by name)
      if (!is.null(model_system$previous_stage_info)) {
        fixed_param_names <- model_system$previous_stage_info$fixed_param_names
        fixed_std_errors <- model_system$previous_stage_info$fixed_std_errors

        if (length(fixed_param_names) > 0 && length(fixed_std_errors) == length(fixed_param_names)) {
          # Match fixed parameter names to indices in current model
          fixed_indices <- match(fixed_param_names, param_names)
          valid_matches <- !is.na(fixed_indices)

          if (any(valid_matches)) {
            std_errors[fixed_indices[valid_matches]] <- fixed_std_errors[valid_matches]
          }
        }
      }

      # Check for problematic standard errors
      if (verbose) {
        n_small_se <- sum(se_free < 1e-6 & se_free > 0)
        n_zero_se <- sum(se_free == 0)
        if (n_small_se > 0) {
          message(sprintf("\n  Warning: %d parameters have very small SE (< 1e-6)", n_small_se))
        }
        if (n_zero_se > 0) {
          message(sprintf("  Warning: %d parameters have zero SE", n_zero_se))
        }
      }
    }

    if (verbose) {
      n_free <- length(free_params)
      n_fixed <- length(fixed_params)
      message(sprintf("Standard errors computed (%d free, %d fixed parameters)", n_free, n_fixed))
    }
  }, error = function(e) {
    if (verbose) {
      warning("Could not compute standard errors: ", e$message)
    }
  })

  # Calculate total optimization time
  optimization_end_time <- Sys.time()
  optimization_time_secs <- as.numeric(difftime(optimization_end_time, optimization_start_time, units = "secs"))

  if (verbose) {
    # Report convergence status based on optimizer
    if ((optimizer %in% c("optim", "nlminb") && convergence == 0) ||
        (optimizer %in% c("L-BFGS", "L-BFGS-B") && convergence == 0) ||
        (optimizer == "trust" && convergence == TRUE)) {
      message(sprintf("Converged successfully in %.1f seconds. Log-likelihood: %.4f",
                     optimization_time_secs, loglik))
    } else {
      # Convergence failed - determine reason
      if (optimizer == "nlminb") {
        reason <- switch(as.character(convergence),
                        "1" = "iteration limit reached",
                        "non-zero convergence code")
      } else if (optimizer %in% c("optim", "L-BFGS", "L-BFGS-B")) {
        reason <- switch(as.character(convergence),
                        "1" = "iteration limit reached",
                        "10" = "degeneracy in Nelder-Mead simplex",
                        "non-zero convergence code")
      } else {
        reason <- "failed to converge"
      }
      message(sprintf("Convergence FAILED (%s) after %.1f seconds. Log-likelihood: %.4f",
                     reason, optimization_time_secs, loglik))
    }
  }

  # Build equality constraint info for param_table
  equality_tied_to <- rep(NA_character_, length(estimates))
  if (length(param_constraints$equality_mapping) > 0) {
    for (derived_idx_str in names(param_constraints$equality_mapping)) {
      derived_idx <- as.integer(derived_idx_str)
      primary_idx <- param_constraints$equality_mapping[[derived_idx_str]]
      equality_tied_to[derived_idx] <- param_metadata$names[primary_idx]
    }
  }

  # Build param_table for easier access to organized parameter info
  param_table <- data.frame(
    name = param_metadata$names,
    type = param_metadata$types,
    component = sapply(param_metadata$component_id, function(id) {
      if (id == 0) "factor" else model_system$components[[id]]$name
    }),
    estimate = estimates,
    std_error = std_errors,
    fixed = param_constraints$param_fixed,
    tied_to = equality_tied_to,
    stringsAsFactors = FALSE
  )

  # Return results with class for print/summary methods
  result <- list(
    estimates = estimates,
    std_errors = std_errors,
    param_names = param_metadata$names,
    param_table = param_table,
    loglik = loglik,
    convergence = convergence,
    n_restarts = n_restarts_used,
    iterations = opt_result$iterations,
    evaluations = opt_result$evaluations,
    optimization_time = optimization_time_secs,
    model_system = model_system,
    optimizer = optimizer,
    equality_constraints = model_system$equality_constraints
  )
  class(result) <- "factorana_result"
  result
}


#' Helper function to convert model_system to format expected by C++
#'
#' @param model_system model_system object
#' @param data Data frame
#' @return List in format expected by initialize_factor_model_cpp
#' @keywords internal
prepare_model_system_for_cpp <- function(model_system, data) {
  # Extract variable names and create index mapping
  var_names <- names(data)
  var_index <- setNames(seq_along(var_names) - 1, var_names)  # 0-indexed

  # Prepare factor model info
  fm <- model_system$factor
  factor_info <- list(
    n_factors = fm$n_factors,
    n_types = fm$n_types,
    n_mixtures = fm$n_mixtures,
    correlation = fm$correlation,
    loading_normalization = fm$loading_normalization
  )

  # Prepare component info
  components_info <- lapply(model_system$components, function(comp) {
    # Get outcome index
    outcome_idx <- var_index[comp$outcome]

    # Get regressor indices
    regressor_idx <- var_index[comp$covariates]

    # Get missing indicator if present
    missing_idx <- if (!is.null(comp$evaluation_indicator)) {
      var_index[comp$evaluation_indicator]
    } else {
      -1L
    }

    list(
      name = comp$name,
      model_type = comp$model_type,
      outcome_idx = outcome_idx,
      missing_idx = missing_idx,
      regressor_idx = as.integer(regressor_idx),
      factor_normalizations = comp$loading_normalization,
      num_choices = if (comp$model_type == "oprobit") {
        length(unique(data[[comp$outcome]]))
      } else {
        2L
      }
    )
  })

  list(
    factor = factor_info,
    components = components_info
  )
}
