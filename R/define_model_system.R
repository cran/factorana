#' Define a model system
#'
#' @param components A named list of model_component objects
#' @param factor A factor_model object
#' @param previous_stage Optional result from a previous estimation stage.
#'   If provided, the previous stage components and parameters will be fixed
#'   and prepended to the new components. This enables sequential/multi-stage
#'   estimation where early stages are held fixed while later stages are optimized.
#' @param weights Optional. Name of a variable in the data containing observation weights.
#'   When specified, each observation's contribution to the log-likelihood is multiplied
#'   by its weight. Useful for survey weights, importance sampling, or giving different
#'   observations different influence on the estimation. Weights should be positive.
#'   The variable is extracted from the data passed to \code{estimate_model_rcpp()}.
#' @param equality_constraints Optional. A list of character vectors, where each vector
#'   specifies parameter names that should be constrained to be equal during estimation.
#'   The first parameter in each group is the "primary" (freely estimated), and all
#'   other parameters in the group are set equal to the primary.
#'   Example: \code{list(c("Y1_loading_1", "Y2_loading_2"), c("Y1_sigma", "Y2_sigma"))}
#'   This is useful for measurement invariance constraints in longitudinal models.
#' @param free_params Optional. A character vector of parameter names from previous_stage
#'   that should remain FREE (not fixed) in the current stage. This allows selectively
#'   freeing specific parameters while keeping others fixed. Commonly used to free
#'   factor variances while keeping measurement loadings/thresholds fixed.
#'   Example: \code{c("factor_var_1", "factor_var_2")} to free factor variances.
#'   Only used when previous_stage is provided.
#'
#' @return An object of class "model_system". A list of model_component objects and one factor_model object.
#' @examples
#' dat <- data.frame(y1 = rnorm(50), y2 = rnorm(50), intercept = 1)
#' fm <- define_factor_model(n_factors = 1)
#' mc1 <- define_model_component("m1", dat, "y1", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = 1)
#' mc2 <- define_model_component("m2", dat, "y2", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = NA_real_)
#' ms <- define_model_system(components = list(mc1, mc2), factor = fm)
#' @export
define_model_system <- function(components, factor, previous_stage = NULL, weights = NULL,
                                equality_constraints = NULL, free_params = NULL) {
  # Validate the inputs:

  if (!is.list(components) || !all(sapply(components, inherits, "model_component"))) {
    stop("Input must be a list of model_component objects.")
  }

  # Auto-name components from their internal name field if list is unnamed
  if (is.null(names(components))) {
    comp_names <- vapply(components, function(c) c$name, character(1))
    if (any(is.na(comp_names) | comp_names == "")) {
      stop("model_component objects must have names (either via list names or internal 'name' field)")
    }
    names(components) <- comp_names
  }

  if (!inherits(factor, "factor_model")) {
    stop("`factor` must be of class 'factor_model'")
  }


  if (!all(sapply(components, function(x) inherits(x, "model_component")))) {
    stop("All elements must be of class 'model_component'")
  }

  #check if each component has the same factor model as factor

 if (!all(vapply(components, function(comp) identical(get_factor(comp), factor), logical(1)))) {
   stop("All model_components must have Factor as the factor_model")
 }

  # Validate dynamic components
  for (comp in components) {
    if (isTRUE(comp$is_dynamic)) {
      if (is.null(comp$outcome_factor) ||
          comp$outcome_factor < 1L ||
          comp$outcome_factor > factor$n_factors) {
        stop(sprintf("Dynamic component '%s' has invalid outcome_factor (%s). Must be between 1 and %d.",
                     comp$name, comp$outcome_factor, factor$n_factors))
      }
    }
  }

  # Validate weights parameter
  if (!is.null(weights)) {
    if (!is.character(weights) || length(weights) != 1) {
      stop("`weights` must be a single character string (variable name in data)")
    }
  }

  # Validate: mixtures and factor_covariates are mutually exclusive
  if (factor$n_mixtures > 1 && !is.null(factor$factor_covariates)) {
    stop("Mixture models (n_mixtures > 1) cannot be combined with factor_covariates. ",
         "These features are mutually exclusive.")
  }

  # Validate equality_constraints parameter
  if (!is.null(equality_constraints)) {
    if (!is.list(equality_constraints)) {
      stop("`equality_constraints` must be a list of character vectors")
    }
    for (i in seq_along(equality_constraints)) {
      constraint <- equality_constraints[[i]]
      if (!is.character(constraint) || length(constraint) < 2) {
        stop(sprintf("Each equality constraint must be a character vector with at least 2 parameter names (constraint %d)", i))
      }
    }
  }

  # Validate free_params
  if (!is.null(free_params) && is.null(previous_stage)) {
    warning("`free_params` is ignored when `previous_stage` is not provided")
    free_params <- NULL
  }
  if (!is.null(free_params)) {
    if (!is.character(free_params)) {
      stop("`free_params` must be a character vector of parameter names")
    }
  }

  # Handle previous_stage if provided
  previous_stage_info <- NULL
  if (!is.null(previous_stage)) {
    # Validate previous_stage
    if (!is.list(previous_stage) ||
        !all(c("model_system", "estimates", "std_errors") %in% names(previous_stage))) {
      stop("`previous_stage` must be a result object from estimate_model_rcpp() with model_system, estimates, and std_errors")
    }

    prev_ms <- previous_stage$model_system
    if (!inherits(prev_ms, "model_system")) {
      stop("`previous_stage$model_system` must be a model_system object")
    }

    # Check factor models match
    # For SE_linear/SE_quadratic Stage 2: allow different factor structures
    # as long as n_factors matches. This enables 2-stage estimation where
    # Stage 1 uses independent factors and Stage 2 uses SE structure.
    se_structures <- c("SE_linear", "SE_quadratic")
    allow_different_structure <- FALSE

    if (!identical(prev_ms$factor, factor)) {
      # Check if this is an SE model with matching n_factors
      if (factor$factor_structure %in% se_structures &&
          prev_ms$factor$n_factors == factor$n_factors) {
        allow_different_structure <- TRUE
        message(sprintf("Two-stage SE estimation: Stage 1 (%s) -> Stage 2 (%s)",
                        prev_ms$factor$factor_structure, factor$factor_structure))
        message("  Measurement parameters (loadings, thresholds) will be fixed from Stage 1")
        message("  Factor distribution parameters will use Stage 2 structure")
      } else {
        stop("Factor model in previous_stage must be identical to current factor model, ",
             "unless Stage 2 uses SE_linear or SE_quadratic with matching n_factors")
      }
    }

    # Validate free_params are in previous_stage estimates
    if (!is.null(free_params)) {
      prev_param_names <- names(previous_stage$estimates)
      invalid_params <- setdiff(free_params, prev_param_names)
      if (length(invalid_params) > 0) {
        stop("free_params contains parameter names not in previous_stage: ",
             paste(invalid_params, collapse = ", "))
      }
    }

    # Mark all previous stage components as having fixed parameters
    prev_components <- prev_ms$components
    for (i in seq_along(prev_components)) {
      prev_components[[i]]$all_params_fixed <- TRUE
    }

    # Validate: no new component may share a name with a previous-stage component.
    # Duplicate names cause parameter-matching failures: match() finds the first
    # occurrence (the fixed Stage 1 copy) and silently skips the second (the user's
    # free copy), leaving those parameters unfixed and producing wild estimates.
    prev_names <- vapply(prev_components, function(c) c$name, character(1))
    new_names <- vapply(components, function(c) c$name, character(1))
    dup_names <- intersect(prev_names, new_names)
    if (length(dup_names) > 0) {
      stop(sprintf(
        paste0("Component name(s) duplicated between previous_stage and new components: %s. ",
               "When using previous_stage, only pass NEW components (e.g., the outcome model). ",
               "The measurement components from Stage 1 are prepended automatically."),
        paste(dup_names, collapse = ", ")))
    }

    # Prepend previous stage components to new components
    components <- c(prev_components, components)

    # Determine which parameters to fix (all except free_params)
    all_param_names <- names(previous_stage$estimates)

    # For SE structure Stage 2: only fix measurement parameters, not factor distribution
    if (allow_different_structure) {
      # Identify factor distribution parameters to exclude from fixing.
      # These include: factor_var_*, se_*, chol_* (correlation params),
      # factor_mean_* (mean covariates), typeprob_*_intercept and
      # type_*_loading_* (type-mixture parameters that describe the factor
      # distribution over types). The type-related patterns only apply to
      # factor-level parameters, not to component-level type intercepts,
      # which have the form "{comp_name}_type_*_intercept" and therefore do
      # not begin with "type_" or "typeprob_".
      factor_dist_patterns <- c("^factor_var", "^se_", "^chol_",
                                "^factor_mean_", "^typeprob_",
                                "^type_[0-9]+_loading_")
      factor_dist_params <- unlist(lapply(factor_dist_patterns, function(p) {
        grep(p, all_param_names, value = TRUE)
      }))

      # Only fix measurement parameters (loadings, thresholds, sigmas, intercepts, betas)
      measurement_params <- setdiff(all_param_names, factor_dist_params)

      # Apply free_params override if specified
      fixed_param_names <- if (!is.null(free_params)) {
        setdiff(measurement_params, free_params)
      } else {
        measurement_params
      }

      message(sprintf("  Fixing %d measurement parameters, ignoring %d factor distribution parameters",
                      length(fixed_param_names), length(factor_dist_params)))
    } else {
      # Standard behavior: fix all except free_params
      fixed_param_names <- if (!is.null(free_params)) {
        setdiff(all_param_names, free_params)
      } else {
        all_param_names
      }
    }

    # Store metadata about previous stage
    previous_stage_info <- list(
      n_components = length(prev_components),
      fixed_param_values = previous_stage$estimates[fixed_param_names],
      fixed_param_names = fixed_param_names,
      free_param_names = free_params,
      fixed_std_errors = previous_stage$std_errors[fixed_param_names],
      n_params_fixed = length(fixed_param_names),
      all_param_values = previous_stage$estimates,  # Keep all for initialization
      allow_different_structure = allow_different_structure
    )

    # Check if factor variance should be fixed
    # For SE structure Stage 2: factor variances are always free (use Stage 2 structure)
    if (allow_different_structure) {
      factor$variance_fixed <- FALSE
    } else {
      # Standard behavior
      # Fix variance only if it's NOT in free_params
      var_param_names <- grep("^factor_var", all_param_names, value = TRUE)
      free_var_params <- intersect(var_param_names, free_params)

      if (length(free_var_params) == 0) {
        # All factor variances are fixed
        factor$variance_fixed <- TRUE
        factor$variance_value <- previous_stage$estimates[var_param_names[1]]
      } else if (length(free_var_params) < length(var_param_names)) {
        # Some variances fixed, some free - complex case
        # For now, mark which specific variances are fixed
        factor$variance_fixed <- FALSE
        factor$variance_partially_fixed <- TRUE
        factor$fixed_variance_params <- setdiff(var_param_names, free_var_params)
        factor$fixed_variance_values <- previous_stage$estimates[factor$fixed_variance_params]
      } else {
        # All factor variances are free
        factor$variance_fixed <- FALSE
      }
    }
  }

  out <- list(
    components = components,
    factor = factor,
    previous_stage_info = previous_stage_info,
    weights = weights,
    equality_constraints = equality_constraints)
  class(out) <- "model_system"
  return(out)
}

#' @export
print.model_system <- function(x, ...) {

  stopifnot(inherits(x, "model_system")) #check that x is a class model_system
  comps <- x$components
  n <- length(comps)
  nms <- names(comps) #captures the list names (check if there are errors in doing it this way)

  cat("Model System\n")
  cat("------------\n")
  cat("Components:", n, "\n")
  if (!is.null(x$weights)) {
    cat("Observation weights:", x$weights, "\n")
  }
  if (!is.null(x$equality_constraints)) {
    cat("Equality constraints:", length(x$equality_constraints), "group(s)\n")
    for (i in seq_along(x$equality_constraints)) {
      cat(sprintf("  [%d] %s\n", i, paste(x$equality_constraints[[i]], collapse = " = ")))
    }
  }

  if (!n) return(invisible(x))

  for (i in seq_along(comps)) { #loops over indices, seq_along safer for empty vectors
    # header for this component
    label <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("<unnamed-", i, ">") #get label for component, otherwise call it unnamed - number
    cat("\n[", i, "] ", label, "\n", sep = "") #print section header
    cat(strrep("-", 2 + nchar(label) + nchar(as.character(i))), "\n", sep = "") #print underline

    # delegate to the component's own print method
    print(comps[[i]], ...)   # passes ... along to the model_component's own print method
  }

  invisible(x) #if 0 components, stop printing and return x invisibly
}

