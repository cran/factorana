#' Fix a factor-distribution parameter at model-definition time
#'
#' Constrains a parameter in the latent-factor distribution (factor variance,
#' SE-equation slope / intercept / residual variance, type-probability
#' intercept, type-probability loading, factor-mean covariate coefficient,
#' or SE covariate coefficient) to a fixed value. The parameter will be held
#' at that value during estimation, excluded from the optimizer's free vector,
#' and reported in `result$estimates` with `result$std_errors == 0` and
#' `result$param_table$fixed == TRUE`.
#'
#' Use this for substantive restrictions known at model-definition time, e.g.,
#' setting a `type_<t>_loading_<k>` to 0 when factor `k` is known not to enter
#' the type-probability model. Estimating such a loading free can leak
#' identification, push the optimizer along a flat ridge, and inflate the
#' standard errors of the remaining type-model parameters.
#'
#' @section Conflicts with `previous_stage` / `free_params`:
#'   When `define_model_system(previous_stage = ..., free_params = ...)` is
#'   used, `fix_factor_param()` always wins. The user-fixed value is the
#'   binding constraint; any value for the same name in `previous_stage$
#'   estimates` is ignored, and inclusion of the same name in `free_params`
#'   is silently honoured (the parameter stays fixed regardless). A warning
#'   is emitted once per conflict so accidental misuse surfaces during
#'   debugging without spamming routine workflows.
#'
#' @section Auto-fixed slots:
#'   For SE structures, `type_<t>_loading_<n_factors>` (the type loading on
#'   the outcome factor) is auto-fixed at 0 because the type probability
#'   model must depend only on input factors. Calling
#'   `fix_factor_param(fm, "type_<t>_loading_<n_factors>", 0)` is a no-op.
#'   Calling it with any non-zero value is an error.
#'
#'   When a factor variance is non-identified (no measurement component
#'   fixes a non-zero loading on that factor), it is normally auto-fixed at
#'   its initial value of 1. An explicit
#'   `fix_factor_param(fm, "factor_var_<k>", v)` overrides the auto-fix and
#'   uses the user-supplied value.
#'
#' @param factor_model A `factor_model` object from `define_factor_model()`.
#' @param name Either a single character string naming the factor-distribution
#'   parameter to fix, or a character vector of names paired element-wise with
#'   `value`, or a named numeric vector whose names are the parameter names
#'   and whose values are the fixed values (in which case `value` must be
#'   omitted).
#' @param value Numeric. The value to fix the parameter to. Pass `NA` to
#'   release a previously fixed parameter (unfix). Length must match `name`
#'   (length 1 is recycled).
#'
#' @return The modified `factor_model` object. The fix is stored in
#'   `factor_model$fixed_params` as a named numeric vector.
#'
#' @examples
#' fm <- define_factor_model(n_factors = 3, n_types = 2,
#'                           factor_structure = "SE_linear")
#'
#' # Single fix
#' fm <- fix_factor_param(fm, "type_2_loading_2", 0.0)
#'
#' # Batch via named numeric
#' fm <- fix_factor_param(fm, c(type_2_loading_2 = 0.0,
#'                              type_2_loading_3 = 0.0))
#'
#' # Unfix
#' fm <- fix_factor_param(fm, "type_2_loading_2", NA_real_)
#'
#' @export
fix_factor_param <- function(factor_model, name, value = NULL) {

  if (!inherits(factor_model, "factor_model")) {
    stop("`factor_model` must be an object of class 'factor_model' from ",
         "define_factor_model().")
  }

  # Resolve the (name, value) pairs from either of the two input shapes.
  if (is.null(value)) {
    if (!is.numeric(name) || is.null(names(name))) {
      stop("When `value` is omitted, `name` must be a NAMED numeric vector ",
           "(e.g., c(type_2_loading_2 = 0)).")
    }
    pairs_names <- names(name)
    pairs_values <- as.numeric(name)
  } else {
    if (!is.character(name) || length(name) == 0L) {
      stop("`name` must be a non-empty character vector of factor-",
           "distribution parameter names.")
    }
    if (!is.numeric(value)) {
      stop("`value` must be numeric (use NA to unfix).")
    }
    if (length(value) == 1L && length(name) > 1L) {
      value <- rep(value, length(name))
    }
    if (length(name) != length(value)) {
      stop("`name` and `value` must have the same length (or `value` of ",
           "length 1 to recycle).")
    }
    pairs_names <- name
    pairs_values <- as.numeric(value)
  }

  if (anyDuplicated(pairs_names)) {
    dups <- unique(pairs_names[duplicated(pairs_names)])
    stop("Duplicate parameter name(s) in `name`: ",
         paste(dups, collapse = ", "))
  }

  valid_names <- .factor_param_valid_names(factor_model)

  for (i in seq_along(pairs_names)) {
    nm <- pairs_names[i]
    val <- pairs_values[i]

    if (!nzchar(nm)) {
      stop("Parameter name must be a non-empty string.")
    }
    if (!(nm %in% valid_names)) {
      stop(sprintf(
        "'%s' is not a valid factor-distribution parameter name for this ",
        nm),
        sprintf(
          "factor model (n_factors=%d, n_types=%d, factor_structure='%s'%s%s).\n",
          factor_model$n_factors, factor_model$n_types,
          factor_model$factor_structure,
          if (!is.null(factor_model$factor_covariates))
            sprintf(", factor_covariates=c(%s)",
                    paste(sprintf('"%s"', factor_model$factor_covariates),
                          collapse = ",")) else "",
          if (!is.null(factor_model$se_covariates))
            sprintf(", se_covariates=c(%s)",
                    paste(sprintf('"%s"', factor_model$se_covariates),
                          collapse = ",")) else ""
        ),
        "Valid names: ",
        paste(valid_names, collapse = ", "))
    }

    # Auto-fixed outcome-factor type loading: idempotent at 0, error otherwise.
    fs <- factor_model$factor_structure
    if (!is.null(fs) && fs %in% c("SE_linear", "SE_quadratic") &&
        factor_model$n_types > 1L) {
      m <- regmatches(nm, regexec("^type_([0-9]+)_loading_([0-9]+)$", nm))[[1]]
      if (length(m) >= 3) {
        k_idx <- as.integer(m[3])
        if (k_idx == factor_model$n_factors && !is.na(val) && abs(val) > 1e-12) {
          stop(sprintf(
            "fix_factor_param(): '%s' is the type loading on the outcome ",
            nm),
            sprintf(
              "factor of an %s model and is auto-fixed at 0 (type ",
              fs),
            "probabilities must depend only on input factors). ",
            sprintf("Cannot fix it at %g; only 0 (or no fix) is allowed.", val))
        }
      }
    }

    if (is.na(val)) {
      # Unfix: drop the entry if present.
      if (!is.null(factor_model$fixed_params) &&
          nm %in% names(factor_model$fixed_params)) {
        factor_model$fixed_params <-
          factor_model$fixed_params[setdiff(names(factor_model$fixed_params), nm)]
        if (length(factor_model$fixed_params) == 0L) {
          factor_model$fixed_params <- NULL
        }
      }
    } else {
      # Fix or replace.
      if (is.null(factor_model$fixed_params)) {
        factor_model$fixed_params <- numeric(0)
      }
      factor_model$fixed_params[nm] <- val
    }
  }

  factor_model
}


# Internal: enumerate all factor-distribution parameter names that
# build_parameter_metadata() would generate for this factor model.
# MUST stay in sync with the layout in optimize_model.R::
# build_parameter_metadata() and FactorModel.cpp + rcpp_interface.cpp.
.factor_param_valid_names <- function(fm) {
  n_factors <- as.integer(fm$n_factors)
  n_types <- if (is.null(fm$n_types)) 1L else as.integer(fm$n_types)
  n_mixtures <- if (is.null(fm$n_mixtures)) 1L else as.integer(fm$n_mixtures)
  fs <- if (is.null(fm$factor_structure)) "independent" else fm$factor_structure
  is_se <- fs %in% c("SE_linear", "SE_quadratic")
  n_var_factors <- if (is_se) (n_factors - 1L) else n_factors

  out <- character(0)

  # Block 1: factor variances (per mixture)
  for (m in seq_len(n_mixtures)) {
    for (k in seq_len(n_var_factors)) {
      out <- c(out,
               if (n_mixtures == 1L) sprintf("factor_var_%d", k)
               else sprintf("mix%d_factor_var_%d", m, k))
    }
  }

  # Correlation parameter (correlation structure only, n_factors == 2)
  if (fs == "correlation" && n_factors == 2L) {
    out <- c(out, "factor_corr_1_2")
  }

  # Block 2: mixture means + log-weights (n_mixtures > 1)
  if (n_mixtures > 1L) {
    for (m in seq_len(n_mixtures - 1L)) {
      for (k in seq_len(n_var_factors)) {
        out <- c(out, sprintf("mix%d_factor_mean_%d", m, k))
      }
    }
    for (m in seq_len(n_mixtures - 1L)) {
      out <- c(out, sprintf("mix%d_logweight", m))
    }
  }

  # Block 3: SE parameters (SE structures only)
  if (is_se) {
    out <- c(out, "se_intercept")
    for (k in seq_len(n_var_factors)) {
      out <- c(out, sprintf("se_linear_%d", k))
    }
    if (fs == "SE_quadratic") {
      for (k in seq_len(n_var_factors)) {
        out <- c(out, sprintf("se_quadratic_%d", k))
      }
    }
    if (n_types > 1L) {
      for (t in 2:n_types) {
        out <- c(out, sprintf("se_intercept_type_%d", t))
      }
    }
    out <- c(out, "se_residual_var")
  }

  # Block 4: type-probability params (typeprob intercepts + type loadings)
  if (n_types > 1L) {
    for (t in 2:n_types) {
      out <- c(out, sprintf("typeprob_%d_intercept", t))
    }
    for (t in 2:n_types) {
      for (k in seq_len(n_factors)) {
        out <- c(out, sprintf("type_%d_loading_%d", t, k))
      }
    }
  }

  # Block 5: factor-mean covariate params
  if (!is.null(fm$factor_covariates) && length(fm$factor_covariates) > 0L) {
    n_fac_with_mean <- if (is_se) (n_factors - 1L) else n_factors
    for (k in seq_len(n_fac_with_mean)) {
      for (cn in fm$factor_covariates) {
        out <- c(out, sprintf("factor_mean_%d_%s", k, cn))
      }
    }
  }

  # Block 6: SE covariate params
  if (!is.null(fm$se_covariates) && length(fm$se_covariates) > 0L) {
    for (cn in fm$se_covariates) {
      out <- c(out, sprintf("se_cov_%s", cn))
    }
  }

  out
}
