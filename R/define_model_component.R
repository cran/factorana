#' Define a model component
#'
#' @param name name of model component ####Long name and s name (like in the C++)
#' @param data data.frame for validation
#' @param outcome Character. Name of the outcome variable.
#' @param factor model. Object of Class factor_model.
#' @param evaluation_indicator Character (optional). Variable used for evaluation subsample.
#' @param covariates List of Character vectors. Names of covariates.
#' @param model_type Character. Type of model (e.g., "linear", "logit", "probit").
#' @param intercept Logical. Whether to include an intercept (default = TRUE).
#' @param num_choices Integer. Number of choices (for multinomial models).
#' @param nrank Integer (optional). Rank for exploded multinomial logit.
#' @param exclude_chosen Logical. For exploded logit, whether to exclude already-chosen
#'   alternatives from later ranks (default TRUE). Set to FALSE for exploded nested logit
#'   where the same nest can be chosen multiple times.
#' @param rankshare_var Character (optional). Column name (or prefix) for rank-share
#'   correction variables. These provide rank-and-choice-specific adjustments to the
#'   linear predictor. Data layout: (num_choices-1) * nrank columns, accessed as
#'   rankshare_var + (num_choices-1)*irank + icat for irank=0..nrank-1, icat=0..num_choices-2.
#' @param loading_normalization Numeric vector of length `n_factors` (optional).
#'   Component-specific loading constraints. Overrides factor model normalization.
#'   - `NA` --> loading is free (estimated).
#'   - numeric value --> loading is fixed at that value (e.g. `1` for identification).
#'   If NULL, uses the factor model's default normalization.
#' @param factor_spec Character. Specification for factor terms in linear predictor.
#'   - `"linear"` (default): Only linear factor terms (lambda * f)
#'   - `"quadratic"`: Linear + quadratic terms (lambda * f + lambda_quad * f^2)
#'   - `"interactions"`: Linear + interaction terms (lambda * f + lambda_inter * f_j * f_k)
#'   - `"full"`: Linear + quadratic + interaction terms
#'   Note: Interaction terms require n_factors >= 2.
#' @param use_types Logical. Whether this component uses type-specific intercepts
#'   when n_types > 1 in the factor model. Default FALSE. When TRUE and n_types > 1,
#'   the component will have (n_types - 1) type-specific intercept parameters that
#'   shift the linear predictor for each non-reference type. This allows types to
#'   affect outcome models while keeping measurement models type-invariant.
#' @param skip_collinearity_check Logical. If TRUE, skip the multicollinearity check
#'   on the design matrix. Useful when many coefficients will be fixed via fix_coefficient()
#'   after component creation, which resolves the collinearity. Default FALSE.
#'
#' @return An object of class "model_component". A list representing the model component
#' @examples
#' dat <- data.frame(y = rnorm(50), intercept = 1)
#' fm <- define_factor_model(n_factors = 1)
#' mc <- define_model_component("Y", dat, "y", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = 1)
#' @export

define_model_component <- function(name,
                                   data,
                                   outcome,
                                   factor,
                                   evaluation_indicator = NULL,
                                   covariates,
                                   model_type = c("linear", "logit", "probit", "oprobit"),
                                   intercept = TRUE,
                                   num_choices = 2,
                                   nrank = NULL,
                                   exclude_chosen = TRUE,
                                   rankshare_var = NULL,
                                   loading_normalization = NULL,
                                   factor_spec = c("linear", "quadratic", "interactions", "full"),
                                   use_types = FALSE,
                                   skip_collinearity_check = FALSE) {
  # ---- 1. Basic argument checks ----
  # Confirm data types and presence of required columns/objects

  # Accept data.frame, data.table, or tibble - convert to plain data.frame for consistent behavior
  if (!is.data.frame(data)) stop("`data` must be a data.frame (or data.table/tibble).")
  data <- as.data.frame(data)

  # Validate outcome - can be single column name or vector of column names (for exploded logit)
  if (!is.character(outcome) || length(outcome) < 1L) {
    stop("`outcome` must be a character vector of column name(s).")
  }
  missing_outcomes <- setdiff(outcome, names(data))
  if (length(missing_outcomes) > 0) {
    stop("Outcome variable(s) not found in data: ", paste(missing_outcomes, collapse = ", "))
  }

  #is the factor model a list?
  if (!is.list(factor)) stop("`factor` must be a list.") #should it be a list or object of S3 class???

  # factor model class + derive k and normalization
  if (!inherits(factor, "factor_model")) {
    stop("`factor` must be an object of class 'factor_model'.")
  }

  # ---- 2. Validate factor model and normalization ----
  # Pull number of factors and ensure valid loading_normalization setup

  k <- as.integer(factor$n_factors)
  if (is.na(k) || k < 0L) stop("factor$n_factors must be a non-negative integer")

  # Component-specific normalization (defaults to all free)
  if (is.null(loading_normalization)) {
    # Default: all loadings are free (NA) if k > 0, empty vector if k = 0
    loading_normalization <- if (k > 0) rep(NA_real_, k) else numeric(0)
  } else {
    # Validate component-specific normalization
    if (!is.numeric(loading_normalization) ||
        length(loading_normalization) != k) {
      if (k == 0) {
        stop("`loading_normalization` must be numeric(0) when n_factors = 0.")
      } else {
        stop("`loading_normalization` must be numeric and length ", k, ".")
      }
    }
  }

  # ---- 3. Validate evaluation indicator ----
  # Ensure indicator is a valid single logical/numeric variable if provided

  if (!is.null(evaluation_indicator)) {
    if (!is.character(evaluation_indicator) || length(evaluation_indicator) != 1L)
      stop("`evaluation_indicator` must be a single column name (character).")
    if (!(evaluation_indicator %in% names(data)))
      stop("`evaluation_indicator` not found in `data`.")
  }

  # ---- 4. Validate covariates ----
  # Check all covariate names exist and are non-missing
  # NULL or character(0) is allowed (factor-only model)

  if (!is.null(covariates)) {
    if (!is.character(covariates) || length(covariates) < 1L) {
      stop("`covariates` must be NULL or a non-empty character vector of column names.")
    }

    missing <- setdiff(covariates, names(data))
    if (length(missing)) {
      stop("Covariates not found in `data`: ", paste(missing, collapse = ", "))
    }
  } else {
    # NULL covariates: set to empty character vector for consistency
    covariates <- character(0)
  }

  # ---- 5. Validate model configuration arguments ----
  # model_type, intercept, num_choices, nrank, factor_spec
  model_type <- match.arg(model_type)
  factor_spec <- match.arg(factor_spec)

  # Validate factor_spec for interactions (requires k >= 2)
  if (factor_spec %in% c("interactions", "full") && k < 2) {
    warning("factor_spec='", factor_spec, "' requires n_factors >= 2. ",
            "Downgrading to '", if (factor_spec == "full") "quadratic" else "linear", "'.")
    factor_spec <- if (factor_spec == "full") "quadratic" else "linear"
  }

  # Validate use_types
  if (!is.logical(use_types) || length(use_types) != 1L) {
    stop("`use_types` must be a single TRUE/FALSE.")
  }
  # Get n_types from factor model
  n_types <- factor$n_types
  if (is.null(n_types)) n_types <- 1L
  if (use_types && n_types < 2L) {
    warning("`use_types = TRUE` has no effect when n_types < 2. Setting to FALSE.")
    use_types <- FALSE
  }

  # intercept: is it a boolean?
  if (!is.logical(intercept) || length(intercept) != 1L) {
    stop("`intercept` must be a single TRUE/FALSE.")
  }

  # ---- Intercept sanity checks ----
  has_intercept_in_covariates <- length(covariates) > 0 &&
    any(covariates %in% c("intercept", "constant"))

  # Check 1: ERROR if an intercept covariate is included in an ordered probit model.
  # Ordered probit absorbs the intercept into the cut points, so a separate
  # intercept covariate is not identified. This silently produces NaN parameters
  # and a singular Hessian, wasting hours of computation.
  if (model_type == "oprobit" && has_intercept_in_covariates) {
    stop(sprintf(
      paste0("Component '%s': ordered probit (oprobit) models must NOT include an ",
             "intercept covariate (found '%s' in covariates). The intercept is absorbed ",
             "into the cutpoint thresholds. Remove the intercept from covariates."),
      name, intersect(covariates, c("intercept", "constant"))[1]))
  }

  # Check 2: WARNING if a linear/probit/logit model has no intercept and the
  # outcome is far from zero-mean. Without an intercept the model is
  # misspecified for non-centered outcomes, causing degenerate convergence.
  # Only fire when the caller has NOT explicitly opted out of an intercept
  # via `intercept = FALSE`: the explicit opt-out is a deliberate choice
  # (common in validation and shape-only tests), not the accidental
  # omission this warning is designed to catch.
  if (isTRUE(intercept) &&
      model_type %in% c("linear", "probit", "logit") &&
      !has_intercept_in_covariates) {
    # Compute outcome statistics on the FULL data (before eval_indicator subsetting)
    y_raw <- data[[outcome[1]]]
    y_raw <- y_raw[is.finite(y_raw)]
    if (length(y_raw) > 10) {
      y_mean <- mean(y_raw)
      y_sd <- sd(y_raw)
      if (y_sd > 0 && abs(y_mean) > 0.1 * y_sd) {
        warning(sprintf(
          paste0("Component '%s': no intercept covariate found for %s model, but ",
                 "outcome '%s' has non-zero mean (%.3g, sd=%.3g). Without an intercept, ",
                 "the model may converge to a degenerate point. Consider adding an ",
                 "intercept: data$intercept <- 1; covariates = c('intercept', ...)."),
          name, model_type, outcome[1], y_mean, y_sd),
          immediate. = TRUE)
      }
    }
  }

  # Warn if intercept=TRUE for models where it has no effect
  # Note: intercept=TRUE is just metadata. Users must add their own intercept column to data
  # and include it in covariates if they need an intercept term.
  if (isTRUE(intercept)) {
    if (model_type == "oprobit") {
      # Already handled above as an error if intercept IS in covariates.
      # This warning is for the intercept=TRUE FLAG, not the covariate.
    } else if (model_type == "linear" && !has_intercept_in_covariates) {
      warning("`intercept = TRUE` does not automatically add an intercept. ",
              "Add a column of 1s to your data (e.g., data$intercept <- 1) ",
              "and include 'intercept' in covariates.")
    } else if (model_type == "probit" && !has_intercept_in_covariates) {
      warning("`intercept = TRUE` does not automatically add an intercept for probit models. ",
              "Add a column of 1s to your data and include 'intercept' in covariates if needed.")
    }
  }

  # num_choices: basic validation (detailed validation against data comes later)
  num_choices <- as.integer(num_choices)
  if (is.na(num_choices) || num_choices < 2L) {
    stop("`num_choices` must be an integer >= 2.")
  }
  if (num_choices > 50L) {
    stop("`num_choices` cannot exceed 50. Found: ", num_choices)
  }

  # n_rank: infer from outcome vector length, or validate explicit parameter
  # For exploded logit, nrank is the number of ranked choices
  if (length(outcome) > 1) {
    # Multiple outcome columns = exploded logit
    if (model_type != "logit") {
      stop("Multiple outcome variables (exploded logit) only supported for model_type='logit'.")
    }
    inferred_nrank <- length(outcome)
    if (!is.null(nrank) && nrank != inferred_nrank) {
      warning("nrank (", nrank, ") does not match length(outcome) (", inferred_nrank, "). Using length(outcome).")
    }
    nrank <- inferred_nrank
  } else if (!is.null(nrank)) {
    nrank <- as.integer(nrank)
    if (!is.finite(nrank) || nrank < 1L) stop("`nrank` must be a positive integer when provided.")
  } else {
    nrank <- 1L  # Default: single outcome (standard logit)
  }

  # Validate exclude_chosen (for exploded logit)
  if (!is.logical(exclude_chosen) || length(exclude_chosen) != 1L) {
    stop("`exclude_chosen` must be a single TRUE/FALSE.")
  }
  if (!exclude_chosen && nrank == 1L) {
    warning("`exclude_chosen=FALSE` has no effect for standard logit (nrank=1).")
  }

  # Validate rankshare_var (for exploded nested logit)
  if (!is.null(rankshare_var)) {
    if (!is.character(rankshare_var) || length(rankshare_var) != 1L) {
      stop("`rankshare_var` must be a single column name (character).")
    }
    if (!(rankshare_var %in% names(data))) {
      stop("`rankshare_var` '", rankshare_var, "' not found in data.")
    }
    if (nrank == 1L) {
      warning("`rankshare_var` has no effect for standard logit (nrank=1).")
    }
    if (model_type != "logit") {
      stop("`rankshare_var` is only supported for model_type='logit'.")
    }
  }

  # ---- 6. Evaluation subset conditioning ----
  # Restrict data to rows where evaluation_indicator = TRUE/1.
  # Panel data with wave- or cohort-specific item availability often has
  # evaluation_indicator that is entirely zero for some (component, wave)
  # combinations. In that case the component contributes nothing to the
  # likelihood and we warn-and-skip rather than erroring: returning NULL
  # so define_model_system() can drop it. Only the empty-input case
  # (data already has zero rows BEFORE filtering) remains a hard error,
  # since that is almost always a caller mistake.

  if (nrow(data) == 0L) stop("Evaluation subset has zero rows")

  idx <- rep(TRUE, nrow(data))  # default: check all rows
  if (!is.null(evaluation_indicator)) {
    ei <- data[[evaluation_indicator]]
    if (is.logical(ei)) {
      idx <- !is.na(ei) & ei
    } else if (is.numeric(ei) || is.integer(ei)) {
      idx <- !is.na(ei) & (ei == 1L)
    } else {
      stop("`evaluation_indicator` must be logical or 0/1 numeric.")
    }
    if (sum(idx) == 0L) {
      warning(sprintf(
        "Component '%s' is skipped: evaluation_indicator '%s' has no TRUE/1 rows.",
        name, evaluation_indicator),
        call. = FALSE)
      return(invisible(NULL))
    }
  }

  data <- data[idx, , drop = FALSE]
  rownames(data) <- NULL

  idx <- rep(TRUE, nrow(data))

  # ---- 7. Ordered probit handling ----
  # Convert numeric or unordered factor outcomes to ordered factor with >= 3 categories
  if (model_type == "oprobit") {
    y_sub <- data[[outcome]]

    if (is.factor(y_sub)) {
      # Make sure it is 'ordered'
      if (!is.ordered(y_sub)) y_sub <- ordered(y_sub)
    } else {
      # Accept integer-like labels (1..J or 0..J-1), otherwise error
      if (!is.numeric(y_sub) && !is.integer(y_sub))
        stop("oprobit outcome must be integer-like or an ordered factor.")
      # map to contiguous 1..J and mark ordered
      u   <- sort(unique(na.omit(as.integer(y_sub))))
      map <- match(as.integer(y_sub), u)
      y_sub <- ordered(map)
    }

    if (nlevels(y_sub) < 3L)
      stop("Ordered probit requires an outcome with >= 3 ordered categories.")

    data[[outcome]] <- y_sub
  }

  # ---- 8. Missing value checks ----
  # Error if NAs found in outcome or covariates (users must handle missing data explicitly)
  # Exception: For exploded logit, missing values in rank outcomes indicate unused ranks (allowed)

  if (length(outcome) == 1) {
    # Standard single outcome
    if (anyNA(data[[outcome]][idx])) {
      stop("Missing values in outcome variable within evaluation subset.")
    }
  }
  # For exploded logit (length(outcome) > 1), missing values in rank columns are allowed
  # They indicate the individual didn't use that rank

  # covariates must not have missing (on the same subset)
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (anyNA(data[[cov]][idx])) {
        stop("Missing values found in covariate: ", cov)
      }
    }
  }


  # ---- 8b. Multicollinearity check ----
  # Build design matrix and check for rank deficiency and near-collinearity

  if (length(covariates) > 0 && !skip_collinearity_check) {
    X <- as.matrix(data[idx, covariates, drop = FALSE])
    col_names <- covariates

    expected_rank <- ncol(X)
    actual_rank <- qr(X)$rank

    if (actual_rank < expected_rank) {
      # Find which columns are linearly dependent
      # Use QR decomposition with pivoting to identify problematic columns
      qr_result <- qr(X, LAPACK = TRUE)
      pivot_order <- qr_result$pivot
      dependent_cols <- col_names[pivot_order[(actual_rank + 1):expected_rank]]

      stop("Design matrix is rank deficient (multicollinearity detected).\n",
           "  Expected rank: ", expected_rank, ", Actual rank: ", actual_rank, "\n",
           "  Problematic column(s): ", paste(dependent_cols, collapse = ", "), "\n",
           "  Remove or combine collinear variables.")
    }

    # Check for near-collinearity using condition number of X'X
    # Condition number > 30 is often considered problematic
    # Condition number > 1000 indicates severe multicollinearity
    if (nrow(X) > ncol(X)) {
      XtX <- crossprod(X)
      eigenvalues <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values

      # Condition number is ratio of largest to smallest eigenvalue
      # Use absolute values to handle numerical precision issues
      eigenvalues <- abs(eigenvalues)
      if (min(eigenvalues) > .Machine$double.eps) {
        condition_number <- max(eigenvalues) / min(eigenvalues)

        if (condition_number > 1e12) {
          warning("Severe multicollinearity detected (condition number = ",
                  format(condition_number, scientific = TRUE, digits = 2), ").\n",
                  "  This may cause numerical instability in estimation.\n",
                  "  Consider removing or combining highly correlated variables.")
        } else if (condition_number > 1e6) {
          warning("Moderate multicollinearity detected (condition number = ",
                  format(condition_number, scientific = TRUE, digits = 2), ").\n",
                  "  Standard errors may be inflated.")
        }
      }
    }
  }


  # ---- 9. Model-type specific validity checks ----
  # For single outcome, extract the column. For multiple outcomes (exploded logit), use first column for probit check
  y <- if (length(outcome) == 1) data[[outcome]] else data[[outcome[1]]]

  if (model_type == "probit") {
    if (!all(y %in% c(0, 1))) {
      stop("Outcome for probit must be coded 0/1.")
    }
    if (num_choices != 2L) {
      stop("Probit model requires num_choices = 2. Found: ", num_choices)
    }
  }

  if (model_type == "logit") {
    # For exploded logit (nrank > 1), validate all outcome columns
    if (nrank > 1) {
      all_vals <- c()
      for (out_col in outcome) {
        y_rank <- data[[out_col]]
        # For exploded logit, missing ranks can have values outside 1..num_choices (e.g., 0, NA, -1)
        # Only validate non-missing values
        valid_vals <- na.omit(y_rank)
        valid_vals <- valid_vals[valid_vals >= 1 & valid_vals <= num_choices]
        all_vals <- c(all_vals, valid_vals)
      }
      unique_vals <- sort(unique(all_vals))
      if (length(unique_vals) == 0) {
        stop("No valid outcome values found in ranked outcome columns.")
      }
      # For exploded logit, we don't require all choices to appear - just validate range
      if (max(unique_vals) > num_choices) {
        stop("Outcome values exceed num_choices (", num_choices, "). Found max: ", max(unique_vals))
      }
    } else {
      # Standard logit: validate single outcome column
      # C++ expects 1-indexed outcomes (1, 2, ..., K) - NOT 0-indexed
      unique_vals <- sort(unique(na.omit(y)))
      min_val <- min(unique_vals)
      max_val <- max(unique_vals)

      # Check outcomes are 1, 2, ..., K (contiguous 1-indexed integers)
      if (min_val < 1) {
        stop("Logit outcomes must be 1-indexed (1, 2, ..., K). Found minimum value: ", min_val,
             "\n  If your data uses 0/1 coding, add 1 to convert: data$outcome <- data$outcome + 1")
      }
      if (!all(unique_vals == seq(min_val, max_val))) {
        stop("Logit outcomes must be contiguous integers. Found: ", paste(unique_vals, collapse = ", "))
      }

      n_unique <- length(unique_vals)
      if (n_unique > 50L) {
        stop("Logit outcome has too many unique values (", n_unique, "). Maximum supported: 50.")
      }
      if (num_choices != n_unique) {
        stop("num_choices (", num_choices, ") does not match detected unique outcome values (", n_unique, ").")
      }
      # Validate that outcomes start at 1 (for proper C++ indexing)
      if (min_val != 1) {
        stop("Logit outcomes must start at 1. Found minimum: ", min_val,
             "\n  Recode your outcomes to use 1, 2, ..., ", n_unique)
      }
    }
  }


  # ---- 10. Final consistency for ordered probit ----
  # (Recheck contiguous ordered categories)

  n_cats <- NULL
  if (model_type == "oprobit") {
    y_sub <- data[[outcome]][idx]
    # Accept integers 1..J, 0..J-1, or an ordered factor; coerce to ordered factor
    if (is.factor(y_sub)) {
      if (!is.ordered(y_sub)) y_sub <- ordered(y_sub) #coerce to an ordered factor here.
    } else {
      if (!is.numeric(y_sub) && !is.integer(y_sub))
        stop("oprobit outcome must be integer-like or ordered factor.")

      u <- sort(unique(na.omit(as.integer(y_sub))))
      # make contiguous levels starting at 1
      map <- match(as.integer(y_sub), u)
      y_sub <- ordered(map)
    }
    n_cats <- length(levels(y_sub))
    if (n_cats < 2L) stop("oprobit needs at least 2 ordered categories.")
    if (n_cats > 50L) {
      stop("Ordered probit outcome has too many categories (", n_cats, "). Maximum supported: 50.")
    }
    if (num_choices != n_cats) {
      stop("num_choices (", num_choices, ") does not match detected outcome categories (", n_cats, ").")
    }
    # (Optional) replace in data to keep consistency downstream:
    data[[outcome]] <- y_sub
  }

  # ---- 11. Build output object ----
  # Assemble metadata, settings, and derived quantities into list

  # Calculate number of FREE factor loadings (not fixed)
  n_free_loadings <- sum(is.na(loading_normalization))

  # Calculate number of second-order factor loadings
  # Quadratic loadings: one per factor (always free)
  n_quadratic_loadings <- if (factor_spec %in% c("quadratic", "full")) k else 0L

  # Interaction loadings: one per unique pair j < k (always free)
  n_interaction_loadings <- if (factor_spec %in% c("interactions", "full") && k >= 2) {
    as.integer(k * (k - 1) / 2)
  } else {
    0L
  }

  # Total second-order loadings
  n_second_order_loadings <- n_quadratic_loadings + n_interaction_loadings

  # Calculate number of type intercept parameters
  n_type_intercepts <- if (use_types && n_types > 1L) {
    if (model_type == "logit" && num_choices > 2) {
      # Multinomial logit: each non-reference choice gets (n_types - 1) type intercepts
      (num_choices - 1L) * (n_types - 1L)
    } else {
      # Other models: (n_types - 1) type intercepts
      n_types - 1L
    }
  } else {
    0L
  }

  # Calculate number of model parameters based on model type
  if (model_type == "logit" && num_choices > 2) {
    # Multinomial logit: each non-reference choice has its own parameters
    nparamchoice <- length(covariates) + n_free_loadings + n_second_order_loadings
    nparam_model <- (num_choices - 1) * nparamchoice + n_type_intercepts
  } else if (model_type == "oprobit") {
    # Ordered probit: shared coefficients + (num_choices - 1) thresholds + type intercepts
    nparam_model <- length(covariates) + n_free_loadings + n_second_order_loadings + (num_choices - 1) + n_type_intercepts
  } else {
    # Binary models (linear, probit, binary logit)
    nparam_model <- length(covariates) + n_free_loadings + n_second_order_loadings + n_type_intercepts
    if (model_type == "linear") nparam_model <- nparam_model + 1  # Add sigma
  }

  # Number of observations used (after evaluation indicator filtering)
  n_obs <- nrow(data)

  # Data is only needed for validation during component creation.
  # It is not stored to avoid memory duplication (data is passed separately to estimate_model_rcpp).
  out <- list(
    name = name,
    outcome = outcome,
    factor = factor,
    evaluation_indicator = evaluation_indicator,
    covariates = covariates,
    model_type = model_type,
    intercept = intercept,
    num_choices = num_choices,
    nrank = nrank,
    exclude_chosen = exclude_chosen,
    rankshare_var = rankshare_var,
    nparam_model = nparam_model,
    n_obs = n_obs,
    k = k,
    loading = rep(NA_real_, k),
    loading_normalization = loading_normalization,
    factor_spec = factor_spec,
    n_quadratic_loadings = n_quadratic_loadings,
    n_interaction_loadings = n_interaction_loadings,
    use_types = use_types,
    n_type_intercepts = n_type_intercepts,
    fixed_coefficients = list()  # List of fixed coefficient constraints
  )

  class(out) <- "model_component"
  return(out)
}

# ---- S3 methods: getters and printers ----


#' Get component name (internal)
#' @param x Object to extract name from (a \code{model_component}).
#' @param ... Additional arguments (not used).
#' @return A single character string with the component name, or
#'   \code{NA_character_} if the component has no name set.
#' @keywords internal
get_component_name <- function(x, ...) UseMethod("get_component_name")

#' @rdname get_component_name
#' @exportS3Method
get_component_name.model_component <- function(x, ...) {
  nm <- x$name
  if (is.null(nm) || is.na(nm) || !nzchar(nm)) NA_character_ else nm
}

#' Get factor from component (internal)
#' @param x Object to extract factor from (a \code{model_component}).
#' @param ... Additional arguments (not used).
#' @return The \code{factor_model} object attached to the component.
#' @keywords internal
get_factor <- function(x, ...) UseMethod("get_factor")

#' @rdname get_factor
#' @exportS3Method
get_factor.model_component <- function(x, ...) {
  return(x$factor)
}


#' Print method for model_component objects
#'
#' @param x An object of class "model_component".
#' @param ... Not used.
#' @return Invisibly returns \code{x}. Called for its side effect of printing
#'   a human-readable summary of the component to the console.
#' @export
print.model_component <- function(x, ...) {
  cat("Model Component\n")
  cat("------------------------------\n")
  cat("Model:                   ", x$name, "\n")
  cat("Outcome variable:        ", x$outcome, "\n")
  cat("Model type:              ", x$model_type, "\n")
  cat("Intercept:               ", ifelse(x$intercept, "Yes", "No"), "\n")
  cat("Loading normalization:   ",
      paste0(x$loading_normalization, collapse = ", "), "\n")
  cat("Factor specification:    ", x$factor_spec, "\n")
  if (x$n_quadratic_loadings > 0) {
    cat("Quadratic loadings:      ", x$n_quadratic_loadings, "\n")
  }
  if (x$n_interaction_loadings > 0) {
    cat("Interaction loadings:    ", x$n_interaction_loadings, "\n")
  }
  if (!is.null(x$use_types) && x$use_types) {
    cat("Use types:               ", "Yes (", x$n_type_intercepts, " type intercepts)\n", sep = "")
  }
  cat("Number of choices:       ", x$num_choices, "\n")
  if (!is.null(x$nrank) && x$nrank > 1) {
    cat("Rank (nrank):            ", x$nrank, "\n")
    cat("Exclude chosen:          ", ifelse(x$exclude_chosen, "Yes", "No"), "\n")
    if (!is.null(x$rankshare_var)) {
      cat("Rankshare variable:      ", x$rankshare_var, "\n")
    }
  }
  if (!is.null(x$evaluation_indicator)) {
    cat("Evaluation indicator:    ", x$evaluation_indicator, "\n")
  }
  if (length(x$covariates) > 0) {
    cat("Covariates:              ", paste(x$covariates, collapse = ", "), "\n")
  } else {
    cat("Covariates:               (none - factor-only model)\n")
  }

  cat("Observations:            ", x$n_obs, "\n")
  cat("Total parameters:        ", x$nparam_model, "\n")
  invisible(x)
}

