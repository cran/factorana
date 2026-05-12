#' Fix a coefficient in a model component
#'
#' Constrains a regression coefficient (beta) to a fixed value during estimation.
#' The coefficient will not be optimized - it remains at the specified value.
#'
#' @param component A model_component object from define_model_component()
#' @param covariate Character. Name of the covariate whose coefficient to fix.
#'   Must be one of the covariates specified when creating the component.
#' @param value Numeric. The value to fix the coefficient to.
#' @param choice Integer (optional). For multinomial logit models with num_choices > 2,
#'   specifies which choice's coefficient to fix. Must be between 1 and (num_choices - 1).
#'   Choice 0 (reference category) has no parameters.
#'   For other model types, this should be NULL (default).
#'
#' @return The modified model_component object with the fixed coefficient constraint added.
#'
#' @details
#' Fixed coefficients are stored in the component and used during:
#' \itemize{
#'   \item Parameter initialization: fixed values are used directly (no estimation)
#'   \item Optimization: fixed parameters are excluded from the optimization
#'   \item Likelihood evaluation: fixed values are inserted into the parameter vector
#' }
#'
#' Note: This function only fixes regression coefficients (betas), not:
#' \itemize{
#'   \item Factor loadings (use loading_normalization in define_model_component)
#'   \item Sigma (for linear models)
#'   \item Thresholds (for ordered probit)
#' }
#'
#' @examples
#' # Build a minimal model component and fix a coefficient
#' dat <- data.frame(intercept = 1, x1 = rnorm(20), y = rnorm(20))
#' fm <- define_factor_model(n_factors = 1)
#' mc <- define_model_component(
#'   name = "y", data = dat, outcome = "y", factor = fm,
#'   covariates = c("intercept", "x1"), model_type = "linear",
#'   loading_normalization = 1
#' )
#'
#' # Fix intercept to 0 and x1 coefficient to 0.1
#' mc <- fix_coefficient(mc, covariate = "intercept", value = 0.0)
#' mc <- fix_coefficient(mc, covariate = "x1", value = 0.1)
#' length(mc$fixed_coefficients)  # 2 constraints stored on the component
#'
#' @export
fix_coefficient <- function(component, covariate, value, choice = NULL) {


  # ---- 1. Validate component ----
  # NULL pass-through: define_model_component() returns NULL with a warning
  # when its evaluation_indicator has no TRUE/1 rows. Pipelines that chain
  # fix_coefficient() onto every component should propagate the skip
  # rather than crash.
  if (is.null(component)) return(invisible(NULL))
  if (!inherits(component, "model_component")) {
    stop("`component` must be an object of class 'model_component'.")
  }

  # ---- 2. Validate covariate ----
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.")
  }

  if (!(covariate %in% component$covariates)) {
    stop("Covariate '", covariate, "' not found in component. ",
         "Available covariates: ", paste(component$covariates, collapse = ", "))
  }

  # ---- 3. Validate value ----
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("`value` must be a single finite numeric value.")
  }

  # ---- 4. Validate choice (for multinomial logit) ----
  if (component$model_type == "logit" && component$num_choices > 2) {
    # Multinomial logit: choice is required or defaults to applying to all choices
    if (!is.null(choice)) {
      choice <- as.integer(choice)
      if (!is.finite(choice) || choice < 1L || choice > (component$num_choices - 1L)) {
        stop("`choice` must be between 1 and ", component$num_choices - 1L,
             " (choice 0 is the reference category with no parameters).")
      }
    }
  } else {
    # Non-multinomial models: choice must be NULL
    if (!is.null(choice)) {
      stop("`choice` argument is only valid for multinomial logit models with num_choices > 2.")
    }
  }

  # ---- 5. Check for duplicate constraints ----
  for (fc in component$fixed_coefficients) {
    if (fc$covariate == covariate && identical(fc$choice, choice)) {
      stop("Coefficient for covariate '", covariate, "'",
           if (!is.null(choice)) paste0(" (choice ", choice, ")") else "",
           " is already fixed. Remove and recreate component to change.")
    }
  }

  # ---- 6. Add the constraint ----
  new_constraint <- list(
    covariate = covariate,
    value = value,
    choice = choice
  )

  component$fixed_coefficients <- c(component$fixed_coefficients, list(new_constraint))

  # ---- 7. Update parameter count ----
  # Reduce nparam_model by 1 for each fixed coefficient
  component$nparam_model <- component$nparam_model - 1L

  return(component)
}


#' Get the number of fixed coefficients in a model component
#'
#' @param component A model_component object
#' @return Integer count of fixed coefficients
#' @keywords internal
n_fixed_coefficients <- function(component) {
  length(component$fixed_coefficients)
}


#' Check if a covariate's coefficient is fixed
#'
#' @param component A model_component object
#' @param covariate Character. Name of the covariate to check.
#' @param choice Integer (optional). For multinomial logit, which choice to check.
#' @return Logical TRUE if fixed, FALSE otherwise
#' @keywords internal
is_coefficient_fixed <- function(component, covariate, choice = NULL) {
  for (fc in component$fixed_coefficients) {
    if (fc$covariate == covariate && identical(fc$choice, choice)) {
      return(TRUE)
    }
  }
  return(FALSE)
}


#' Get the fixed value for a covariate's coefficient
#'
#' @param component A model_component object
#' @param covariate Character. Name of the covariate.
#' @param choice Integer (optional). For multinomial logit, which choice.
#' @return The fixed value, or NULL if not fixed
#' @keywords internal
get_fixed_coefficient_value <- function(component, covariate, choice = NULL) {
  for (fc in component$fixed_coefficients) {
    if (fc$covariate == covariate && identical(fc$choice, choice)) {
      return(fc$value)
    }
  }
  return(NULL)
}


#' Fix type-specific intercepts to zero for a model component
#'
#' For models with n_types > 1, each component has type-specific intercepts
#' that shift the linear predictor for each non-reference type. This function
#' constrains those intercepts to zero, effectively removing the type-specific
#' shift for this component.
#'
#' @param component A model_component object from define_model_component()
#' @param types Integer vector (optional). Which types to fix. Must be between 2 and n_types
#'   (type 1 is the reference with no intercept). Default is NULL, meaning all
#'   non-reference types.
#' @param choice Integer (optional). For multinomial logit models with num_choices > 2,
#'   specifies which choice's type intercepts to fix. Must be between 1 and (num_choices - 1).
#'   Default is NULL, meaning all choices.
#'
#' @return The modified model_component object with the fixed type intercept constraints added.
#'
#' @details
#' When n_types > 1, the model includes a latent type structure where different
#' types can have different intercepts in each measurement equation. By default,
#' these intercepts are freely estimated. Fixing them to zero constrains the
#' outcome equation to have no type-specific shifts (though the type model loadings
#' at the factor level may still differ).
#'
#' This is useful when you want the type model to affect outcomes only through
#' factor loadings, not through direct type-specific intercepts.
#'
#' Fixed type intercepts are stored in the component and used during:
#' \itemize{
#'   \item Parameter initialization: fixed values are used directly (no estimation)
#'   \item Optimization: fixed parameters are excluded from the optimization
#'   \item Likelihood evaluation: fixed values are inserted into the parameter vector
#' }
#'
#' @examples
#' # Build a component with n_types = 2 and fix the type-2 intercept
#' dat <- data.frame(intercept = 1, y = rnorm(20))
#' fm <- define_factor_model(n_factors = 1, n_types = 2)
#' mc <- define_model_component(
#'   name = "y", data = dat, outcome = "y", factor = fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = 1, use_types = TRUE
#' )
#' mc <- fix_type_intercepts(mc)
#' length(mc$fixed_type_intercepts)
#'
#' @export
fix_type_intercepts <- function(component, types = NULL, choice = NULL) {

  # ---- 1. Validate component ----
  if (!inherits(component, "model_component")) {
    stop("`component` must be an object of class 'model_component'.")
  }

  # ---- 2. Get n_types from component's factor model ----
  n_types <- component$factor$n_types
  if (is.null(n_types) || n_types < 2L) {
    stop("Component's factor model has n_types < 2. Type intercepts only exist when n_types >= 2.")
  }

  # ---- 3. Validate types argument ----
  if (is.null(types)) {
    # Default: fix all non-reference types
    types <- 2L:n_types
  } else {
    types <- as.integer(types)
    if (any(!is.finite(types)) || any(types < 2L) || any(types > n_types)) {
      stop("`types` must be integers between 2 and ", n_types, " (type 1 is the reference).")
    }
    types <- unique(types)
  }

  # ---- 4. Validate choice (for multinomial logit) ----
  if (component$model_type == "logit" && component$num_choices > 2) {
    # Multinomial logit: determine which choices to fix
    if (!is.null(choice)) {
      choice <- as.integer(choice)
      if (!is.finite(choice) || choice < 1L || choice > (component$num_choices - 1L)) {
        stop("`choice` must be between 1 and ", component$num_choices - 1L,
             " (choice 0 is the reference category with no parameters).")
      }
      choices_to_fix <- choice
    } else {
      # Default: fix all choices
      choices_to_fix <- 1L:(component$num_choices - 1L)
    }
  } else {
    # Non-multinomial models: choice must be NULL
    if (!is.null(choice)) {
      stop("`choice` argument is only valid for multinomial logit models with num_choices > 2.")
    }
    choices_to_fix <- NULL
  }

  # ---- 5. Initialize fixed_type_intercepts list if needed ----
  if (is.null(component$fixed_type_intercepts)) {
    component$fixed_type_intercepts <- list()
  }

  # ---- 6. Add constraints and count new ones ----
  n_new_constraints <- 0L

  if (!is.null(choices_to_fix)) {
    # Multinomial logit: add constraint for each choice × type combination
    for (ch in choices_to_fix) {
      for (typ in types) {
        # Check for duplicate constraints
        already_fixed <- FALSE
        for (fti in component$fixed_type_intercepts) {
          if (fti$type == typ && identical(fti$choice, ch)) {
            already_fixed <- TRUE
            break
          }
        }
        if (!already_fixed) {
          new_constraint <- list(
            type = typ,
            value = 0.0,
            choice = ch
          )
          component$fixed_type_intercepts <- c(component$fixed_type_intercepts, list(new_constraint))
          n_new_constraints <- n_new_constraints + 1L
        }
      }
    }
  } else {
    # Non-multinomial: add constraint for each type
    for (typ in types) {
      # Check for duplicate constraints
      already_fixed <- FALSE
      for (fti in component$fixed_type_intercepts) {
        if (fti$type == typ && is.null(fti$choice)) {
          already_fixed <- TRUE
          break
        }
      }
      if (!already_fixed) {
        new_constraint <- list(
          type = typ,
          value = 0.0,
          choice = NULL
        )
        component$fixed_type_intercepts <- c(component$fixed_type_intercepts, list(new_constraint))
        n_new_constraints <- n_new_constraints + 1L
      }
    }
  }

  # Note: We do NOT reduce nparam_model here because the C++ layer

  # always expects type intercept parameters in the vector when n_types > 1.
  # The fixed type intercepts will be handled via parameter constraints in
  # optimize_model.R's setup_parameter_constraints() function.

  return(component)
}


#' Check if a type intercept is fixed
#'
#' @param component A model_component object
#' @param type Integer. Type number to check (2 to n_types).
#' @param choice Integer (optional). For multinomial logit, which choice to check.
#' @return Logical TRUE if fixed, FALSE otherwise
#' @keywords internal
is_type_intercept_fixed <- function(component, type, choice = NULL) {
  if (is.null(component$fixed_type_intercepts)) {
    return(FALSE)
  }
  for (fti in component$fixed_type_intercepts) {
    if (fti$type == type && identical(fti$choice, choice)) {
      return(TRUE)
    }
  }
  return(FALSE)
}


#' Get the number of fixed type intercepts in a model component
#'
#' @param component A model_component object
#' @return Integer count of fixed type intercepts
#' @keywords internal
n_fixed_type_intercepts <- function(component) {
  if (is.null(component$fixed_type_intercepts)) {
    return(0L)
  }
  length(component$fixed_type_intercepts)
}
