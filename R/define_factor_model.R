#' Define latent factor model structure
#'
#' Creates an object of class `"factor_model"` that specifies the structure of the
#' unobserved latent factors. This includes the number of factors and mixture components.
#' Loading constraints are specified at the component level via `define_model_component()`.
#' Numerical integration settings (quadrature points) are specified in `define_estimation_control()`.
#'
#' @param n_factors Integer. Number of latent factors (>=0). Use 0 for models without latent factors.
#' @param n_types Integer. Number of types (>=1)
#' @param factor_structure Character. Structure of factor dependencies. Options:
#'   - `"independent"` (default): Factors are independent
#'   - `"correlation"`: Correlated factors via Cholesky decomposition (2 factors only)
#'   - `"SE_linear"`: Structural equation f_k = alpha + alpha_1*f_1 + ... + epsilon
#'   - `"SE_quadratic"`: Adds quadratic terms: f_k = alpha + alpha_1*f_1 + alpha_q1*f_1^2 + ... + epsilon
#' @param n_mixtures Integer. Number of discrete mixtures (default = 1, allowed: 1-3)
#' @param factor_covariates Character vector. Names of covariates that shift factor means.
#'   When specified, the factor distribution becomes f_i ~ N((X_i - mean(X)) * gamma, Sigma)
#'   where X_i is the covariate vector for observation i and gamma is a matrix of coefficients
#'   (one column per factor). This allows factor means to vary by observed characteristics.
#'   The covariates must be present in the data passed to \code{estimate_model_rcpp()}.
#'   **Note:** All covariates are automatically demeaned internally to ensure the overall
#'   factor mean is zero (required for identification). Including an intercept/constant
#'   has no effect since it becomes zero after demeaning.
#'   For SE_linear/SE_quadratic structures, covariates only affect input factor means
#'   (not the outcome factor which is determined by the structural equation).
#' @param se_covariates Character vector. Names of covariates that directly affect the
#'   outcome factor in SE_linear/SE_quadratic models. When specified, the structural
#'   equation becomes f_k = intercept + sum(alpha_j * f_j) + sum(beta_m * X_m) + epsilon.
#'   Only valid for factor_structure = "SE_linear" or "SE_quadratic".
#'   All covariates are automatically demeaned internally for identification.
#'
#' @return An object of class "factor_model"
#' @examples
#' # Single factor model
#' fm <- define_factor_model(n_factors = 1)
#'
#' # Two-factor structural equation model
#' fm_se <- define_factor_model(n_factors = 2, factor_structure = "SE_linear")
#' @export
define_factor_model <- function(n_factors,
                                n_types = 1,
                                factor_structure = "independent",
                                n_mixtures = 1,
                                factor_covariates = NULL,
                                se_covariates = NULL) {

  # ---- 1. Input validation ----
  # Check all arguments are the correct type and within supported range.

  if (!is.numeric(n_factors) || n_factors < 0) stop("n_factors must be a non-negative integer.")
  if (!is.numeric(n_types) || n_types < 1) stop("n_types must be a positive integer")
  if (!n_mixtures %in% 1:3) stop("n_mixtures should be between 1-3") #currently this is the case, might change later

  # Validate factor_structure
  valid_structures <- c("independent", "correlation", "SE_linear", "SE_quadratic")
  if (!factor_structure %in% valid_structures) {
    stop("factor_structure must be one of: ", paste(valid_structures, collapse = ", "))
  }

  # Validate factor_structure requirements
  if (factor_structure == "correlation" && n_factors > 2) {
    stop("Correlated factor models are currently only supported for n_factors = 2. ",
         "Models with 3+ correlated factors are not yet implemented.")
  }
  if (factor_structure %in% c("SE_linear", "SE_quadratic") && n_factors < 2) {
    stop("Structural equation models require at least 2 factors (1 input + 1 outcome).")
  }

  # Validate factor_covariates
  if (!is.null(factor_covariates)) {
    if (!is.character(factor_covariates)) {
      stop("factor_covariates must be a character vector of covariate names")
    }
    if (length(factor_covariates) == 0) {
      factor_covariates <- NULL
    } else {
      # Note: All covariates are automatically demeaned internally to ensure
      # E[f_k] = 0 for identification (factor mean intercept is not separately
      # identified from measurement equation intercepts).
      #
      # If user includes 'intercept' or 'constant', we REMOVE it because after
      # demeaning it becomes a column of zeros, making its coefficient unidentified
      # (gradient is always zero, Hessian is singular for that parameter).
      intercept_names <- c("intercept", "constant", "Intercept", "Constant", "(Intercept)", "1")
      intercept_in_covs <- factor_covariates %in% intercept_names
      if (any(intercept_in_covs)) {
        removed <- factor_covariates[intercept_in_covs]
        factor_covariates <- factor_covariates[!intercept_in_covs]
        message("Removing intercept/constant from factor_covariates: ",
                paste(removed, collapse = ", "), ". ",
                "All covariates are automatically demeaned for identification, ",
                "so a constant term becomes zero (unidentified coefficient).")
        if (length(factor_covariates) == 0) {
          factor_covariates <- NULL
        }
      }
    }
  }

  # Validate se_covariates
  if (!is.null(se_covariates)) {
    if (!factor_structure %in% c("SE_linear", "SE_quadratic")) {
      stop("se_covariates is only valid for factor_structure = 'SE_linear' or 'SE_quadratic'")
    }
    if (!is.character(se_covariates)) {
      stop("se_covariates must be a character vector of covariate names")
    }
    if (length(se_covariates) == 0) {
      se_covariates <- NULL
    } else {
      # Remove intercept/constant terms (become zero after demeaning)
      intercept_names <- c("intercept", "constant", "Intercept", "Constant", "(Intercept)", "1")
      intercept_in_covs <- se_covariates %in% intercept_names
      if (any(intercept_in_covs)) {
        removed <- se_covariates[intercept_in_covs]
        se_covariates <- se_covariates[!intercept_in_covs]
        message("Removing intercept/constant from se_covariates: ",
                paste(removed, collapse = ", "), ". ",
                "SE covariates are automatically demeaned, ",
                "so a constant term becomes zero (unidentified coefficient).")
        if (length(se_covariates) == 0) {
          se_covariates <- NULL
        }
      }
    }
  }

  # Compute number of factor mean parameters
  # For each covariate, we estimate one coefficient per factor
  # For SE models, only input factors have mean parameters (outcome factor mean
  # is determined by the structural equation)
  n_factor_cov <- if (!is.null(factor_covariates)) length(factor_covariates) else 0L
  if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
    # Only input factors get covariate effects on mean
    n_factors_with_mean <- n_factors - 1L
  } else {
    n_factors_with_mean <- n_factors
  }
  n_factor_mean_param <- as.integer(n_factor_cov * n_factors_with_mean)

  # ---- 2. Compute number of variance/covariance parameters ----
  # Depends on factor_structure:
  # - "independent": k variances (one per factor)
  # - "correlation": k*(k+1)/2 (variances + covariances via Cholesky)
  # - "SE_linear", "SE_quadratic": (k-1) input factor variances only (outcome variance is derived)

  if (factor_structure == "correlation") {
    f_nvariance <- n_factors * (n_factors + 1L) / 2L
  } else if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
    # For SE models: only input factors (first k-1) have variance parameters
    # The last factor's distribution is determined by the structural equation
    f_nvariance <- n_factors - 1L
  } else {
    f_nvariance <- n_factors
  }

  # ---- 3. Compute SE parameters ----
  # For SE_linear with k factors (k-1 input, 1 outcome):
  # - 1 intercept (alpha)
  # - (k-1) linear coefficients (alpha_1, ..., alpha_{k-1})
  # - 1 residual variance (sigma^2_epsilon)
  # Total: 1 + (k-1) + 1 = k + 1 parameters, e.g., 3 for 2-factor model
  #
  # For SE_quadratic: adds (k-1) quadratic coefficients (for f_1^2, f_2^2, ...)
  # Total: 1 + (k-1) + (k-1) + 1 = 2k parameters, e.g., 4 for 2-factor model

  if (factor_structure == "SE_linear") {
    n_input_factors <- n_factors - 1L
    nse_param <- 1L + n_input_factors + 1L  # intercept + linear coefs + residual var
  } else if (factor_structure == "SE_quadratic") {
    n_input_factors <- n_factors - 1L
    nse_param <- 1L + n_input_factors + n_input_factors + 1L  # intercept + linear + quadratic + residual var
  } else {
    nse_param <- 0L
  }

  # Compute number of SE covariate parameters (one coefficient per covariate)
  n_se_cov <- if (!is.null(se_covariates)) length(se_covariates) else 0L
  n_se_covariate_param <- as.integer(n_se_cov)

  # ---- 4. Compute total number of parameters for factor distribution ----
  # Combines variances/covariances with mixture-related parameters and SE parameters.
  #
  # Parameter structure for mixtures (nmix > 1):
  # - f_nvariance * nmix variances (one per factor per mixture)
  # - (nmix - 1) * n_factors_for_mixture means (for non-reference mixtures)
  # - (nmix - 1) log-weights (for non-reference mixtures)
  #
  # For SE models, mixture means only apply to input factors (n_factors - 1)
  # since the outcome factor is determined by the structural equation.
  n_factors_for_mixture <- if (factor_structure %in% c("SE_linear", "SE_quadratic")) {
    n_factors - 1L
  } else {
    n_factors
  }

  nfac_param <- as.integer(f_nvariance * n_mixtures +
                             (n_mixtures - 1L) * n_factors_for_mixture +
                             (n_mixtures - 1L) +
                             nse_param +
                             n_factor_mean_param +
                             n_se_covariate_param)

  # ---- 5. Compute type model parameter count ----
  # For n_types > 1, the type probability model has:
  # - (n_types - 1) * n_factors factor loadings (type 1 is reference)
  # Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
  # Note: Type-specific intercepts are added per model component, not here.
  ntyp_param <- if (n_types > 1L) as.integer((n_types - 1L) * n_factors) else 0L

  # ---- 6. Construct the factor_model object ----
  # Bundle all parameters and metadata into a list, and assign class.

  out <- list(
    n_factors = as.integer(n_factors),
    n_types = as.integer(n_types),
    factor_structure = factor_structure,
    correlation = (factor_structure == "correlation"),  # For backward compatibility
    n_mixtures = as.integer(n_mixtures),
    nfac_param = nfac_param,
    nse_param = as.integer(nse_param),
    ntyp_param = ntyp_param,
    factor_covariates = factor_covariates,
    n_factor_covariates = n_factor_cov,
    n_factor_mean_param = n_factor_mean_param,
    se_covariates = se_covariates,
    n_se_covariates = n_se_cov,
    n_se_covariate_param = n_se_covariate_param,
    params = rep(0.0, nfac_param)
  )

  class(out) <- "factor_model"
  return(out)
}


#' @export
print.factor_model <- function(x, ...) {
  cat("Factor Model\n")
  cat("------------\n")
  cat("Number of latent factors:", x$n_factors, "\n")
  cat("Number of types:", x$n_types, "\n")
  cat("Factor structure:", x$factor_structure, "\n")
  if (x$factor_structure == "SE_linear") {
    cat("  Structural equation: f_k = alpha + alpha_1*f_1 + ... + epsilon\n")
    cat("  SE parameters:", x$nse_param, "(intercept, linear coefs, residual var)\n")
  } else if (x$factor_structure == "SE_quadratic") {
    cat("  Structural equation: f_k = alpha + alpha_1*f_1 + alpha_q1*f_1^2 + ... + epsilon\n")
    cat("  SE parameters:", x$nse_param, "(intercept, linear coefs, quadratic coefs, residual var)\n")
  }
  if (!is.null(x$factor_covariates) && length(x$factor_covariates) > 0) {
    cat("Factor mean covariates:", paste(x$factor_covariates, collapse = ", "), "\n")
    cat("  Factor mean parameters:", x$n_factor_mean_param, "\n")
    cat("  Factor distribution: f_i ~ N(X_i * gamma, Sigma)\n")
  }
  cat("Number of mixtures:", x$n_mixtures, "\n")
  cat("Number of parameters in latent factor distribution:", x$nfac_param, "\n")
  if (x$n_types > 1L) {
    cat("Number of parameters in type probability model:", x$ntyp_param, "\n")
  }
  invisible(x)
}

