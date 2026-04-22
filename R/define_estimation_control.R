#' Define estimation control settings
#'
#' @param n_quad_points Integer. Number of Gauss-Hermite quadrature points for numerical integration (default = 16)
#' @param num_cores Integer. Number of processes to use for parallel estimation (default = 1)
#' @param cluster_type Character. Type of parallel cluster to use: "auto" (default), "FORK", or "PSOCK".
#'   "auto" uses FORK on Unix (faster, shared memory) and PSOCK on Windows.
#'   "FORK" forces fork-based parallelism (Unix only, faster due to shared memory).
#'   "PSOCK" forces socket-based parallelism (works on all platforms, more overhead).
#' @param adaptive_integration Logical. Whether to use adaptive integration in second-stage estimation
#'   (default = FALSE). When TRUE, the number of quadrature points per observation is determined
#'   based on the standard error of factor scores from a previous estimation stage.
#' @param adapt_int_thresh Numeric. Threshold for adaptive integration (default = 0.5).
#'   Smaller values use more integration points. The formula is:
#'   \code{n_quad_obs = 1 + 2 * floor(factor_se / factor_var / adapt_int_thresh)}.
#'   When factor_se is small relative to factor variance, fewer quadrature points are used.
#'   Legacy code default is 0.5.
#'
#' @details
#' Adaptive integration is useful in two-stage estimation where factor scores have been
#' estimated in a first stage. The key insight is that observations with well-identified
#' factor scores (small SE) need fewer integration points, while observations with
#' poorly-identified factors (large SE) need full quadrature.
#'
#' To use adaptive integration:
#' \enumerate{
#'   \item Estimate Stage 1 model to get factor scores via \code{estimate_factorscores_rcpp()}
#'   \item Initialize FactorModel for Stage 2 via \code{initialize_factor_model_cpp()}
#'   \item Call \code{set_adaptive_quadrature_cpp()} with factor scores, SEs, and variances
#'   \item Run estimation - the adaptive settings are applied automatically
#' }
#'
#' The diagnostic output from \code{set_adaptive_quadrature_cpp(verbose=TRUE)} shows:
#' \itemize{
#'   \item Distribution of integration points across observations
#'   \item Average integration points vs standard (shows computational savings)
#'   \item Computational reduction percentage
#' }
#'
#' @return An object of class estimation_control containing control settings
#' @export
define_estimation_control <- function(n_quad_points = 16, num_cores = 1,
                                      cluster_type = "auto",
                                      adaptive_integration = FALSE,
                                      adapt_int_thresh = 0.5) {
  if (!is.numeric(n_quad_points) || n_quad_points < 1) {
    stop("n_quad_points must be a positive integer.")
  }
  if (!is.numeric(num_cores) || num_cores < 1) {
    stop("num_cores must be a positive integer.")
  }
  cluster_type <- match.arg(cluster_type, c("auto", "FORK", "PSOCK"))
  if (cluster_type == "FORK" && .Platform$OS.type != "unix") {
    warning("FORK clusters are only available on Unix. Using PSOCK instead.")
    cluster_type <- "PSOCK"
  }
  if (!is.logical(adaptive_integration)) {
    stop("adaptive_integration must be TRUE or FALSE.")
  }
  if (!is.numeric(adapt_int_thresh) || adapt_int_thresh <= 0) {
    stop("adapt_int_thresh must be a positive number.")
  }

  out <- list(
    n_quad_points = as.integer(n_quad_points),
    num_cores = as.integer(num_cores),
    cluster_type = cluster_type,
    adaptive_integration = adaptive_integration,
    adapt_int_thresh = adapt_int_thresh
  )

  class(out) <- "estimation_control"
  return(out)
}


# TODO: should be a function that calculates and outputs how many parameters
# placeholder below
#
# component_param_count <- function(mc) {
#   # If you already computed/stored it, trust that:
#   if (!is.null(mc$nparam_model)) return(as.integer(mc$nparam_model))
#
#   # Fallback: derive something sensible (adjust if your spec changes)
#   k_cov  <- if (!is.null(mc$covariates)) length(mc$covariates) else 0L
#   k_int  <- if (isTRUE(mc$intercept)) 1L else 0L
#   k_fac  <- if (!is.null(mc$factor$n_factors)) as.integer(mc$factor$n_factors) else 0L
#   k_types_minus_1 <- if (!is.null(mc$factor$n_types)) max(as.integer(mc$factor$n_types) - 1L, 0L) else 0L
#   as.integer(k_cov + k_int + k_fac + k_types_minus_1)
# }
#
# param_counts <- function(model_system, factor_model) {
#   stopifnot(inherits(model_system, "model_system"),
#             inherits(factor_model, "factor_model"))
#
#   comps <- model_system$components
#   per_component <- vapply(comps, component_param_count, integer(1))
#   names(per_component) <- names(comps)
#
#   total_component_params <- sum(per_component)
#   factor_dist_params     <- as.integer(factor_model$nfac_param)
#   total_params           <- as.integer(total_component_params + factor_dist_params)
#
#   list(
#     per_component            = per_component,                 # named int vector
#     total_component_params   = total_component_params,        # int
#     factor_dist_params       = factor_dist_params,            # int
#     total_params             = total_params                   # int
#   )
# }





#' @export
print.estimation_control <- function(x, ...) {
  cat("Estimation Control\n")
  cat("------------------\n")
  cat("Number of quadrature points:", x$n_quad_points, "\n")
  cat("Number of cores for parallelization:", x$num_cores, "\n")
  cat("Cluster type:", x$cluster_type, "\n")
  if (isTRUE(x$adaptive_integration)) {
    cat("Adaptive integration: enabled (threshold =", x$adapt_int_thresh, ")\n")
  } else {
    cat("Adaptive integration: disabled\n")
  }
}
