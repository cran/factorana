#' Estimate Factor Scores
#'
#' Estimates factor scores for each observation after model estimation.
#' The model parameters are held fixed at their estimated values, and
#' factor scores are computed as the posterior mode for each observation.
#'
#' @param result A factorana_result object from estimate_model_rcpp()
#' @param data Data frame containing all variables (same as used in estimation)
#' @param control Optional estimation control object. If NULL, uses default.
#' @param parallel Whether to use parallel computation (default FALSE).
#'   When TRUE, uses the num_cores setting from control.
#' @param verbose Whether to print progress (default TRUE)
#' @param include_prior Whether to include the factor prior in likelihood/SE
#'   computation (default FALSE). When FALSE, matches legacy C++ behavior where
#'   SEs are based only on observation likelihood. When TRUE, includes the prior
#'   contribution to the Hessian, resulting in smaller SEs.
#' @param id_var Optional character string specifying the name of an ID variable
#'   in data. If provided, this variable is included in the output data frame for
#'   easier merging with other datasets. The ID values are taken from the original
#'   data in the order observations were processed. Note: The ID column must be
#'   numeric (not character) since data is converted to a numeric matrix. For
#'   character IDs, use the obs_id column to merge with your original data.
#'
#' @details
#' Factor scores are estimated by maximizing the posterior density:
#' \deqn{L(f|y_i, \theta) = p(y_i|f, \theta) \cdot \phi(f|0, \sigma^2)}
#'
#' where:
#' \itemize{
#'   \item \eqn{f} is the vector of factor values for observation i
#'   \item \eqn{y_i} is the observed data for observation i
#'   \item \eqn{\theta} are the fixed model parameters (from previous estimation)
#'   \item \eqn{\phi(f|0, \sigma^2)} is the normal prior on factors
#' }
#'
#' Standard errors are computed from the diagonal of the inverse Hessian
#' of the log-posterior at the mode. By default (include_prior=FALSE), the SE
#' is based only on the observation likelihood Hessian, matching the legacy
#' C++ implementation. Set include_prior=TRUE to include the prior's Hessian
#' contribution (-1/sigma^2) which produces smaller SEs.
#'
#' When parallel=TRUE, observations are distributed across cores using
#' doParallel/foreach, with each worker processing a subset of observations.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{obs_id} - Observation index (1-based)
#'   \item \code{<id_var>} - ID variable values (if id_var was specified)
#'   \item \code{factor_1, factor_2, ...} - Estimated factor scores
#'   \item \code{se_factor_1, se_factor_2, ...} - Standard errors
#'   \item \code{converged} - Whether optimization converged for this observation
#'   \item \code{log_posterior} - Log-posterior value at the mode
#' }
#'
#' @examples
#' \donttest{
#' # Estimate a small one-factor model, then recover factor scores
#' set.seed(1); n <- 200
#' f <- rnorm(n)
#' dat <- data.frame(intercept = 1,
#'                   y1 = 1.0 * f + rnorm(n, 0, 0.5),
#'                   y2 = 0.8 * f + rnorm(n, 0, 0.5))
#' fm <- define_factor_model(n_factors = 1)
#' mc1 <- define_model_component("m1", dat, "y1", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = 1)
#' mc2 <- define_model_component("m2", dat, "y2", fm,
#'   covariates = "intercept", model_type = "linear",
#'   loading_normalization = NA_real_)
#' ms <- define_model_system(components = list(mc1, mc2), factor = fm)
#' ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
#' fit <- estimate_model_rcpp(ms, dat, control = ctrl,
#'   optimizer = "nlminb", parallel = FALSE, verbose = FALSE)
#'
#' fscores <- estimate_factorscores_rcpp(fit, dat, control = ctrl,
#'                                        parallel = FALSE, verbose = FALSE)
#' head(fscores)
#' }
#'
#' @export
estimate_factorscores_rcpp <- function(result, data, control = NULL,
                                        parallel = FALSE, verbose = TRUE,
                                        include_prior = FALSE,
                                        id_var = NULL) {

  # Validate input
  if (!inherits(result, "factorana_result")) {
    stop("result must be a factorana_result object from estimate_model_rcpp()")
  }

  # Validate id_var if provided
  id_values <- NULL
  if (!is.null(id_var)) {
    if (!is.character(id_var) || length(id_var) != 1) {
      stop("id_var must be a single character string")
    }
    if (!id_var %in% names(data)) {
      stop(sprintf("id_var '%s' not found in data", id_var))
    }
    id_values <- data[[id_var]]
  }

  # Get model system and estimates from result
  model_system <- result$model_system
  estimates <- result$estimates

  # Deep copy model_system to avoid C++ reuse issues
  model_system <- unserialize(serialize(model_system, NULL))

  # Get control from result if not provided
  if (is.null(control)) {
    control <- define_estimation_control()
  }

  n_factors <- model_system$factor$n_factors
  if (n_factors == 0) {
    stop("No factors to estimate (n_factors = 0)")
  }

  # Convert data to matrix (requires all columns to be numeric)
  data_mat <- as.matrix(data)
  n_obs <- nrow(data_mat)

  if (verbose) {
    mode_str <- if (parallel) sprintf("parallel (%d cores)", control$num_cores) else "serial"
    message(sprintf("Estimating factor scores for %d observations with %d factor(s) [%s]...",
                    n_obs, n_factors, mode_str))
  }

  # Use parallel or serial implementation
  if (parallel && control$num_cores > 1) {
    result_df <- estimate_factorscores_parallel(
      model_system, data_mat, estimates, n_factors, control, verbose, include_prior
    )
  } else {
    result_df <- estimate_factorscores_serial(
      model_system, data_mat, estimates, n_factors, control, verbose, include_prior
    )
  }

  # Add ID variable if specified
  if (!is.null(id_var) && !is.null(id_values)) {
    # Insert ID column as the first column after obs_id
    result_df <- cbind(
      obs_id = result_df$obs_id,
      setNames(data.frame(id_values), id_var),
      result_df[, -1, drop = FALSE]
    )
  }

  return(result_df)
}


#' Serial factor score estimation (internal)
#' @noRd
estimate_factorscores_serial <- function(model_system, data_mat, estimates,
                                          n_factors, control, verbose,
                                          include_prior = FALSE) {
  n_obs <- nrow(data_mat)
  n_quad <- control$n_quad_points

  # Initialize the C++ FactorModel
  fm_ptr <- initialize_factor_model_cpp(model_system, data_mat, n_quad)

  # Get parameter info
  param_info <- get_parameter_info_cpp(fm_ptr)

  if (length(estimates) != param_info$n_param) {
    stop(sprintf("Parameter count mismatch: result has %d, C++ expects %d",
                 length(estimates), param_info$n_param))
  }

  # Initialize results storage
  factor_scores <- matrix(NA_real_, nrow = n_obs, ncol = n_factors)
  factor_ses <- matrix(NA_real_, nrow = n_obs, ncol = n_factors)
  converged <- logical(n_obs)
  log_posteriors <- numeric(n_obs)

  # Progress tracking
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_obs, style = 3)
  }

  # Loop over observations
  for (iobs in seq_len(n_obs)) {
    obs_result <- estimate_single_factorscore(
      fm_ptr, iobs, estimates, n_factors, include_prior
    )

    factor_scores[iobs, ] <- obs_result$factors
    factor_ses[iobs, ] <- obs_result$ses
    converged[iobs] <- obs_result$converged
    log_posteriors[iobs] <- obs_result$log_posterior

    if (verbose) {
      setTxtProgressBar(pb, iobs)
    }
  }

  if (verbose) {
    close(pb)
    n_converged <- sum(converged)
    message(sprintf("Factor score estimation complete: %d/%d converged (%.1f%%)",
                    n_converged, n_obs, 100 * n_converged / n_obs))
  }

  # Create and return result data frame
  create_factorscore_dataframe(factor_scores, factor_ses, converged,
                                log_posteriors, n_factors)
}


#' Parallel factor score estimation (internal)
#' @noRd
estimate_factorscores_parallel <- function(model_system, data_mat, estimates,
                                            n_factors, control, verbose,
                                            include_prior = FALSE) {
  n_obs <- nrow(data_mat)
  n_cores <- control$num_cores

  # Set up parallel backend
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is required for parallel factor score estimation")
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' is required for parallel factor score estimation")
  }

  # Create cluster
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # Ensure cluster is stopped on exit
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)

  # Export library paths to workers (ensures workers can find factorana)
  current_lib_paths <- .libPaths()
  parallel::clusterExport(cl, "current_lib_paths", envir = environment())
  parallel::clusterEvalQ(cl, .libPaths(current_lib_paths))

  # Split observations across workers
  obs_splits <- split(seq_len(n_obs), cut(seq_len(n_obs), n_cores, labels = FALSE))

  if (verbose) {
    message(sprintf("  Distributing %d observations across %d workers...",
                    n_obs, n_cores))
  }

  # Serialize model_system for worker transfer
  # Note: model_system no longer contains data (stripped in define_model_component)
  model_system_raw <- serialize(model_system, NULL)
  n_quad <- control$n_quad_points

  # Run parallel computation
  `%dopar%` <- foreach::`%dopar%`

  worker_results <- foreach::foreach(
    obs_idx = obs_splits,
    .packages = c("factorana"),
    .export = c("estimate_single_factorscore", "include_prior")
  ) %dopar% {
    # Deserialize model_system in worker
    ms_local <- unserialize(model_system_raw)

    # Initialize C++ model for this worker's data
    # We still need all data because observations reference column indices
    fm_ptr <- initialize_factor_model_cpp(ms_local, data_mat, n_quad)

    # Process assigned observations
    n_local <- length(obs_idx)
    local_factors <- matrix(NA_real_, nrow = n_local, ncol = n_factors)
    local_ses <- matrix(NA_real_, nrow = n_local, ncol = n_factors)
    local_converged <- logical(n_local)
    local_log_post <- numeric(n_local)

    for (i in seq_along(obs_idx)) {
      iobs <- obs_idx[i]
      obs_result <- estimate_single_factorscore(
        fm_ptr, iobs, estimates, n_factors, include_prior
      )

      local_factors[i, ] <- obs_result$factors
      local_ses[i, ] <- obs_result$ses
      local_converged[i] <- obs_result$converged
      local_log_post[i] <- obs_result$log_posterior
    }

    list(
      obs_idx = obs_idx,
      factors = local_factors,
      ses = local_ses,
      converged = local_converged,
      log_posterior = local_log_post
    )
  }

  # Combine results from all workers
  factor_scores <- matrix(NA_real_, nrow = n_obs, ncol = n_factors)
  factor_ses <- matrix(NA_real_, nrow = n_obs, ncol = n_factors)
  converged <- logical(n_obs)
  log_posteriors <- numeric(n_obs)

  for (wr in worker_results) {
    idx <- wr$obs_idx
    factor_scores[idx, ] <- wr$factors
    factor_ses[idx, ] <- wr$ses
    converged[idx] <- wr$converged
    log_posteriors[idx] <- wr$log_posterior
  }

  if (verbose) {
    n_converged <- sum(converged)
    message(sprintf("Factor score estimation complete: %d/%d converged (%.1f%%)",
                    n_converged, n_obs, 100 * n_converged / n_obs))
  }

  # Create and return result data frame
  create_factorscore_dataframe(factor_scores, factor_ses, converged,
                                log_posteriors, n_factors)
}


#' Estimate factor scores for a single observation (internal)
#' @noRd
estimate_single_factorscore <- function(fm_ptr, iobs, estimates, n_factors,
                                         include_prior = FALSE) {
  # Define optimization functions
  # Note: For optimization, we always include the prior to find the posterior mode.
  # The include_prior flag only affects SE computation (observation-only vs posterior).
  objective <- function(fac_values) {
    res <- evaluate_factorscore_likelihood_cpp(
      fm_ptr, iobs - 1,  # Convert to 0-based index
      fac_values, estimates,
      compute_gradient = FALSE,
      compute_hessian = FALSE,
      include_prior = TRUE  # Always include prior for finding the mode
    )
    return(-res$logLikelihood)  # Negative for minimization
  }

  gradient <- function(fac_values) {
    res <- evaluate_factorscore_likelihood_cpp(
      fm_ptr, iobs - 1,
      fac_values, estimates,
      compute_gradient = TRUE,
      compute_hessian = FALSE,
      include_prior = TRUE  # Always include prior for finding the mode
    )
    return(-res$gradient)  # Negative for minimization
  }

  # Extract factor variances for clipping bounds
  # Clip at ±4.688738939305818364688 * sqrt(2) * sigma ≈ ±6.63 * sigma
  # This matches legacy C++ behavior (TMinLkhd.cc:4816)
  clip_multiplier <- 4.688738939305818364688 * sqrt(2)
  factor_sds <- numeric(n_factors)
  for (k in seq_len(n_factors)) {
    var_name <- paste0("factor_var_", k)
    if (var_name %in% names(estimates)) {
      factor_sds[k] <- sqrt(estimates[var_name])
    } else {
      factor_sds[k] <- 1.0  # Default if not found
    }
  }
  clip_bounds <- clip_multiplier * factor_sds

  # Use nlminb with multiple starting points for all models
  # nlminb uses the gradient for efficient convergence
  best_result <- NULL
  best_obj <- Inf
  best_grad_norm <- Inf

  # Starting points for optimization - include more points for robustness
  if (n_factors == 1) {
    # Multiple starting points for 1-factor models
    start_points <- list(0, -1, 1, -2, 2, -3, 3)
  } else if (n_factors == 2) {
    start_points <- list(c(0, 0), c(-1, 0), c(1, 0), c(0, -1), c(0, 1),
                         c(-2, 0), c(2, 0), c(0, -2), c(0, 2))
  } else {
    # For 3+ factors, use origin plus perturbations along each axis
    start_points <- list(rep(0, n_factors))
    for (k in seq_len(n_factors)) {
      for (val in c(1, -1, 2, -2)) {
        vec <- rep(0, n_factors)
        vec[k] <- val
        start_points <- c(start_points, list(vec))
      }
    }
  }

  # Gradient tolerance for accepting a solution
  # This catches pathological cases (gradient ~10^7) while allowing reasonable solutions
  grad_tol <- 1000  # Accept if max|gradient| < 1000

  for (init_fac in start_points) {
    init_fac <- as.numeric(init_fac)

    # Use nlminb with R defaults - this is a well-scaled, smooth problem
    # DO NOT CHANGE THESE CONVERGENCE SETTINGS without user permission
    opt_result <- tryCatch({
      stats::nlminb(
        start = init_fac,
        objective = objective,
        gradient = gradient
        # Using R defaults: iter.max=150, eval.max=200, rel.tol=1e-10
      )
    }, error = function(e) {
      NULL
    })

    if (!is.null(opt_result) && is.finite(opt_result$objective)) {
      # Check gradient magnitude at solution
      grad_at_solution <- gradient(opt_result$par)
      grad_norm <- max(abs(grad_at_solution))

      # Accept this solution if it's better and has reasonable gradient
      if (opt_result$objective < best_obj && grad_norm < best_grad_norm) {
        best_obj <- opt_result$objective
        best_result <- opt_result
        best_grad_norm <- grad_norm
      }

      # If we have a good solution (small gradient), we can stop
      if (grad_norm < grad_tol && opt_result$convergence == 0) {
        break
      }
    }
  }

  if (is.null(best_result) || !is.finite(best_obj)) {
    return(list(
      factors = rep(NA_real_, n_factors),
      ses = rep(NA_real_, n_factors),
      converged = FALSE,
      log_posterior = NA_real_
    ))
  }

  # Get the factor scores
  factor_scores <- best_result$par

  # Check if any data contributed for this observation. If all eval indicators
  # are 0, the observation-only log-likelihood at any factor value is 0 (all
  # model densities are 1). In that case the factor score is undefined — there
  # is simply no information. Return NA for both score and SE.
  obs_only <- evaluate_factorscore_likelihood_cpp(
    fm_ptr, iobs - 1, factor_scores, estimates,
    compute_gradient = FALSE, compute_hessian = FALSE,
    include_prior = FALSE
  )
  if (abs(obs_only$logLikelihood) < 1e-12) {
    return(list(
      factors = rep(NA_real_, n_factors),
      ses = rep(NA_real_, n_factors),
      converged = FALSE,
      log_posterior = NA_real_
    ))
  }

  # Clip extreme factor scores to ±6.63*sigma (matches legacy C++ behavior)
  # This prevents extreme values that are outside the Gauss-Hermite quadrature range
  clipped <- FALSE
  for (k in seq_len(n_factors)) {
    if (abs(factor_scores[k]) > clip_bounds[k]) {
      factor_scores[k] <- sign(factor_scores[k]) * clip_bounds[k]
      clipped <- TRUE
    }
  }

  # Consider converged only if gradient is small
  # A large gradient indicates the optimizer stopped prematurely
  converged <- best_grad_norm < grad_tol

  # Compute standard errors from Hessian at optimum (or clipped value)
  # include_prior controls whether SE is based on observation-only or posterior Hessian
  ses <- compute_factorscore_ses(fm_ptr, iobs, factor_scores, estimates, n_factors,
                                  include_prior)

  # When include_prior=TRUE, the posterior SE is mathematically bounded above by the
  # prior SD (data always adds information, never removes it, for log-concave
  # likelihoods). If numerical issues push an SE above the prior SD, clamp it.
  if (include_prior) {
    for (k in seq_len(n_factors)) {
      if (!is.na(ses[k]) && ses[k] > factor_sds[k]) {
        ses[k] <- factor_sds[k]
      }
    }
  }

  list(
    factors = factor_scores,
    ses = ses,
    converged = converged,
    log_posterior = -best_result$objective
  )
}


#' Compute standard errors from Hessian (internal)
#' @noRd
compute_factorscore_ses <- function(fm_ptr, iobs, fac_values, estimates, n_factors,
                                     include_prior = FALSE) {
  hess_result <- evaluate_factorscore_likelihood_cpp(
    fm_ptr, iobs - 1,
    fac_values, estimates,
    compute_gradient = FALSE,
    compute_hessian = TRUE,
    include_prior = include_prior  # Controls whether prior Hessian is included
  )

  # Convert upper triangle to full matrix
  hess_vec <- hess_result$hessian
  hess_mat <- matrix(0, n_factors, n_factors)
  idx <- 1
  for (i in 1:n_factors) {
    for (j in i:n_factors) {
      hess_mat[i, j] <- hess_vec[idx]
      hess_mat[j, i] <- hess_vec[idx]
      idx <- idx + 1
    }
  }

  # Invert Hessian to get covariance matrix
  tryCatch({
    inv_hess <- solve(-hess_mat)
    sqrt(pmax(0, diag(inv_hess)))
  }, error = function(e) {
    # If inversion fails, try pseudoinverse
    tryCatch({
      svd_result <- svd(-hess_mat)
      tol <- max(dim(hess_mat)) * max(svd_result$d) * .Machine$double.eps * 100
      pos_idx <- svd_result$d > tol
      d_inv <- rep(0, length(svd_result$d))
      d_inv[pos_idx] <- 1 / svd_result$d[pos_idx]
      inv_hess <- svd_result$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(svd_result$u)
      sqrt(pmax(0, diag(inv_hess)))
    }, error = function(e2) {
      rep(NA_real_, n_factors)
    })
  })
}


#' Create factor score result data frame (internal)
#' @noRd
create_factorscore_dataframe <- function(factor_scores, factor_ses, converged,
                                          log_posteriors, n_factors) {
  n_obs <- nrow(factor_scores)
  result_df <- data.frame(obs_id = seq_len(n_obs))

  # Add factor scores
  for (k in seq_len(n_factors)) {
    result_df[[sprintf("factor_%d", k)]] <- factor_scores[, k]
  }

  # Add standard errors
  for (k in seq_len(n_factors)) {
    result_df[[sprintf("se_factor_%d", k)]] <- factor_ses[, k]
  }

  # Add convergence and log-posterior
  result_df$converged <- converged
  result_df$log_posterior <- log_posteriors

  # Add class
  class(result_df) <- c("factorana_factorscores", "data.frame")

  result_df
}


#' Print method for factorana_factorscores
#'
#' @param x A factorana_factorscores object
#' @param n Number of rows to print (default 10)
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns \code{x}. Called for its side effect of printing
#'   a summary of the factor score estimates (convergence rate, per-factor
#'   summary statistics, and the first \code{n} rows) to the console.
#' @export
print.factorana_factorscores <- function(x, n = 10, ...) {
  cat("Factor Score Estimates\n")
  cat("======================\n\n")

  # Count factors
  factor_cols <- grep("^factor_[0-9]+$", names(x), value = TRUE)
  n_factors <- length(factor_cols)
  n_obs <- nrow(x)
  n_converged <- sum(x$converged, na.rm = TRUE)

  cat(sprintf("Observations: %d\n", n_obs))
  cat(sprintf("Factors: %d\n", n_factors))
  cat(sprintf("Converged: %d (%.1f%%)\n\n", n_converged, 100 * n_converged / n_obs))

  # Summary statistics for each factor
  cat("Summary Statistics:\n")
  for (k in seq_len(n_factors)) {
    fac_col <- sprintf("factor_%d", k)
    se_col <- sprintf("se_factor_%d", k)
    fac_vals <- x[[fac_col]][x$converged]
    se_vals <- x[[se_col]][x$converged]

    cat(sprintf("  Factor %d: mean=%.3f, sd=%.3f, mean_se=%.3f\n",
                k, mean(fac_vals, na.rm = TRUE),
                sd(fac_vals, na.rm = TRUE),
                mean(se_vals, na.rm = TRUE)))
  }

  cat(sprintf("\nFirst %d rows:\n", min(n, n_obs)))
  print.data.frame(head(x, n))

  invisible(x)
}
