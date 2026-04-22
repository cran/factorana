#' Run estimation and write standard output files
#'
#' Calls estimate_model(), then writes the four expected output files:
#' - model_config.csv
#' - meas_par.csv
#' - system_inits_long.csv
#' - simulated_data.csv
#'
#' @param model_system model_system object
#' @param factor_model factor_model object
#' @param control estimation_control object
#' @param data Optional data frame (to save as simulated_data.csv)
#' @param results_dir Directory for writing output files. This argument is
#'   required; the function writes nothing without an explicit path. Use
#'   \code{tempdir()} in examples or tests.
#' @return Invisibly returns a list with two components: \code{results} (a
#'   data frame of initial parameter values per component, as produced by
#'   \code{estimate_model()}) and \code{packed} (a list with \code{values} and
#'   \code{ses} giving the packed parameter vector and standard errors written
#'   to \code{meas_par.csv}). Called primarily for the side effect of writing
#'   \code{model_config.csv}, \code{meas_par.csv}, \code{system_inits_long.csv},
#'   and optionally \code{simulated_data.csv} into \code{results_dir}.
#' @export
estimate_and_write <- function(model_system, factor_model, control,
                               data = NULL, results_dir) {
  if (missing(results_dir) || is.null(results_dir)) {
    stop("'results_dir' is required. Pass an explicit directory path; use tempdir() in examples.")
  }
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  cat("=== Running estimation ===\n")

  # 1. Estimate model (prints init_df)
  res <- estimate_model(model_system, control)

  # 2. Collect parameter and SE info
  inits  <- lapply(model_system$components, initialize_parameters)
  packed <- pack_values_with_ses(model_system, inits, factor_var_first = 1.0)

  # 3. Write standard outputs
  write_model_config_csv(model_system, factor_model, control,
                         file.path(results_dir, "model_config.csv"))

  param_df <- data.frame(
    index = seq_along(packed$values),
    value = packed$values,
    se    = ifelse(is.na(packed$ses), 1, packed$ses)
  )
  write.csv(param_df, file.path(results_dir, "meas_par.csv"), row.names = FALSE)

  write.csv(res, file.path(results_dir, "system_inits_long.csv"), row.names = FALSE)

  if (!is.null(data))
    write.csv(data, file.path(results_dir, "simulated_data.csv"), row.names = FALSE)

  cat("Wrote all output files to", normalizePath(results_dir), "\n")
  invisible(list(results = res, packed = packed))
}
