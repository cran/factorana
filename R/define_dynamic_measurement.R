#' Define a dynamic measurement system for longitudinal factor models
#'
#' Builds a Stage 1 measurement model for a single latent construct
#' observed at two or more time points. Measurement invariance is imposed
#' on loadings and residual sigmas (tied across periods via
#' \code{equality_constraints}), while measurement intercepts are left
#' period-specific. This is the recommended Stage 1 setup for an SE_linear
#' or SE_quadratic structural model in the second stage.
#'
#' Why period-specific intercepts? Under the usual factor-model
#' identification convention \code{E[f_k] = 0} for every factor \code{k},
#' pooling measurement intercepts across periods biases them by
#' \eqn{\lambda_m \cdot E[f_2]/2} when the period-mean drifts between
#' waves, an artefact that propagates into an under-estimate of
#' \code{se_intercept} in Stage 2. Leaving the intercepts period-specific
#' and carrying the wave-1 intercepts into Stage 2 (via
#' \code{\link{build_dynamic_previous_stage}}) sidesteps the bias.
#'
#' @param data Wide-format data frame. Must contain a column named
#'   \code{paste0(prefix, item)} for every combination of
#'   \code{period_prefixes} and \code{items}, plus any columns named in
#'   \code{covariates} and \code{evaluation_indicator}.
#' @param items Character vector of item names (e.g.,
#'   \code{c("m1", "m2", "m3")}).
#' @param period_prefixes Character vector of column prefixes, one per
#'   period. Column names are assembled as \code{paste0(prefix, item)}.
#'   E.g., \code{c("Y_t1_", "Y_t2_")} yields data columns \code{"Y_t1_m1"},
#'   \code{"Y_t2_m1"}, etc. Length of \code{period_prefixes} determines
#'   the number of latent factors in Stage 1 (one per period).
#' @param model_type One of \code{"linear"}, \code{"oprobit"},
#'   \code{"probit"}, \code{"logit"}. The same type is used for every
#'   item and every period.
#' @param n_categories Required for \code{model_type = "oprobit"}; the
#'   number of ordered categories (shared across items and periods).
#' @param covariates Character vector of covariate column names; default
#'   \code{"intercept"}. Same covariates for every component.
#' @param evaluation_indicator Name of a column with 0/1 values
#'   indicating which observations contribute to each component's
#'   likelihood; \code{NULL} to use all rows.
#'
#' @return An object of class \code{"dynamic_measurement"}: a list with
#'   \itemize{
#'     \item \code{model_system}: a \code{model_system} object ready to
#'       pass to \code{\link{estimate_model_rcpp}}.
#'     \item \code{factor_model}: the underlying \code{factor_model}.
#'     \item \code{items}, \code{period_prefixes}, \code{model_type},
#'       \code{n_categories}, \code{covariates},
#'       \code{evaluation_indicator}: the inputs, kept for use by
#'       \code{\link{build_dynamic_previous_stage}}.
#'   }
#' @seealso \code{\link{build_dynamic_previous_stage}} constructs the
#'   Stage 2 \code{previous_stage} object from the Stage 1 estimation
#'   result.
#' @examples
#' \donttest{
#' # Simulate a simple dynamic single-factor model
#' set.seed(1); n <- 500
#' f1 <- rnorm(n); eps <- rnorm(n, 0, sqrt(0.5))
#' f2 <- 0.4 + 0.6 * f1 + eps
#' dat <- data.frame(
#'   intercept = 1, eval = 1L,
#'   Y_t1_m1 = 1.5 + f1 + rnorm(n, 0, 0.7),
#'   Y_t1_m2 = 1.0 + 0.9 * f1 + rnorm(n, 0, 0.75),
#'   Y_t1_m3 = 0.8 + 1.1 * f1 + rnorm(n, 0, 0.65),
#'   Y_t2_m1 = 1.5 + f2 + rnorm(n, 0, 0.7),
#'   Y_t2_m2 = 1.0 + 0.9 * f2 + rnorm(n, 0, 0.75),
#'   Y_t2_m3 = 0.8 + 1.1 * f2 + rnorm(n, 0, 0.65)
#' )
#' dyn <- define_dynamic_measurement(
#'   data = dat, items = c("m1", "m2", "m3"),
#'   period_prefixes = c("Y_t1_", "Y_t2_"),
#'   model_type = "linear", evaluation_indicator = "eval"
#' )
#' ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
#' s1 <- estimate_model_rcpp(dyn$model_system, dat, control = ctrl,
#'                           optimizer = "nlminb", parallel = FALSE,
#'                           verbose = FALSE)
#' prev <- build_dynamic_previous_stage(dyn, s1, dat)   # Stage 2 input
#' }
#' @export
define_dynamic_measurement <- function(
    data,
    items,
    period_prefixes,
    model_type           = "linear",
    n_categories         = NULL,
    covariates           = "intercept",
    evaluation_indicator = NULL) {

  # ---- Validate ----
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is.character(items) || length(items) < 1L) {
    stop("`items` must be a non-empty character vector.")
  }
  if (!is.character(period_prefixes) || length(period_prefixes) < 2L) {
    stop("`period_prefixes` must be a character vector of length >= 2.")
  }
  allowed_types <- c("linear", "oprobit", "probit", "logit")
  if (!model_type %in% allowed_types) {
    stop("`model_type` must be one of: ", paste(allowed_types, collapse = ", "))
  }
  if (model_type == "oprobit" &&
      (is.null(n_categories) || !is.numeric(n_categories) || n_categories < 3)) {
    stop("`n_categories` must be an integer >= 3 for model_type = 'oprobit'.")
  }

  # Check all required columns exist
  needed <- unlist(lapply(period_prefixes, function(p) paste0(p, items)))
  missing_cols <- setdiff(needed, colnames(data))
  if (length(missing_cols) > 0L) {
    stop("Data is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!is.null(evaluation_indicator) &&
      !evaluation_indicator %in% colnames(data)) {
    stop("evaluation_indicator '", evaluation_indicator, "' not found in data.")
  }

  n_periods <- length(period_prefixes)
  n_items   <- length(items)

  # Oprobit absorbs its intercept into the cutpoint thresholds, so
  # factorana rejects an "intercept" covariate for oprobit components.
  # The wrapper's default covariates = "intercept" is correct for linear
  # / probit / logit; strip it for oprobit so the default works for all
  # model types without requiring the user to tailor covariates.
  if (model_type == "oprobit") {
    covariates <- setdiff(covariates, c("intercept", "constant"))
    if (length(covariates) == 0L) covariates <- NULL
  }

  # ---- Build factor model: one factor per period ----
  fm <- define_factor_model(
    n_factors        = n_periods,
    n_types          = 1L,
    factor_structure = "independent"
  )

  # ---- Build components ----
  # For each period p (1..n_periods) and each item i:
  #   - component name = paste0(period_prefixes[p], items[i])
  #   - outcome column = same
  #   - loading_normalization: zero on every factor except p; on factor p,
  #     fixed to 1 for i=1, free (NA_real_) otherwise.
  components <- list()
  for (p in seq_len(n_periods)) {
    for (i in seq_len(n_items)) {
      norm <- rep(0, n_periods)
      norm[p] <- if (i == 1L) 1 else NA_real_

      comp_name <- paste0(period_prefixes[p], items[i])

      comp_args <- list(
        name                  = comp_name,
        data                  = data,
        outcome               = comp_name,
        factor                = fm,
        covariates            = covariates,
        model_type            = model_type,
        loading_normalization = norm,
        use_types             = FALSE
      )
      if (!is.null(evaluation_indicator)) {
        comp_args$evaluation_indicator <- evaluation_indicator
      }
      if (model_type == "oprobit") {
        comp_args$num_choices <- n_categories
      }

      components[[length(components) + 1L]] <- do.call(
        define_model_component, comp_args
      )
    }
  }

  # ---- Build equality constraints ----
  # Tie shared measurement parameters across factors. Which parameter
  # names are free to tie depends on model_type:
  #   linear  : loadings (from item 2 on), intercepts excluded, sigmas
  #   probit  : loadings, NO sigma (fixed = 1)
  #   oprobit : loadings, thresholds (n_categories - 1 per item)
  #   logit   : loadings (binary) or per-choice loadings (multi)
  # We tie: loadings (for items 2..n_items) and sigmas/thresholds as
  # appropriate. Intercepts are never tied (period-specific by design).
  #
  # Parameter names follow the factorana convention:
  #   <comp_name>_loading_<factor_idx>
  #   <comp_name>_sigma         (linear, logit's binary sigma N/A)
  #   <comp_name>_thresh_<k>    (oprobit; k = 1..n_categories-1)
  eq_constraints <- list()

  # Tie loadings across periods for items 2..n_items (item 1's loading
  # is fixed to 1 on its factor slot).
  if (n_items >= 2L) {
    for (i in 2:n_items) {
      tied <- character(n_periods)
      for (p in seq_len(n_periods)) {
        comp_name <- paste0(period_prefixes[p], items[i])
        tied[p]   <- paste0(comp_name, "_loading_", p)
      }
      eq_constraints[[length(eq_constraints) + 1L]] <- tied
    }
  }

  # Tie sigmas (linear) or thresholds (oprobit) for every item across
  # periods. Linear: every item has a sigma parameter. Probit/logit
  # (binary): no sigma parameter (scale fixed). Oprobit: threshold
  # parameters instead of sigma.
  if (model_type == "linear") {
    for (i in seq_len(n_items)) {
      tied <- character(n_periods)
      for (p in seq_len(n_periods)) {
        comp_name <- paste0(period_prefixes[p], items[i])
        tied[p]   <- paste0(comp_name, "_sigma")
      }
      eq_constraints[[length(eq_constraints) + 1L]] <- tied
    }
  } else if (model_type == "oprobit") {
    # Tie THRESHOLD INCREMENTS (k = 2..n_categories-1) across periods, and
    # leave the FIRST threshold period-specific. factorana parameterises
    # ordered-probit cutpoints as an increment vector:
    #   cutpoint_k = thresh_1 + thresh_2 + ... + thresh_k
    # so thresh_1 is the location (analog of a linear intercept) and the
    # later increments are the category spacing (analog of the scale).
    # Tying only the increments mirrors the linear strategy of tying sigmas
    # while leaving intercepts period-specific, so the wave-specific
    # population-mean shift in the latent factor can be absorbed by
    # wave-specific thresh_1 values. Tying every cutpoint (the earlier
    # behaviour) forced the factor variances to contort to fit the mean
    # drift and led to conv=1 in Stage 1 whenever the population mean
    # shifted across periods.
    n_thresh <- as.integer(n_categories) - 1L
    if (n_thresh >= 2L) {
      for (i in seq_len(n_items)) {
        for (k in 2:n_thresh) {
          tied <- character(n_periods)
          for (p in seq_len(n_periods)) {
            comp_name <- paste0(period_prefixes[p], items[i])
            tied[p]   <- paste0(comp_name, "_thresh_", k)
          }
          eq_constraints[[length(eq_constraints) + 1L]] <- tied
        }
      }
    }
  }
  # probit / logit binary: no sigma or threshold ties

  # ---- Build model_system ----
  ms <- define_model_system(
    components           = components,
    factor               = fm,
    equality_constraints = eq_constraints
  )

  structure(
    list(
      model_system         = ms,
      factor_model         = fm,
      items                = items,
      period_prefixes      = period_prefixes,
      model_type           = model_type,
      n_categories         = n_categories,
      covariates           = covariates,
      evaluation_indicator = evaluation_indicator,
      equality_constraints = eq_constraints
    ),
    class = "dynamic_measurement"
  )
}


#' Build a Stage 2 previous_stage object from a dynamic measurement fit
#'
#' Constructs a dummy \code{previous_stage} result that plugs the
#' anchor-period measurement intercepts into every factor slot, pairs
#' them with the (tied) shared loadings and residual sigmas / thresholds,
#' and is ready to pass as \code{previous_stage} to a Stage 2
#' \code{SE_linear} or \code{SE_quadratic} model via
#' \code{\link{define_model_system}}.
#'
#' @param dyn A \code{dynamic_measurement} object from
#'   \code{\link{define_dynamic_measurement}}.
#' @param stage1_result The result object from
#'   \code{estimate_model_rcpp(dyn$model_system, ...)}.
#' @param data The same data frame passed to
#'   \code{\link{define_dynamic_measurement}} (needed to rebuild the
#'   dummy model system's components; components do not retain the
#'   data they were defined on).
#' @param anchor_period Integer index into \code{dyn$period_prefixes}
#'   giving the period whose measurement intercepts should be carried
#'   into Stage 2. The recommended choice is 1: under
#'   \code{E[f_1] = 0} this period's intercepts identify the true DGP
#'   intercepts. Default: 1.
#'
#' @return A list suitable as \code{previous_stage} in
#'   \code{\link{define_model_system}}: has fields \code{model_system},
#'   \code{estimates}, \code{std_errors}, \code{convergence},
#'   \code{loglik}. Every per-component measurement parameter is the
#'   corresponding Stage 1 estimate, with the intercepts overwritten to
#'   the anchor-period values for all periods.
#' @seealso \code{\link{define_dynamic_measurement}}.
#' @export
build_dynamic_previous_stage <- function(dyn, stage1_result, data,
                                         anchor_period = 1L) {
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    stop("`data` must be a data frame (pass the same frame used in ",
         "define_dynamic_measurement()).")
  }

  if (!inherits(dyn, "dynamic_measurement")) {
    stop("`dyn` must be a dynamic_measurement object from ",
         "define_dynamic_measurement().")
  }
  if (!is.list(stage1_result) ||
      !all(c("estimates", "std_errors") %in% names(stage1_result))) {
    stop("`stage1_result` must have components `estimates` and `std_errors` ",
         "(the return value from estimate_model_rcpp()).")
  }
  anchor_period <- as.integer(anchor_period)
  n_periods <- length(dyn$period_prefixes)
  if (is.na(anchor_period) || anchor_period < 1L ||
      anchor_period > n_periods) {
    stop(sprintf("`anchor_period` must be an integer in 1..%d.", n_periods))
  }

  s1_est <- stage1_result$estimates
  items  <- dyn$items
  prefixes <- dyn$period_prefixes
  n_items  <- length(items)
  mt <- dyn$model_type
  n_cat <- dyn$n_categories

  # Pull anchor-period intercepts (per item). Oprobit and probit/logit
  # absorb the intercept into the cutpoints / link, so there are no
  # explicit "_intercept" parameters; in that case anchor_intercepts is
  # left empty and the per-component assembly below skips intercept slots.
  anchor_pref <- prefixes[anchor_period]
  has_intercept_param <- mt == "linear"
  if (has_intercept_param) {
    anchor_intercepts <- vapply(items, function(it) {
      nm <- paste0(anchor_pref, it, "_intercept")
      if (!nm %in% names(s1_est)) {
        stop("Stage 1 result missing expected parameter: ", nm)
      }
      unname(s1_est[nm])
    }, numeric(1))
    names(anchor_intercepts) <- items
  } else {
    anchor_intercepts <- setNames(numeric(0), character(0))
  }

  # Pull tied loadings (items 2..n_items, value equal across all periods)
  tied_loadings <- if (n_items >= 2L) {
    vapply(items[-1L], function(it) {
      # Any period works because they are tied; use anchor.
      nm <- paste0(anchor_pref, it, "_loading_", anchor_period)
      unname(s1_est[nm])
    }, numeric(1))
  } else {
    numeric(0)
  }
  if (length(tied_loadings) > 0L) names(tied_loadings) <- items[-1L]

  # Pull tied sigmas (linear) or thresholds (oprobit)
  tied_sigmas <- NULL
  tied_thresh <- NULL
  if (mt == "linear") {
    tied_sigmas <- vapply(items, function(it) {
      nm <- paste0(anchor_pref, it, "_sigma")
      unname(s1_est[nm])
    }, numeric(1))
    names(tied_sigmas) <- items
  } else if (mt == "oprobit") {
    # For every item, pull the anchor-period threshold vector. Increments
    # (k >= 2) are tied across periods in Stage 1; thresh_1 is period-
    # specific. Here we carry the anchor period's thresh_1 into every
    # period's slot so that Stage 2 inherits a common location (analog of
    # the linear wrapper carrying wave-1 intercepts).
    n_thresh <- as.integer(n_cat) - 1L
    tied_thresh <- lapply(items, function(it) {
      vapply(seq_len(n_thresh), function(k) {
        nm <- paste0(anchor_pref, it, "_thresh_", k)
        if (!nm %in% names(s1_est)) {
          stop("Stage 1 result missing expected parameter: ", nm)
        }
        unname(s1_est[nm])
      }, numeric(1))
    })
    names(tied_thresh) <- items
  }

  # Pull covariate coefficients (other than intercept): these are NOT
  # tied across periods in this wrapper; we simply pick anchor-period
  # values per component. In most use cases `covariates = "intercept"`
  # only and this branch does nothing.
  non_intercept_covs <- setdiff(dyn$covariates, "intercept")

  # ---- Assemble combined estimate vector in canonical order expected
  #      by build_parameter_metadata(): factor_var_1..factor_var_K, then
  #      per-component (intercept, [loadings], [cov-betas], [sigma or
  #      thresholds]).
  vals  <- numeric(0)
  names_ <- character(0)

  for (k in seq_len(n_periods)) {
    vals   <- c(vals,   unname(s1_est[paste0("factor_var_", k)]))
    names_ <- c(names_, paste0("factor_var_", k))
  }

  for (p in seq_len(n_periods)) {
    for (i in seq_len(n_items)) {
      comp_name <- paste0(prefixes[p], items[i])

      # Intercept: anchor period's value for EVERY period. Skip for model
      # types that do not have an explicit intercept (oprobit / probit /
      # logit), where the location is absorbed into cutpoints / link.
      if (has_intercept_param) {
        vals   <- c(vals,   anchor_intercepts[[items[i]]])
        names_ <- c(names_, paste0(comp_name, "_intercept"))
      }

      # Free loading (item >= 2) on factor p
      if (i >= 2L) {
        vals   <- c(vals,   tied_loadings[[items[i]]])
        names_ <- c(names_, paste0(comp_name, "_loading_", p))
      }

      # Covariate betas (other than intercept): anchor-period values
      for (cv in non_intercept_covs) {
        nm <- paste0(prefixes[anchor_period], items[i], "_", cv)
        if (nm %in% names(s1_est)) {
          vals   <- c(vals,   unname(s1_est[nm]))
          names_ <- c(names_, paste0(comp_name, "_", cv))
        }
      }

      # Sigma / thresholds
      if (mt == "linear") {
        vals   <- c(vals,   tied_sigmas[[items[i]]])
        names_ <- c(names_, paste0(comp_name, "_sigma"))
      } else if (mt == "oprobit") {
        for (k in seq_along(tied_thresh[[items[i]]])) {
          vals   <- c(vals,   tied_thresh[[items[i]]][k])
          names_ <- c(names_, paste0(comp_name, "_thresh_", k))
        }
      }
    }
  }
  names(vals) <- names_

  # Rebuild a clean model_system for the dummy (independent k-factor,
  # no equality_constraints: Stage 2 will fix these via previous_stage).
  fm_dummy <- define_factor_model(
    n_factors        = n_periods,
    n_types          = 1L,
    factor_structure = "independent"
  )
  dummy_data <- data
  comps_dummy <- list()
  for (p in seq_len(n_periods)) {
    for (i in seq_len(n_items)) {
      norm <- rep(0, n_periods)
      norm[p] <- if (i == 1L) 1 else NA_real_
      comp_args <- list(
        name                  = paste0(prefixes[p], items[i]),
        data                  = dummy_data,
        outcome               = paste0(prefixes[p], items[i]),
        factor                = fm_dummy,
        covariates            = dyn$covariates,
        model_type            = mt,
        loading_normalization = norm,
        use_types             = FALSE
      )
      if (!is.null(dyn$evaluation_indicator)) {
        comp_args$evaluation_indicator <- dyn$evaluation_indicator
      }
      if (mt == "oprobit") {
        comp_args$num_choices <- n_cat
      }
      comps_dummy[[length(comps_dummy) + 1L]] <- do.call(
        define_model_component, comp_args
      )
    }
  }
  ms_dummy <- define_model_system(components = comps_dummy, factor = fm_dummy)

  list(
    model_system = ms_dummy,
    estimates    = vals,
    std_errors   = setNames(rep(0, length(vals)), names_),
    convergence  = 0L,
    loglik       = 0.0
  )
}
