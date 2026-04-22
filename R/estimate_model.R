
# ------- Helper function: initialize parameters for a single model component ---------

# Internal: Initialize parameters for a single model component
# Runs regression without factors and sets initial parameter values
# @param mc A model_component object
# @return a list of initialized parameters
init_single_component <- function(mc) {
  df0 <- as.data.frame(mc$data)

  # ---- 1. Apply evaluation indicator & drop missing outcomes ----
  # Keep only rows where evaluation_indicator == TRUE / 1, and drop rows with missing outcome.

  if (!is.null(mc$evaluation_indicator)) {
    ei <- df0[[mc$evaluation_indicator]]
    if (is.null(ei)) stop("evaluation_indicator '", mc$evaluation_indicator, "' not found in data.")
    if (is.logical(ei)) {
      keep <- !is.na(ei) & ei
    } else if (is.numeric(ei) || is.integer(ei)) {
      keep <- !is.na(ei) & (ei == 1L)
    } else {
      stop("`evaluation_indicator` must be logical or 0/1 numeric.")
    }
  } else {
    keep <- rep(TRUE, nrow(df0))
  }

  # Check that the outcome exists and is non-missing
  # For exploded logit (multiple outcomes), use first outcome column
  outcome_col <- if (length(mc$outcome) > 1) mc$outcome[1] else mc$outcome
  if (!outcome_col %in% names(df0)) stop("Outcome '", outcome_col, "' not in data.")
  keep <- keep & !is.na(df0[[outcome_col]])
  df0 <- df0[keep, , drop = FALSE]
  if (nrow(df0) == 0) stop("No rows to estimate for ", mc$name, " after conditioning on eval & non-missing outcome.")
  # -----------------------------------------------------------

  # ---- 2. Prepare regression data ----
  # Build y and X matrices, combining into a clean df for fitting
  # For exploded logit, use first outcome for initialization purposes

  y <- df0[[outcome_col]]
  covars <- unlist(mc$covariates, use.names = FALSE)
  X <- if (length(covars)) df0[, covars, drop = FALSE] else NULL
  df <- if (is.null(X)) data.frame(y = y) else data.frame(y = y, X, check.names = FALSE)

  # ---- 3. Retrieve model & factor information ----
  model_type <- mc$model_type

  k <- if (!is.null(mc$k)) mc$k else as.integer(mc$factor$n_factors)
  if (is.na(k) || k < 1L) stop("initialize_parameters: invalid k (number of factors).")

  # ---- 4. Normalize factor loadings ----
  # Apply normalization constraints: NA = free, numeric = fixed value.

  norm_vec <- mc$loading_normalization
  if (!is.numeric(norm_vec) || length(norm_vec) != k) {
    stop("initialize_parameters: loading_normalization must be numeric length k.")
  }

  # k-length default loadings, then apply constraints (NA = free; numeric = fixed)
  init_loading <- rep(0.3, k)
  fixed_idx <- which(!is.na(norm_vec))
  if (length(fixed_idx)) init_loading[fixed_idx] <- norm_vec[fixed_idx]


  # ---- 5. Fit model-type specific regression ----
  # Estimate initial intercepts/betas via simple models, ignoring latent factors.
  out <- NULL

  if (model_type == "linear") {
    fit <- lm(y ~ ., data = df)
    coefs <- coef(fit)
    out <- list(
      intercept = unname(coefs[1]),
      betas     = unname(coefs[-1])
#      loading   = 0.1 * sd(y)
    )

  } else if (model_type == "probit") {
    fit <- glm(y ~ ., data = df, family = binomial(link = "probit"))
    coefs <- coef(fit)
    out <- list(
      intercept = unname(coefs[1]),
      betas     = unname(coefs[-1])
 #     loading   = 0.1
    )

  } else if (model_type == "logit") {
    # (If this is truly binary logit, glm(binomial(link="logit")) is more typical.)
    fit <- nnet::multinom(y ~ ., data = df, trace = FALSE)
    coefs <- coef(fit)
    out <- list(
      intercept = unname(coefs[1]),
      betas     = unname(coefs[-1])
 #     loading   = 0.1
    )

  } else if (model_type == "oprobit") {
    # use df0 (conditioned), not mc$data
    # create threshold vector evenly spaced in [-1, 1]
    y_op <- df0[[ mc$outcome ]]
    if (!is.ordered(y_op)) {
      stop("initialize_parameters(oprobit): outcome must be an ordered factor. Got: ",
           paste(class(y_op), collapse = "/"))
    }
    n_cats <- nlevels(y_op)
    if (n_cats < 3L) {
      stop("Ordered probit needs >= 3 categories; got ", n_cats,
           ". Did the evaluation subset collapse categories?")
    }
    n_thresh <- n_cats - 1L

    out <- list(
      intercept  = 0,
      betas      = rep(0, length(mc$covariates)),
 #     loading    = 1.0,
      thresholds = seq(-1, 1, length.out = n_thresh)
    )

  } else {
    stop("Unsupported model type: ", model_type)
  }

  # ---- 6. Append defaults common to all components ----
  # Add initialized loadings, factor variance, and correlation placeholders.
  out$loading <- init_loading
  out$factor_var <- 1
  out$factor_cor <- 0


  return(out)
}


# ---- Estimate model (currently just initializes parameters) ----

#' Estimate model
#'
#' @param ms an object of class model_system
#' @param control an object of class estimation_control
#'
#' @return description

estimate_model <- function(ms, control){

  stopifnot(inherits(ms, "model_system"))
  stopifnot(inherits(control, "estimation_control"))

  factor <- ms$factor
  components <-ms$components

  #initialize parameters for each component
  inits <- lapply(components, initialize_parameters)

  #EDIT: flatten matrices
  init_df <- do.call(rbind, lapply(seq_along(inits), function(i) {
    comp <- components[[i]]
    init <- inits[[i]]

    intercept   <- if (!is.null(init$intercept)) init$intercept else NA
    betas_str   <- if (!is.null(init$betas)) paste(init$betas, collapse = ";") else ""
    loading_str <- if (!is.null(init$loading)) paste(init$loading, collapse = ";") else ""
    fvar_str    <- if (!is.null(init$factor_var)) paste(init$factor_var, collapse = ";") else ""
    fcor_str    <- if (!is.null(init$factor_cor)) paste(as.vector(init$factor_cor), collapse = ";") else ""

    data.frame(
      component = comp$name,
      intercept = intercept,
      betas     = betas_str,
      loading   = loading_str,
      factor_var = fvar_str,
      factor_cor = fcor_str,
      stringsAsFactors = FALSE
    )
  }))

  # Print to console for now
  print(init_df)

  # Later: write.csv/init_df, pass to optimizer
  invisible(init_df)
}

# ---- utilities ----------------------(might need to move this to the utils.R file)

#' @keywords internal
.to_chr <- function(x) {
  # Normalize scalars to character for CSV; leave JSON strings alone.
  if (is.null(x) || length(x) == 0) return("")        # <- fix so can handle NAs
  if (length(x) > 1) return(paste0(x, collapse = ","))# only used for scalars here, joins multilength values w commas
  if (is.logical(x)) return(ifelse(x, "TRUE", "FALSE")) #normalize type to string
  if (is.numeric(x)) return(as.character(x)) #normalize type to string
  if (is.character(x)) return(x)
  # fallback for other scalars
  paste0(x)
}

#' @keywords internal
.to_json <- function(x) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    # minimal fallback if jsonlite isn't installed
    return(paste(x, collapse = "|"))
  }
  jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")
}

# --- helper: filter rows like in estimation ---
.eval_keep <- function(df, eval_indicator, outcome) {
  keep <- if (is.null(eval_indicator)) rep(TRUE, nrow(df)) else {
    ei <- df[[eval_indicator]]
    if (is.logical(ei)) !is.na(ei) & ei else (!is.na(ei) & (ei == 1L))
  }
  keep & !is.na(df[[outcome]])
}

# --- compute SEs for ONE component by refitting without factors ---
compute_se_for_component <- function(mc) {
  df0 <- as.data.frame(mc$data)
  keep <- .eval_keep(df0, mc$evaluation_indicator, mc$outcome)
  df0 <- df0[keep, , drop = FALSE]
  if (nrow(df0) == 0) return(list(se_intercept = NA_real_, se_betas = numeric(0), se_loading = NA_real_))

  y <- df0[[mc$outcome]]
  covars <- unlist(mc$covariates, use.names = FALSE)
  X <- if (length(covars)) df0[, covars, drop = FALSE] else NULL
  df <- if (is.null(X)) data.frame(y = y) else data.frame(y = y, X, check.names = FALSE)

  mt <- mc$model_type
  if (mt == "linear") {
    fit <- stats::lm(y ~ ., data = df)
    sevec <- summary(fit)$coefficients[, "Std. Error"]
  } else if (mt == "probit") {
    fit <- stats::glm(y ~ ., data = df, family = binomial(link = "probit"))
    sevec <- summary(fit)$coefficients[, "Std. Error"]
  } else if (mt == "logit") {
    fit <- nnet::multinom(y ~ ., data = df, trace = FALSE)
    seobj <- summary(fit)$standard.errors
    sevec <- if (is.matrix(seobj)) seobj[1, ] else as.numeric(seobj)
  } else {
    # not implemented here (e.g., oprobit): return NAs
    return(list(se_intercept = NA_real_, se_betas = numeric(0), se_loading = NA_real_))
  }

  list(
    se_intercept = unname(sevec[1]),
    se_betas     = if (length(sevec) > 1) unname(sevec[-1]) else numeric(0),
    se_loading   = NA_real_  # loading SE not computed in this simple refit
  )
}

pack_values_with_ses <- function(ms, inits, factor_var_first = 1.0) {
  vals <- numeric(0)
  ses  <- numeric(0)

  # 1) factor variance first
  vals <- c(vals, factor_var_first)
  ses  <- c(ses,  NA_real_)   # no SE for fixed normalization

  # 2) for each component: intercept, betas, loading
  for (i in seq_along(inits)) {
    init <- inits[[i]]
    mc   <- ms$components[[i]]
    sei  <- compute_se_for_component(mc)

    # intercept
    vals <- c(vals, unname(init$intercept))
    ses  <- c(ses,  sei$se_intercept)

    # betas
    if (length(init$betas)) {
      vals <- c(vals, unname(init$betas))
      # align lengths safely
      if (length(sei$se_betas) == length(init$betas)) {
        ses <- c(ses, sei$se_betas)
      } else {
        ses <- c(ses, rep(NA_real_, length(init$betas)))
      }
    }

    # loading (SE not computed here)
    vals <- c(vals, unname(init$loading))
    ses  <- c(ses, rep(sei$se_loading, length(init$loading)))
  }

  # --- sanity check ---
  if (length(vals) != length(ses)) {
    stop("pack_values_with_ses: mismatch between values (", length(vals),
         ") and ses (", length(ses), ")")
  }

  list(values = vals, ses = ses)
}

#
# # --- pack values + computed SEs; factor variance first ---
# pack_values_with_ses <- function(ms, inits, factor_var_first = 1.0) {
#   vals <- numeric(0); ses <- numeric(0)
#
#   # 1) factor variance first
#
#   vals <- c(vals, unname(ini$loading))
#   ses  <- c(ses, rep(sei$se_loading, length(ini$loading)))
#
#   # vals <- c(vals, factor_var_first)
#   # ses  <- c(ses,  NA_real_)   # no SE for fixed normalization
#
#   # 2) for each component: intercept, betas, loading
#   for (i in seq_along(inits)) {
#     ini <- inits[[i]]
#     mc  <- ms$components[[i]]
#     sei <- compute_se_for_component(mc)
#
#     # intercept
#     vals <- c(vals, unname(ini$intercept))
#     ses  <- c(ses,  sei$se_intercept)
#
#     # betas
#     if (length(ini$betas)) {
#       vals <- c(vals, unname(ini$betas))
#       # align lengths safely
#       if (length(sei$se_betas) == length(ini$betas)) {
#         ses <- c(ses, sei$se_betas)
#       } else {
#         ses <- c(ses, rep(NA_real_, length(ini$betas)))
#       }
#     }
#
#     # loading (SE not computed here)
#     vals <- c(vals, unname(ini$loading))
#     ses  <- c(ses,  sei$se_loading)
#   }
#
#   list(values = vals, ses = ses)
# }

# --- writer unchanged ---
write_meas_par <- function(values, ses, path) {
  if (missing(path) || is.null(path)) {
    stop("'path' is required; pass an explicit file path.")
  }
  stopifnot(length(values) == length(ses))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  idx <- seq_along(values) - 1L
  tab <- data.frame(idx, values = as.numeric(values), se = as.numeric(ses))
  utils::write.table(tab, file = path, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(path)
}


# ---- generic: flatten to key-value rows -----------------------------------

#' Convert object to key-value format (internal)
#' @param x Object to convert (an \code{estimation_control}, \code{factor_model},
#'   \code{model_component}, or \code{model_system}).
#' @param ... Additional arguments (not used).
#' @return A data frame of key-value rows with columns \code{section},
#'   \code{component}, \code{key}, \code{value}, and \code{dtype}, one row per
#'   configuration field. Used internally by \code{\link{write_model_config_csv}}.
#' @keywords internal
as_kv <- function(x, ...) UseMethod("as_kv")


# which fields define factor-model identity
.factor_keys <- c("n_factors","n_types","correlation","n_mixtures","nfac_param")

# ---- estimation_control ----------------------------------------------------

#' @rdname as_kv
#' @exportS3Method
as_kv.estimation_control <- function(x, ...) {
  # adjust keys to whatever your estimation_control actually stores in future if change
  keys   <- c("n_quad_points", "num_cores")
  values <- c(.to_chr(x$n_quad_points), .to_chr(x$num_cores))
  dtype  <- c("int", "int")
  data.frame(
    section   = "estimation_control",
    component = "",
    key = keys, value = values, dtype = dtype,
    stringsAsFactors = FALSE
  )
}

# ---- factor_model ----------------------------------------------------------

#' @rdname as_kv
#' @exportS3Method
as_kv.factor_model <- function(x, ...) {
  keys <- c("n_factors","n_types","correlation","n_mixtures","nfac_param")
  vals <- c(.to_chr(x$n_factors),
            .to_chr(x$n_types),
            .to_chr(x$correlation),
            .to_chr(x$n_mixtures),
            .to_chr(x$nfac_param))
  dtype <- c("int","int","bool","int","int")
  data.frame(
    section   = "factor_model",
    component = "",
    key = keys, value = vals, dtype = dtype,
    stringsAsFactors = FALSE
  )
}

# ---- model_component -------------------------------------------------------

#' @rdname as_kv
#' @exportS3Method
as_kv.model_component <- function(x, ...) {
  comp <- if (!is.null(x$name) && nzchar(x$name) && !is.na(x$name)) x$name else ""
  rows <- list(
    data.frame(section="model_component", component=comp, key="name",
               value=.to_chr(comp), dtype="string"),
    data.frame(section="model_component", component=comp, key="outcome",
               value=.to_chr(x$outcome), dtype="string"),
    data.frame(section="model_component", component=comp, key="model_type",
               value=.to_chr(x$model_type), dtype="string"),
    data.frame(section="model_component", component=comp, key="intercept",
               value=.to_chr(x$intercept), dtype="bool"),
    data.frame(section="model_component", component=comp, key="num_choices",
               value=.to_chr(x$num_choices), dtype="int"),
    data.frame(section="model_component", component=comp, key="nrank",
               value=.to_chr(if (is.null(x$nrank)) "" else x$nrank), dtype="int"),
    data.frame(section="model_component", component=comp, key="evaluation_indicator",
               value=.to_chr(if (is.null(x$evaluation_indicator)) "" else x$evaluation_indicator), dtype="string"),
    data.frame(section="model_component", component=comp, key="covariates",
               value=.to_json(x$covariates), dtype="json"),
    data.frame(section="model_component", component=comp, key="nparam_model",
               value=.to_chr(x$nparam_model), dtype="int")
  )
  do.call(rbind, rows)
}

# ---- model_system aggregator + writer -------------------------------------

#' @rdname as_kv
#' @exportS3Method
as_kv.model_system <- function(x, ...) {
  stopifnot(inherits(x, "model_system"))
  comps <- x$components
  comp_names <- names(comps)

  sys_rows <- data.frame(
    section   = "model_system",
    component = "",
    key       = c("n_components", "component_names"),
    value     = c(.to_chr(length(comps)), .to_json(comp_names)),
    dtype     = c("int", "json"),
    stringsAsFactors = FALSE
  )

  comp_rows <- do.call(
    rbind,
    lapply(comps, function(mc) as_kv(mc))
  )

  # if no components, comp_rows is NULL
  if (is.null(comp_rows)) return(sys_rows)
  rbind(sys_rows, comp_rows)
}

#' Write a single CSV with all configuration rows
#' @param model_system a model_system object
#' @param factor_model a factor_model object
#' @param estimation_control an estimation_control object
#' @param file Path to CSV to write. Required; pass an explicit path (use
#'   \code{tempdir()} in examples or tests).
#' @return Invisibly returns the data frame of configuration rows that was
#'   written to \code{file} (columns \code{section}, \code{component},
#'   \code{key}, \code{value}, \code{dtype}). Called primarily for its side
#'   effect of writing a CSV.
#' @export
write_model_config_csv <- function(model_system, factor_model, estimation_control, file) {
  if (missing(file) || is.null(file)) {
    stop("'file' is required; pass an explicit file path.")
  }
  stopifnot(inherits(model_system, "model_system"),
            inherits(factor_model, "factor_model"),
            inherits(estimation_control, "estimation_control"))

  df <- rbind(
    as_kv(estimation_control),
    as_kv(factor_model),
    as_kv(model_system)
  )
  # Write; keep strings as-is
  utils::write.csv(df, file = file, row.names = FALSE, na = "")
  invisible(df)
}
