#Helper functions. Right now mainly for checking data types

# Check if vector is binary (0/1)
is_binary <- function(v) {
  u <- unique(v[!is.na(v)])
  length(u) <= 2 && all(u %in% c(0, 1))
}

# Check if vector is numeric or integer
is_numeric_or_integer <- function(v) {
  is.numeric(v) || is.integer(v)
}

# Check for zero-variance columns
zero_variance_cols <- function(df) {
  names(df)[vapply(df, function(col) {
    u <- unique(col[!is.na(col)])
    length(u) <= 1
  }, logical(1))]
}

# Check for all-NA columns
all_na_cols <- function(df) {
  names(df)[vapply(df, function(col) all(is.na(col)), logical(1))]
}

# Check NA rate per column and flag columns exceeding threshold
high_na_cols <- function(df, threshold = 0.5) {
  na_rates <- vapply(df, function(col) mean(is.na(col)), numeric(1))
  names(na_rates)[na_rates > threshold]
}

# Safe NULL default operator (like `%||%` in rlang)
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Apply loading constraints (internal)
#' @keywords internal
#' @noRd
apply_loading_constraints <- function(loading, normalization) {
  stopifnot(length(loading) == length(normalization))
  fixed <- which(!is.na(normalization))
  if (length(fixed)) loading[fixed] <- normalization[fixed]
  loading
}

