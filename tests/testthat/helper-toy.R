test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

make_toy <- function(n = 300L, seed = 42L){
  set.seed(seed)
  f <- rnorm(n)
  X1 <- rnorm(n); X0 <- rnorm(n); Z1 <- rnorm(n)
  Y1 <- 2 + 1.0*X1 + 0.8*f + rnorm(n)
  Y0 <- 1 + 0.5*X0 + 0.3*f + rnorm(n)
  D  <- as.integer(0.3*Z1 + 0.7*f + rnorm(n) > 0)
  Y  <- ifelse(D==1, Y1, Y0)
  T  <- 0.5 + 1.2*f + rnorm(n)
  dat <- data.frame(Y, D, X1, X0, Z1, T)
  dat$eval_y1 <- as.integer(D == 1L)
  dat$eval_y0 <- as.integer(D == 0L)
  dat
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Minimal packer/writer used in output tests (placeholders for SEs)
pack_values_with_dummy_ses <- function(ms, inits, factor_var_first = 1.0, se_fill = 0.0) {
  vals <- numeric(0)
  vals <- c(vals, factor_var_first)  # factor variance first
  for (i in seq_along(inits)) {
    ini <- inits[[i]]
    vals <- c(vals,
              unname(ini$intercept),
              if (length(ini$betas)) unname(ini$betas),
              if (!is.null(ini$cutpoints) && length(ini$cutpoints)) unname(ini$cutpoints),
              unname(ini$loading))
  }
  ses <- rep(se_fill, length(vals))
  list(values = vals, ses = ses)
}

write_meas_par <- function(values, ses, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  idx <- seq_along(values) - 1L  # 0-based index
  tab <- data.frame(idx, values = as.numeric(values), se = as.numeric(ses))
  utils::write.table(tab, file = path, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(path)
}
