## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(factorana)

set.seed(1)
n <- 300

# Latent factors (true values, unobserved in practice)
f_cog    <- rnorm(n, mean = 0, sd = 1.0)
f_noncog <- rnorm(n, mean = 0, sd = 0.8)

# Cognitive indicators: loadings (1.0, 0.9, 0.7); error sd = 0.5
cog1 <- 0.0 + 1.0 * f_cog + rnorm(n, 0, 0.5)
cog2 <- 0.2 + 0.9 * f_cog + rnorm(n, 0, 0.5)
cog3 <- 0.1 + 0.7 * f_cog + rnorm(n, 0, 0.5)

# Non-cognitive indicators: loadings (1.0, 1.1, 0.8); error sd = 0.5
nc1 <- 0.0 + 1.0 * f_noncog + rnorm(n, 0, 0.5)
nc2 <- 0.1 + 1.1 * f_noncog + rnorm(n, 0, 0.5)
nc3 <- 0.0 + 0.8 * f_noncog + rnorm(n, 0, 0.5)

dat <- data.frame(
  intercept = 1,
  cog1 = cog1, cog2 = cog2, cog3 = cog3,
  nc1  = nc1,  nc2  = nc2,  nc3  = nc3
)
head(dat)

## -----------------------------------------------------------------------------
fm <- define_factor_model(n_factors = 2, factor_structure = "independent")

## -----------------------------------------------------------------------------
# Cognitive indicators: load on factor 1 only
mc_cog1 <- define_model_component(
  name = "cog1", data = dat, outcome = "cog1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(1, 0)        # factor 1 loading = 1, factor 2 loading = 0
)
mc_cog2 <- define_model_component(
  name = "cog2", data = dat, outcome = "cog2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(NA_real_, 0) # factor 1 loading free, factor 2 loading = 0
)
mc_cog3 <- define_model_component(
  name = "cog3", data = dat, outcome = "cog3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(NA_real_, 0)
)

# Non-cognitive indicators: load on factor 2 only
mc_nc1 <- define_model_component(
  name = "nc1", data = dat, outcome = "nc1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, 1)
)
mc_nc2 <- define_model_component(
  name = "nc2", data = dat, outcome = "nc2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, NA_real_)
)
mc_nc3 <- define_model_component(
  name = "nc3", data = dat, outcome = "nc3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = c(0, NA_real_)
)

## -----------------------------------------------------------------------------
ms <- define_model_system(
  components = list(mc_cog1, mc_cog2, mc_cog3, mc_nc1, mc_nc2, mc_nc3),
  factor = fm
)

## -----------------------------------------------------------------------------
ctrl <- define_estimation_control(n_quad_points = 6, num_cores = 1)

fit <- estimate_model_rcpp(
  model_system = ms,
  data         = dat,
  control      = ctrl,
  optimizer    = "nlminb",
  parallel     = FALSE,
  verbose      = FALSE
)

fit$convergence  # 0 indicates successful convergence

## -----------------------------------------------------------------------------
# Tidy table of parameter estimates with standard errors
components_table(fit, digits = 3)

## -----------------------------------------------------------------------------
fscores <- estimate_factorscores_rcpp(
  fit, dat, control = ctrl, parallel = FALSE, verbose = FALSE
)
head(fscores[, c("obs_id", "factor_1", "factor_2",
                 "se_factor_1", "se_factor_2", "converged")])

# Correlation of estimated factor scores with the true (unobserved) factors
cor(fscores$factor_1, f_cog)
cor(fscores$factor_2, f_noncog)

