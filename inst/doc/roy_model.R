## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(factorana)

set.seed(108)
n <- 500

# Observed covariates and the unobserved ability factor
x1 <- rnorm(n)           # shifts wages
x2 <- rnorm(n)           # shifts wages and sector choice
f  <- rnorm(n, sd = 1)   # latent ability (unobserved in practice)

# Test scores measure ability with error (shared scale across tests)
T1 <- 2.0 + 1.0 * f + rnorm(n, 0, 0.5)
T2 <- 1.5 + 1.2 * f + rnorm(n, 0, 0.6)
T3 <- 1.0 + 0.8 * f + rnorm(n, 0, 0.4)

# Potential wages in each sector
wage0 <- 2.0 + 0.5 * x1 + 0.3 * x2 + 0.5 * f + rnorm(n, 0, 0.6)
wage1 <- 2.5 + 0.6 * x1 +             1.0 * f + rnorm(n, 0, 0.7)

# Sector choice: higher ability -> more likely to pick sector 1
z_sector <- 0.0 + 0.4 * x2 + 0.8 * f
sector <- as.numeric(runif(n) < pnorm(z_sector))

# Only the wage in the chosen sector is observed
wage <- ifelse(sector == 1, wage1, wage0)

dat <- data.frame(
  intercept = 1,
  x1 = x1, x2 = x2,
  T1 = T1, T2 = T2, T3 = T3,
  wage = wage,
  sector = sector,
  eval_tests = 1L,            # always observe all three tests
  eval_wage0 = 1L - sector,   # wage0 observed iff sector = 0
  eval_wage1 = sector,        # wage1 observed iff sector = 1
  eval_sector = 1L            # sector choice always observed
)

## -----------------------------------------------------------------------------
fm <- define_factor_model(n_factors = 1, n_types = 1)

## -----------------------------------------------------------------------------
mc_T1 <- define_model_component(
  name = "T1", data = dat, outcome = "T1", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = 1.0,
  evaluation_indicator = "eval_tests"
)
mc_T2 <- define_model_component(
  name = "T2", data = dat, outcome = "T2", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_tests"
)
mc_T3 <- define_model_component(
  name = "T3", data = dat, outcome = "T3", factor = fm,
  covariates = "intercept", model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_tests"
)

## -----------------------------------------------------------------------------
mc_wage0 <- define_model_component(
  name = "wage0", data = dat, outcome = "wage", factor = fm,
  covariates = c("intercept", "x1", "x2"), model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_wage0"
)
mc_wage1 <- define_model_component(
  name = "wage1", data = dat, outcome = "wage", factor = fm,
  covariates = c("intercept", "x1"), model_type = "linear",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_wage1"
)

## -----------------------------------------------------------------------------
mc_sector <- define_model_component(
  name = "sector", data = dat, outcome = "sector", factor = fm,
  covariates = c("intercept", "x2"), model_type = "probit",
  loading_normalization = NA_real_,
  evaluation_indicator = "eval_sector"
)

## -----------------------------------------------------------------------------
ms <- define_model_system(
  components = list(mc_T1, mc_T2, mc_T3, mc_wage0, mc_wage1, mc_sector),
  factor     = fm
)

## -----------------------------------------------------------------------------
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)

fit <- estimate_model_rcpp(
  model_system = ms,
  data         = dat,
  control      = ctrl,
  optimizer    = "nlminb",
  parallel     = FALSE,
  verbose      = FALSE
)

fit$convergence   # 0 == strict convergence
fit$loglik

## -----------------------------------------------------------------------------
components_table(fit, digits = 3)

## -----------------------------------------------------------------------------
fscores <- estimate_factorscores_rcpp(
  fit, dat, control = ctrl, parallel = FALSE, verbose = FALSE
)
head(fscores[, c("obs_id", "factor_1", "se_factor_1", "converged")])

# Correlation of estimated scores with the true (simulated) ability
cor(fscores$factor_1, f)

