## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## -----------------------------------------------------------------------------
library(factorana)

set.seed(41)
n <- 1500

# Structural parameters (what we want to recover)
true_var_f1   <- 1.0
true_alpha    <- 0.4
true_beta     <- 0.6
true_sigma_e2 <- 0.5

# Shared measurement parameters
item_int   <- c(1.5, 1.0, 0.8)
item_load  <- c(1.0, 0.9, 1.1)   # first loading fixed to 1 in the model
item_sigma <- c(0.7, 0.75, 0.65)

f1  <- rnorm(n, 0, sqrt(true_var_f1))
eps <- rnorm(n, 0, sqrt(true_sigma_e2))
f2  <- true_alpha + true_beta * f1 + eps

gen_Y <- function(f, i) {
  item_int[i] + item_load[i] * f + rnorm(length(f), 0, item_sigma[i])
}

dat <- data.frame(
  intercept = 1,
  eval      = 1L,
  Y_t1_m1 = gen_Y(f1, 1), Y_t1_m2 = gen_Y(f1, 2), Y_t1_m3 = gen_Y(f1, 3),
  Y_t2_m1 = gen_Y(f2, 1), Y_t2_m2 = gen_Y(f2, 2), Y_t2_m3 = gen_Y(f2, 3)
)

## -----------------------------------------------------------------------------
dyn <- define_dynamic_measurement(
  data                 = dat,
  items                = c("m1", "m2", "m3"),
  period_prefixes      = c("Y_t1_", "Y_t2_"),
  model_type           = "linear",
  evaluation_indicator = "eval"
)
# The wrapper generates 5 equality constraints: 2 for loadings (items
# m2, m3; item m1's loading is fixed to 1 on its factor slot) and 3
# for sigmas.
length(dyn$equality_constraints)

## -----------------------------------------------------------------------------
ctrl <- define_estimation_control(n_quad_points = 8, num_cores = 1)
result_stage1 <- estimate_model_rcpp(
  model_system = dyn$model_system,
  data         = dat,
  control      = ctrl,
  optimizer    = "nlminb",
  parallel     = FALSE,
  verbose      = FALSE
)
result_stage1$convergence

## -----------------------------------------------------------------------------
est <- result_stage1$estimates
tab <- data.frame(
  m        = 1:3,
  DGP_tau  = item_int,
  wave_1   = round(c(est["Y_t1_m1_intercept"],
                     est["Y_t1_m2_intercept"],
                     est["Y_t1_m3_intercept"]), 3),
  wave_2   = round(c(est["Y_t2_m1_intercept"],
                     est["Y_t2_m2_intercept"],
                     est["Y_t2_m3_intercept"]), 3)
)
knitr::kable(tab, row.names = FALSE)

## -----------------------------------------------------------------------------
dummy <- build_dynamic_previous_stage(
  dyn           = dyn,
  stage1_result = result_stage1,
  data          = dat,
  anchor_period = 1L
)

fm_stage2 <- define_factor_model(
  n_factors        = 2,
  n_types          = 1,
  factor_structure = "SE_linear"
)
ms_stage2 <- define_model_system(
  components     = list(),       # measurement components prepended from previous_stage
  factor         = fm_stage2,
  previous_stage = dummy
)

init_s2 <- initialize_parameters(ms_stage2, dat, verbose = FALSE)
init_s2$init_params["factor_var_1"]    <- unname(dummy$estimates["factor_var_1"])
init_s2$init_params["se_intercept"]    <- 0.0
init_s2$init_params["se_linear_1"]     <- 0.5
init_s2$init_params["se_residual_var"] <- 0.5

result_stage2 <- estimate_model_rcpp(
  model_system = ms_stage2,
  data         = dat,
  init_params  = init_s2$init_params,
  control      = ctrl,
  optimizer    = "nlminb",
  parallel     = FALSE,
  verbose      = FALSE
)
result_stage2$convergence

## -----------------------------------------------------------------------------
est <- result_stage2$estimates
se  <- result_stage2$std_errors
ps  <- c("factor_var_1", "se_intercept", "se_linear_1", "se_residual_var")

tab <- data.frame(
  param = ps,
  true  = c(true_var_f1, true_alpha, true_beta, true_sigma_e2),
  est   = round(unname(est[ps]), 4),
  se    = round(unname(se[ps]),  4)
)
tab$z <- round((tab$est - tab$true) / tab$se, 2)
knitr::kable(tab, row.names = FALSE)

