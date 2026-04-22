# factorana 1.2.0

## Bug fixes

* Stage-2 SE workflow (`define_model_system(..., previous_stage = ...)`
  with `factor_structure = "SE_linear"` or `"SE_quadratic"`) was
  silently fixing `typeprob_*_intercept` and `type_*_loading_*`
  parameters when the previous stage already had `n_types > 1`. The
  `factor_dist_patterns` regex list in `define_model_system()`
  enumerated only `factor_var_*`, `se_*`, `chol_*`, and
  `factor_mean_*` as factor-distribution parameters, so typeprob and
  type-loading slots fell through into `measurement_params` and were
  forced to their Stage 1 values. The C++ initialization correctly
  left those parameters free. The disagreement scrambled the mapping
  between the C++ free-parameter vector and the R-side `free_idx`:
  `evaluate_likelihood_rcpp()` extracted 7 free gradient entries from
  C++ and scattered them into a 5-slot R map, shifting every value
  past position 5. The result was a Hessian whose SE x SE sub-block
  looked permuted (max rel_err ~1.67, previously flagged as the "TEST
  4 known issue"). Added `^typeprob_` and `^type_[0-9]+_loading_` to
  the pattern list; now the two sides agree that those are
  factor-level free parameters. Also matches `TEST 1`'s explicit
  expectation that typeprob and input-factor type loadings stay free
  in Stage 2.

* Type-probability mixture contribution to the Hessian on the
  `d^2 L_mix / d(sigma^2_k)^2` diagonal was missing two terms in
  `FactorModel::CalcLkhd` (section 3b):
  1. The cross term `2 * dpi_t/d(sigma^2_k) * dL_t/d(sigma^2_k)` was
     added only once rather than twice. The off-diagonal code
     correctly added both orderings of the cross term; the `k == l`
     branch kept a "we only add once due to symmetry in sum" comment
     that was wrong. On the diagonal the two orderings coincide so
     the correct behavior is a factor of 2.
  2. The chain-rule second derivative of the type probability through
     the GH scaling of the input factor itself: when `f_k = sigma_k *
     x_q`, the full diagonal second derivative of `pi_t` w.r.t.
     `sigma^2_k` is
     `d^2(pi)/d(f_k)^2 * (df_k/d(sigma^2_k))^2 +
      d(pi)/d(f_k) * d^2(f_k)/d(sigma^2_k)^2`.
     Only the first term was being accumulated; the
     `d(pi)/d(f_k) * d^2(f_k)/d(sigma^2_k)^2` term is now added when
     `k == l`.
  In combination with the free/fixed mapping fix, the Stage-1
  with-types to Stage-2 SE_linear Hessian FD max rel_err dropped from
  1.67 to ~8e-6 at `n_quad = 12`, reaching machine precision at
  `n_quad >= 24`. The previously skipped TEST 4 in
  `test-two-stage-se-types.R` has been re-enabled as a regression
  guard.

# factorana 1.1.8

## Minor changes

* `define_model_component()`: the no-intercept sanity warning now fires
  only when the caller has NOT explicitly opted out of an intercept via
  `intercept = FALSE`. The explicit opt-out is a deliberate choice
  (common in validation and shape-only tests and in ordered-probit
  setups where the intercept is absorbed into the cutpoints); treating
  it as a candidate for the misspecification warning produced noise
  without signal. Accidental omissions (`intercept` left at its default
  of `TRUE`, but no intercept covariate provided) still trigger the
  warning.

* Attempted re-enable of the Stage-1-with-types / Stage-2-SE_linear
  Hessian FD placeholder (`test-two-stage-se-types.R` TEST 4): still
  fails after the v1.1.7 Hessian accumulation fix (max err ~1.7 on
  the SE x SE sub-block, gradient passes). The remaining mismatch is
  a separate issue in how the type-probability model interacts with
  SE_linear under `previous_stage`, not covered by the general
  accumulation fix. Re-skipped with an updated comment.

# factorana 1.1.7

## Bug fixes

* Analytical Hessian now accumulates correctly under equality
  constraints. Previously `FactorModel::CalcLkhd` iterated over
  `freeparlist` for its Hessian contribution loops, which skipped
  equality-tied (derived) parameters. The subsequent
  `ExtractFreeHessian` aggregation had no tied-position values to sum
  into the primaries, so the analytical Hessian at primary x primary
  positions missed the `d^2L / d(primary) d(derived)`, `d^2L /
  d(derived) d(primary)`, and `d^2L / d(derived) d(derived)` terms.
  This produced analytical zeros at positions where finite differences
  (with equality enforced at reinitialisation) reported magnitudes of
  10 to 50, meaning SE estimates under equality constraints were
  biased. Fix: iterate over `gradparlist` (free plus tied) in the
  Hessian loops and symmetrise `full_hessL` before
  `ExtractFreeHessian`. Validated by the now-passing Stage 1 FD tests
  in `test-dynamic-single-factor.R` (linear max_err ~2.5e-6, oprobit
  Hessian max_err ~3.3e-6).

* The C++ name-to-index map for ordered-probit thresholds was looking
  up the field `n_categories` instead of the actual component field
  `num_choices`, so threshold equality constraints were silently
  dropped. Fixed in `rcpp_interface.cpp`. This was the root cause of
  the conv = 1 + "items to replace" warnings observed in the Mental
  Health Trap simulation's Stage 1 oprobit estimation.

## New features

* `define_dynamic_measurement()` for `model_type = "oprobit"` now
  ties only the threshold INCREMENTS (`_thresh_k` for
  `k = 2..K-1`) across periods and leaves `_thresh_1`
  period-specific, mirroring the linear strategy (tie sigmas, free
  intercepts). Patch contributed by an external agent on the MH Trap
  project.

* `define_dynamic_measurement()` now silently strips `"intercept"` /
  `"constant"` from `covariates` when `model_type = "oprobit"` (ordered
  probit absorbs the intercept into the cutpoints). The default
  `covariates = "intercept"` now works for every supported model type.

## New tests

* `test-se-models.R`: `.build_se_type_model` gains an
  `indicator_type` argument (`"linear"` or `"oprobit"`) and two new
  tests exercise single-stage `SE_linear + n_types = 2` FD (gradient
  and Hessian) with ordered-probit indicators.

* `test-two-stage-se-types.R`: TEST 3 adds an oprobit variant of the
  two-stage FD test at Stage 2 (analog of the linear TEST 2). The old
  TEST 3 known-issue placeholder is renumbered to TEST 4.

* `test-dynamic-single-factor.R`:
  * TEST 4: structural test for the oprobit wrapper path (equality
    constraints, no estimation).
  * TEST 5: oprobit Stage 1 converges cleanly with tied thresholds
    (regression guard for the threshold-name-lookup bug).
  * TEST 6-7: Stage 1 FD (gradient + Hessian) for the dynamic-measurement
    wrapper at a non-MLE parameter point, linear and oprobit. These
    tests actively exercise the Hessian-accumulation fix.
  * TEST 8: oprobit two-stage structural parameter recovery. Became
    feasible once the two C++ bugs were fixed; recovery of
    factor_var_1, se_intercept, se_linear_1, se_residual_var is
    within tolerance at n = 2500.

# factorana 1.1.6

## Bug fixes

* `define_dynamic_measurement()` with `model_type = "oprobit"` now ties
  only the threshold INCREMENTS across periods (`_thresh_k` for
  `k = 2..K-1`) and leaves `_thresh_1` period-specific. Previously all
  thresholds were tied, which forced any wave-to-wave shift in the
  latent factor mean into the factor variances and produced Stage 1
  convergence code 1 (false convergence) plus boundary warnings. The
  new behaviour mirrors the linear case: tie the scale (sigmas /
  threshold increments), leave the location (intercepts / thresh_1)
  period-specific. Patch contributed by an external agent working on
  the Mental Health Trap simulation, where Stage 1 convergence went
  from conv = 1 with `factor_var_1 = 4.06` / `factor_var_2 = 0.98` to
  conv = 0 with balanced variances.

* `define_dynamic_measurement()` now silently strips `"intercept"` or
  `"constant"` from `covariates` when `model_type = "oprobit"`. Ordered
  probit absorbs the intercept into the cutpoints, so factorana rejects
  an intercept covariate on oprobit components. The wrapper's default
  `covariates = "intercept"` now works for every supported model type
  without requiring the user to tailor it by model type.

## Tests

* `test-dynamic-single-factor.R` gains a structural test for the
  oprobit wrapper path that verifies: the intercept covariate is
  stripped, equality constraints contain only threshold increments
  `k = 2..K-1` (not `thresh_1`), `thresh_1` appears as a free
  parameter in every period, and `build_dynamic_previous_stage()`
  produces a dummy with no `_intercept` parameters and with the
  anchor-period `thresh_1` carried into every period's slot.
  Estimation-level recovery for oprobit is not asserted: the oprobit
  dynamic model is empirically fragile at moderate n and identification
  is tracked separately.

# factorana 1.1.5

## Bug fix

* `build_dynamic_previous_stage()` now handles oprobit, probit, and
  logit model types. Those types do not have explicit `_intercept`
  parameters (location is absorbed into cutpoints or the link
  function); the previous version would error when looking up a
  missing intercept. Linear behaviour is unchanged.

## Tests

* `test-dynamic-single-factor.R` gains a third test that exercises the
  wrapper plus Stage 2 `SE_linear` with `n_types = 2`. Types shift the
  period-2 factor mean via `se_intercept_type_2`. Recovery of the
  well-identified structural parameters (factor_var_1, se_linear_1,
  se_intercept_type_2, se_residual_var) is verified on a simulated
  DGP, with a more generous tolerance on `se_intercept` (which trades
  off against `typeprob_2_intercept` and `type_2_loading_1` when
  measurement information density is low).

# factorana 1.1.4

## New features

* New helper functions `define_dynamic_measurement()` and
  `build_dynamic_previous_stage()` encapsulate the standard workflow for
  estimating an SE_linear or SE_quadratic structural model on a single
  latent construct observed at two or more time points:
  - `define_dynamic_measurement()` builds the Stage 1 measurement model
    (a k-factor independent system with loadings and residual sigmas
    tied across periods via `equality_constraints`, measurement
    intercepts left period-specific).
  - `build_dynamic_previous_stage()` constructs a Stage 2
    `previous_stage` object that plugs the anchor-period (wave 1 by
    default) intercepts into every factor slot. This anchors the
    measurement level under the factor-identification convention
    `E[f_k] = 0` and lets the observed period-to-period mean shift in
    Y identify the structural intercept (alpha) in Stage 2.
  - Supports `model_type` of "linear", "oprobit", "probit", and
    "logit" with appropriate tying of thresholds (oprobit) or sigmas
    (linear).
* Refactored `tests/testthat/test-dynamic-single-factor.R` to use the
  wrapper; recovery results are unchanged.
* Refactored `vignettes/dynamic_structural.Rmd` to use the wrapper.

# factorana 1.1.3

## New tests and documentation

* New `tests/testthat/test-dynamic-single-factor.R` covers the standard
  workflow for estimating an SE_linear dynamic structural equation on a
  single latent construct measured at two time points:
  - Stage 1 fits a 2-factor independent measurement model with
    `equality_constraints` tying factor loadings and residual sigmas
    across periods but leaving measurement intercepts period-specific.
  - A dummy 2-factor `previous_stage` object is built that carries the
    wave-1 intercepts into both factor slots (discarding the wave-2
    intercepts, which absorb the structural-equation mean shift).
  - Stage 2 fits `SE_linear` and recovers the structural intercept
    `se_intercept` (alpha), slope `se_linear_1` (beta), residual
    variance `se_residual_var`, and input factor variance
    `factor_var_1`.
* New vignette `vignettes/dynamic_structural.Rmd` walks through the
  same workflow with executable code and explains the motivation for
  using wave-1 intercepts in Stage 2.
* Removed the parameter-recovery test from
  `test-two-stage-se-types.R` (its Stage-1-no-types ->
  Stage-2-with-types setup is not a canonical workflow); the shape
  test, FD gradient/Hessian test, and the skipped
  Stage-1-with-types known-issue placeholder remain in place.

# factorana 1.1.2

## Bug fixes

* Two-stage estimation with `previous_stage + SE_linear/SE_quadratic` and
  `n_types > 1` now correctly builds the Stage 2 parameter vector. Prior
  versions omitted the `typeprob_*` and `type_*_loading_*` slots, causing
  either a crash in `setup_parameter_constraints()` or silently mis-fixed
  parameters. The measurement-parameter filter in
  `initialize_parameters()` was also tightened so that factor-level type
  parameters from Stage 1 are no longer duplicated into the measurement
  block. Discovered during analysis of a structural model where types are
  introduced at Stage 2.

## Tests

* New `tests/testthat/test-two-stage-se-types.R` adds:
  - a shape test that verifies Stage 2 SE_linear + `n_types = 2` produces
    the expected parameter vector aligned with `build_parameter_metadata()`,
  - a finite-difference gradient and Hessian check at the DGP parameters,
  - a structural parameter recovery test (`se_linear_1`,
    `se_intercept_type_2`, `se_residual_var`, `factor_var_1`) with init
    `se_intercept = -0.5` (Stage 1 absorbs E[f2] into the measurement
    intercepts, so the MLE `se_intercept` is negative even if the DGP
    constant is 0), and
  - a skipped placeholder documenting a known Hessian-FD mismatch in the
    Stage-1-with-types -> Stage-2-SE_linear variant (not the common
    workflow; tracked for a future fix).

# factorana 1.1.1

## CRAN resubmission fixes
* Added method references (with DOIs) to the `DESCRIPTION` field: Heckman,
  Humphries & Veramendi (2016, 2018) and Humphries, Joensen & Veramendi (2024).
* Added `\value` (via `@return`) to all exported functions that were missing it,
  including `as_kv`, `estimate_and_write`, `write_model_config_csv`, the
  adaptive-quadrature and observation-weight setters, and every `print` method.
* Replaced `\dontrun{}` with runnable examples (`fix_coefficient`,
  `fix_type_intercepts`) or `\donttest{}` blocks (`results_table`,
  `results_to_latex`, `components_table`, `estimate_factorscores_rcpp`,
  `cleanup_parallel_workers`).
* `estimate_and_write()` and `write_model_config_csv()` no longer write to a
  default path; `results_dir` / `file` is now required (use `tempdir()` in
  examples and tests).
* Replaced the non-executing `introduction` vignette with two executable
  vignettes: `measurement_system` (two-factor CFA) and `roy_model` (sector
  choice with a latent ability factor).

## Bug fixes
* Fixed systematic-suite test that constructed an ordered-probit component
  with `intercept` in its covariates; this configuration is rejected (the
  intercept is absorbed into the cut points).

# factorana 1.0.2

## Improvements
* CRAN compliance fixes
* Added introductory vignette
* Improved test coverage for SE models, equality constraints, observation weights
* Documentation improvements: fixed adaptive integration formula in README
* Updated SE_linear example to use larger sample size for reliable convergence

# factorana 1.0.1

## Bug Fixes
* Fix binary logit initialization and add dedicated test

# factorana 1.0.0

## New Features
* Structural equation models (SE_linear, SE_quadratic) for causal factor relationships
* Mixture of normals factor distribution (n_mixtures = 1, 2, or 3)
* Equality constraints for measurement invariance via `equality_constraints` parameter
* Component-level type control via `use_types` parameter
* Observation weights for survey weights or importance sampling
* Checkpointing for long-running estimations via `checkpoint_file` parameter
* Exploded multinomial logit for ranked choice models
* Exploded nested logit with `exclude_chosen = FALSE`
* Rank-share corrections via `rankshare_var` parameter

## Core Features
* Multi-factor models with flexible loading normalization
* Linear, probit, ordered probit, and multinomial logit model components
* Analytical gradients and Hessians for fast convergence
* Parallel estimation via doParallel for large datasets
* Multi-stage (sequential) estimation with fixed early-stage parameters
* Adaptive integration for efficient two-stage estimation
* Factor interaction terms (quadratic and cross-product) via `factor_spec`
* Correlated two-factor models via `factor_structure = "correlation"`
