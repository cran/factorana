#include "FactorModel.h"
#include "distributions.h"
#include "gauss_hermite.h"
#include <R.h>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <set>

FactorModel::FactorModel(int n_obs, int n_var, int n_fac, int n_typ,
                         int n_mix, bool correlated, int n_quad)
    : nobs(n_obs), nvar(n_var), nfac(n_fac), ntyp(n_typ),
      nmix(n_mix), fac_corr(correlated),
      factor_structure(correlated ? FactorStructure::CORRELATION : FactorStructure::INDEPENDENT),
      n_factors_for_mixture(n_fac), n_variance_per_mixture(n_fac),
      n_input_factors(0), n_outcome_factors(0), se_param_start(-1), nse_param(0),
      nse_type_intercepts(0),
      nquad_points(n_quad), use_weights(false), use_adaptive(false), adapt_threshold(0.3),
      use_factor_mean_covariates(false), n_factor_mean_covariates(0), n_factors_with_mean(0),
      factor_mean_param_start(-1),
      use_se_covariates(false), n_se_covariates(0), se_covariate_param_start(-1)
{
    data.resize(nobs * nvar);

    // Initialize with factor parameters
    // For nmix mixtures: nmix * nfac factor variances
    nparam = nfac * nmix;  // Start with factor variances

    // Add correlation parameters if correlated factors
    // For 2-factor correlated model: 1 correlation parameter
    if (fac_corr && nfac == 2) {
        n_variance_per_mixture = nfac + 1;
        nparam = nfac + 1;  // Reset - correlation not supported with mixtures
    }

    // Mixture parameters: means and log-weights for non-reference mixtures
    if (nmix > 1 && !fac_corr) {
        nparam += (nmix - 1) * nfac;  // mixture means
        nparam += (nmix - 1);          // mixture log-weights
    }

    // Compute type model parameters
    // Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    // For ntyp > 1: (ntyp - 1) intercepts + (ntyp - 1) * nfac loading parameters
    if (ntyp > 1) {
        n_typeprob_intercepts = ntyp - 1;
        ntyp_param = n_typeprob_intercepts + (ntyp - 1) * nfac;
        type_param_start = nparam;  // Type params come after factor params
        nparam += ntyp_param;
    } else {
        n_typeprob_intercepts = 0;
        ntyp_param = 0;
        type_param_start = -1;
    }

    nparam_free = nparam;  // All factor parameters are free by default
}

// Constructor with factor_structure support
FactorModel::FactorModel(int n_obs, int n_var, int n_fac, int n_typ,
                         int n_mix, FactorStructure fac_struct, int n_quad)
    : nobs(n_obs), nvar(n_var), nfac(n_fac), ntyp(n_typ),
      nmix(n_mix), fac_corr(fac_struct == FactorStructure::CORRELATION),
      factor_structure(fac_struct),
      n_factors_for_mixture(n_fac), n_variance_per_mixture(n_fac),
      n_input_factors(0), n_outcome_factors(0), se_param_start(-1), nse_param(0),
      nse_type_intercepts(0),
      nquad_points(n_quad), use_weights(false), use_adaptive(false), adapt_threshold(0.3),
      use_factor_mean_covariates(false), n_factor_mean_covariates(0), n_factors_with_mean(0),
      factor_mean_param_start(-1),
      use_se_covariates(false), n_se_covariates(0), se_covariate_param_start(-1)
{
    data.resize(nobs * nvar);

    // Initialize factor parameters based on structure
    if (factor_structure == FactorStructure::SE_LINEAR) {
        // For SE models: input factors have variance parameters, outcome factors don't
        // Integration is over (f_1, ..., f_{k-1}, epsilon_k)
        // Mixtures only apply to input factors
        n_input_factors = nfac - 1;
        n_outcome_factors = 1;
        n_factors_for_mixture = n_input_factors;
        n_variance_per_mixture = n_input_factors;

        // Factor variance parameters: only for input factors (per mixture)
        nparam = n_input_factors * nmix;

        // Mixture parameters: means and log-weights for non-reference mixtures
        if (nmix > 1) {
            nparam += (nmix - 1) * n_input_factors;  // mixture means
            nparam += (nmix - 1);                     // mixture log-weights
        }

        // SE parameters: intercept + n_input_factors linear coefficients + residual variance
        // When ntyp > 1, also add (ntyp-1) type-specific SE intercepts (between linear and residual var)
        se_param_start = nparam;
        nse_type_intercepts = (ntyp > 1) ? (ntyp - 1) : 0;
        nse_param = 1 + n_input_factors + nse_type_intercepts + 1;  // intercept + linear + type intercepts + residual var
        nparam += nse_param;

    } else if (factor_structure == FactorStructure::SE_QUADRATIC) {
        // SE_QUADRATIC: f_k = alpha + alpha_1*f_1 + alpha_q1*f_1^2 + ... + epsilon
        // Integration is over (f_1, ..., f_{k-1}, epsilon_k)
        // Mixtures only apply to input factors
        n_input_factors = nfac - 1;
        n_outcome_factors = 1;
        n_factors_for_mixture = n_input_factors;
        n_variance_per_mixture = n_input_factors;

        // Factor variance parameters: only for input factors (per mixture)
        nparam = n_input_factors * nmix;

        // Mixture parameters: means and log-weights for non-reference mixtures
        if (nmix > 1) {
            nparam += (nmix - 1) * n_input_factors;  // mixture means
            nparam += (nmix - 1);                     // mixture log-weights
        }

        // SE parameters: intercept + linear coefs + quadratic coefs + residual variance
        // When ntyp > 1, also add (ntyp-1) type-specific SE intercepts (between quadratic and residual var)
        se_param_start = nparam;
        nse_type_intercepts = (ntyp > 1) ? (ntyp - 1) : 0;
        nse_param = 1 + n_input_factors + n_input_factors + nse_type_intercepts + 1;  // intercept + linear + quadratic + type intercepts + residual var
        nparam += nse_param;

    } else if (factor_structure == FactorStructure::CORRELATION && nfac == 2) {
        // Correlated factors: nfac variances + 1 correlation
        // Note: correlation structure not supported with mixtures
        n_factors_for_mixture = nfac;
        n_variance_per_mixture = nfac + 1;  // variances + correlation
        nparam = nfac + 1;

    } else {
        // Independent factors: nfac variances (per mixture)
        n_factors_for_mixture = nfac;
        n_variance_per_mixture = nfac;
        nparam = nfac * nmix;

        // Mixture parameters: means and log-weights for non-reference mixtures
        if (nmix > 1) {
            nparam += (nmix - 1) * nfac;  // mixture means
            nparam += (nmix - 1);          // mixture log-weights
        }
    }

    // Compute type model parameters
    // Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    // For ntyp > 1: (ntyp - 1) intercepts + (ntyp - 1) * nfac loading parameters
    if (ntyp > 1) {
        n_typeprob_intercepts = ntyp - 1;
        ntyp_param = n_typeprob_intercepts + (ntyp - 1) * nfac;
        type_param_start = nparam;
        nparam += ntyp_param;
    } else {
        n_typeprob_intercepts = 0;
        ntyp_param = 0;
        type_param_start = -1;
    }

    nparam_free = nparam;
}

void FactorModel::AddModel(std::shared_ptr<Model> model, int nparams)
{
    models.push_back(model);
    param_model_start.push_back(nparam);
    param_model_count.push_back(nparams);
    nparam += nparams;
}

void FactorModel::SetData(const std::vector<double>& dat)
{
    if (dat.size() != nobs * nvar) {
        throw std::runtime_error("Data size mismatch");
    }
    data = dat;
}

void FactorModel::SetData(const Eigen::MatrixXd& dat)
{
    if (dat.rows() != nobs || dat.cols() != nvar) {
        throw std::runtime_error("Data dimension mismatch");
    }
    data.resize(nobs * nvar);
    for (int i = 0; i < nobs; i++) {
        for (int j = 0; j < nvar; j++) {
            data[i * nvar + j] = dat(i, j);
        }
    }
}

void FactorModel::SetQuadrature(const std::vector<double>& nodes,
                                const std::vector<double>& weights)
{
    quad_nodes = nodes;
    quad_weights = weights;
}

void FactorModel::SetObservationWeights(const std::vector<double>& weights)
{
    if (weights.size() != nobs) {
        throw std::runtime_error("Observation weights size mismatch: expected " +
                                 std::to_string(nobs) + " but got " +
                                 std::to_string(weights.size()));
    }
    obs_weights = weights;
    use_weights = true;
}

void FactorModel::SetAdaptiveQuadrature(const std::vector<std::vector<double>>& factor_scores,
                                        const std::vector<std::vector<double>>& factor_ses,
                                        const std::vector<double>& factor_vars,
                                        double threshold,
                                        int max_quad)
{
    // Validate inputs
    if (factor_scores.size() != nobs) {
        throw std::runtime_error("factor_scores size mismatch: expected " +
                                 std::to_string(nobs) + " observations");
    }
    if (factor_ses.size() != nobs) {
        throw std::runtime_error("factor_ses size mismatch: expected " +
                                 std::to_string(nobs) + " observations");
    }
    if (factor_vars.size() != nfac) {
        throw std::runtime_error("factor_vars size mismatch: expected " +
                                 std::to_string(nfac) + " factors");
    }

    // Store settings
    adapt_threshold = threshold;
    adapt_factor_var = factor_vars;

    // Resize storage
    obs_nquad.resize(nobs);
    obs_fac_center.resize(nobs);
    obs_fac_se.resize(nobs);
    obs_weights.resize(nobs, 1.0);

    // Collect all unique nquad values needed
    std::set<int> needed_nquad;

    // Compute per-observation quadrature settings
    for (int iobs = 0; iobs < nobs; iobs++) {
        if (factor_scores[iobs].size() != nfac || factor_ses[iobs].size() != nfac) {
            throw std::runtime_error("Factor score/SE dimensions mismatch for observation " +
                                     std::to_string(iobs));
        }

        obs_nquad[iobs].resize(nfac);
        obs_fac_center[iobs].resize(nfac);
        obs_fac_se[iobs].resize(nfac);

        for (int ifac = 0; ifac < nfac; ifac++) {
            double f_score = factor_scores[iobs][ifac];
            double f_se = factor_ses[iobs][ifac];
            double f_var = factor_vars[ifac];
            double f_sd = std::sqrt(f_var);

            // Handle NaN (from R NA) factor scores or SEs: fall back to full
            // standard quadrature. This occurs for zero-observation individuals
            // where factor scores are undefined (R returns NA_real_).
            if (std::isnan(f_score) || std::isnan(f_se)) {
                obs_nquad[iobs][ifac] = max_quad;
                obs_fac_center[iobs][ifac] = 0.0;  // Center at prior mean
                obs_fac_se[iobs][ifac] = f_sd;     // Use full SD as spread
                needed_nquad.insert(max_quad);
                continue;
            }

            // Handle SE=0 or very small SE: use standard quadrature
            // SE=0 would cause jacobian=0 in IS correction, making likelihood 0
            if (f_se < 1e-10) {
                obs_nquad[iobs][ifac] = max_quad;
                obs_fac_center[iobs][ifac] = 0.0;  // Center at prior mean
                obs_fac_se[iobs][ifac] = f_sd;     // Use full SD as spread
                needed_nquad.insert(max_quad);
                continue;  // Skip rest of this factor's processing
            }

            // Store factor score center and SE
            obs_fac_center[iobs][ifac] = f_score;
            obs_fac_se[iobs][ifac] = f_se;

            // Determine number of quadrature points based on SE
            // Formula from legacy code: nquad = 1 + 2 * floor(se / sd / threshold)
            // When SE is small relative to SD, use 1 point (just the factor score)
            // When SE is large relative to SD, use more points
            double ratio = f_se / f_sd / threshold;
            int nq = 1 + 2 * static_cast<int>(std::floor(ratio));

            // Clamp to reasonable range
            if (nq < 1) nq = 1;
            if (nq > max_quad) nq = max_quad;

            // If SE is very large (> sqrt(factor_var)), use default quadrature centered at 0
            // and use full SD as the spread (essentially standard quadrature)
            if (f_se > f_sd) {
                nq = max_quad;
                obs_fac_center[iobs][ifac] = 0.0;  // Center at prior mean
                obs_fac_se[iobs][ifac] = f_sd;     // Use full SD as spread
            }

            obs_nquad[iobs][ifac] = nq;
            needed_nquad.insert(nq);
        }

        // Note: Importance sampling weights are computed per integration point in CalcLkhd,
        // not per observation, because they depend on the quadrature node values.
        obs_weights[iobs] = 1.0;  // Default weight
    }

    // Pre-compute GH quadrature nodes/weights for all needed sizes
    adapt_nodes.clear();
    adapt_weights.clear();

    double sqrt2 = std::sqrt(2.0);
    double sqrt_pi = std::sqrt(M_PI);

    for (int nq : needed_nquad) {
        std::vector<double> nodes, weights;
        calcgausshermitequadrature(nq, nodes, weights);

        // Scale for standard normal: x_scaled = sqrt(2) * x, w_scaled = w / sqrt(pi)
        for (int i = 0; i < nq; i++) {
            nodes[i] *= sqrt2;
            weights[i] /= sqrt_pi;
        }

        adapt_nodes[nq] = nodes;
        adapt_weights[nq] = weights;
    }

    use_weights = true;
    use_adaptive = true;
}

void FactorModel::DisableAdaptiveQuadrature()
{
    use_adaptive = false;
    use_weights = false;  // Also reset weights since they may have been set by adaptive mode
    obs_nquad.clear();
    obs_fac_center.clear();
    obs_fac_se.clear();
    obs_weights.clear();
    adapt_factor_var.clear();
    adapt_nodes.clear();
    adapt_weights.clear();
}

void FactorModel::SetFactorMeanCovariates(const std::vector<std::vector<double>>& covariate_data)
{
    // Validate input dimensions
    if (covariate_data.size() != nobs) {
        throw std::runtime_error("Factor mean covariate data has wrong number of observations: expected " +
                                 std::to_string(nobs) + " but got " + std::to_string(covariate_data.size()));
    }

    if (covariate_data.empty() || covariate_data[0].empty()) {
        use_factor_mean_covariates = false;
        return;
    }

    n_factor_mean_covariates = static_cast<int>(covariate_data[0].size());

    // Determine which factors get mean covariates
    // For SE models, only input factors get covariates (outcome factor mean is from SE)
    if (factor_structure == FactorStructure::SE_LINEAR ||
        factor_structure == FactorStructure::SE_QUADRATIC) {
        n_factors_with_mean = n_input_factors;
    } else {
        n_factors_with_mean = nfac;
    }

    // Compute parameter index for factor mean coefficients
    // They come after variance parameters and any SE/correlation parameters
    // The order is: factor_var_1..k, [se_params or corr_params], factor_mean_1_cov1, factor_mean_1_cov2, ...
    if (factor_structure == FactorStructure::SE_LINEAR ||
        factor_structure == FactorStructure::SE_QUADRATIC) {
        // After input factor variances and SE params
        factor_mean_param_start = n_input_factors + nse_param;
    } else if (factor_structure == FactorStructure::CORRELATION && nfac == 2) {
        // After factor variances and correlation
        factor_mean_param_start = nfac + 1;
    } else {
        // After factor variances
        factor_mean_param_start = nfac;
    }

    // Store the pre-processed (demeaned) covariate data
    // Demeaning and variance checking is done in R before calling this
    factor_mean_covariate_data = covariate_data;

    // Update total parameter count to include factor mean coefficients
    // n_factors_with_mean * n_factor_mean_covariates additional parameters
    int n_mean_params = n_factors_with_mean * n_factor_mean_covariates;
    nparam += n_mean_params;

    use_factor_mean_covariates = true;
}

void FactorModel::SetSECovariates(const std::vector<std::vector<double>>& covariate_data)
{
    // Only valid for SE models
    if (factor_structure != FactorStructure::SE_LINEAR &&
        factor_structure != FactorStructure::SE_QUADRATIC) {
        throw std::runtime_error("SE covariates are only valid for SE_LINEAR and SE_QUADRATIC factor structures");
    }

    // Validate input dimensions
    if (covariate_data.size() != nobs) {
        throw std::runtime_error("SE covariate data has wrong number of observations: expected " +
                                 std::to_string(nobs) + " but got " + std::to_string(covariate_data.size()));
    }

    if (covariate_data.empty() || covariate_data[0].empty()) {
        use_se_covariates = false;
        return;
    }

    n_se_covariates = static_cast<int>(covariate_data[0].size());

    // SE covariate parameters come after all other factor-level parameters
    // Order: [input_factor_vars | se_params | factor_mean_params | se_covariate_params | model_params]
    se_covariate_param_start = nparam;  // Current nparam is where SE cov params start

    // Store the pre-processed (demeaned) covariate data
    se_covariate_data = covariate_data;

    // Update total parameter count
    nparam += n_se_covariates;

    use_se_covariates = true;
}

void FactorModel::SetParameterConstraints(const std::vector<bool>& fixed)
{
    if (fixed.size() != nparam) {
        throw std::runtime_error("Parameter constraint size mismatch");
    }
    param_fixed = fixed;
    nparam_free = 0;
    freeparlist.clear();
    freeparlist.reserve(nparam);
    for (int i = 0; i < nparam; i++) {
        if (!fixed[i]) {
            freeparlist.push_back(i);
            nparam_free++;
        }
    }

    // Initialize gradparlist to freeparlist (will be updated if equality constraints are set)
    gradparlist = freeparlist;

    // Compute model-local free parameter indices for Hessian optimization
    // This allows Model::Eval() to skip Hessian computation for fixed parameters
    model_free_indices.clear();
    model_free_indices.resize(models.size());
    for (size_t imod = 0; imod < models.size(); imod++) {
        int start = param_model_start[imod];
        int count = param_model_count[imod];
        model_free_indices[imod].clear();
        model_free_indices[imod].reserve(count);
        for (int j = 0; j < count; j++) {
            int global_idx = start + j;
            if (!fixed[global_idx]) {
                model_free_indices[imod].push_back(j);  // Store model-local index
            }
        }
    }

    // Initialize parameter vector with default values
    param.resize(nparam, 0.0);

    // Initialize factor variances to 1.0 by default
    for (int ifac = 0; ifac < nfac; ifac++) {
        param[ifac] = 1.0;
    }
}

void FactorModel::SetParameterConstraints(const std::vector<bool>& fixed,
                                          const std::vector<double>& fixed_values)
{
    if (fixed.size() != nparam) {
        throw std::runtime_error("Parameter constraint size mismatch");
    }
    if (fixed_values.size() != nparam) {
        throw std::runtime_error("Fixed values size mismatch");
    }

    param_fixed = fixed;
    nparam_free = 0;
    freeparlist.clear();
    freeparlist.reserve(nparam);
    for (int i = 0; i < nparam; i++) {
        if (!fixed[i]) {
            freeparlist.push_back(i);
            nparam_free++;
        }
    }

    // Initialize gradparlist to freeparlist (will be updated if equality constraints are set)
    gradparlist = freeparlist;

    // Compute model-local free parameter indices for Hessian optimization
    // This allows Model::Eval() to skip Hessian computation for fixed parameters
    model_free_indices.clear();
    model_free_indices.resize(models.size());
    for (size_t imod = 0; imod < models.size(); imod++) {
        int start = param_model_start[imod];
        int count = param_model_count[imod];
        model_free_indices[imod].clear();
        model_free_indices[imod].reserve(count);
        for (int j = 0; j < count; j++) {
            int global_idx = start + j;
            if (!fixed[global_idx]) {
                model_free_indices[imod].push_back(j);  // Store model-local index
            }
        }
    }

    // Initialize parameter vector with provided values
    // For free parameters, these are defaults (will be overwritten by optimizer)
    // For fixed parameters, these values are used permanently
    param = fixed_values;
}

void FactorModel::SetEqualityConstraints(const std::vector<int>& equality_map)
{
    if (equality_map.size() != nparam) {
        throw std::runtime_error("Equality mapping size mismatch: expected " +
                                 std::to_string(nparam) + " but got " +
                                 std::to_string(equality_map.size()));
    }
    equality_mapping = equality_map;

    // Build gradparlist: params that need gradient computation (free + tied)
    // Tied params need gradients so they can be aggregated to their primaries
    gradparlist.clear();
    std::set<int> grad_params_set;

    // First add all free params
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) {
            grad_params_set.insert(i);
        }
    }

    // Then add all tied params (they're fixed but need gradients)
    for (int i = 0; i < nparam; i++) {
        if (equality_mapping[i] >= 0) {
            grad_params_set.insert(i);
        }
    }

    // Convert to sorted vector
    gradparlist.assign(grad_params_set.begin(), grad_params_set.end());
}

void FactorModel::MapFreeToFull(const std::vector<double>& free_params)
{
    if (free_params.size() != nparam_free) {
        throw std::runtime_error("Free parameter size mismatch");
    }

    int ifree = 0;
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) {
            param[i] = free_params[ifree++];
        }
        // Fixed parameters keep their current values in param vector
    }

    // For equality constraints, set tied parameter values to their primary's value
    if (!equality_mapping.empty()) {
        for (int i = 0; i < nparam; i++) {
            if (equality_mapping[i] >= 0) {
                param[i] = param[equality_mapping[i]];
            }
        }
    }
}

void FactorModel::ExtractFreeGradient(const std::vector<double>& full_grad,
                                      std::vector<double>& free_grad)
{
    // If equality constraints exist, first aggregate tied gradients to primaries
    std::vector<double> aggregated_grad = full_grad;

    if (!equality_mapping.empty()) {
        for (int i = 0; i < nparam; i++) {
            if (equality_mapping[i] >= 0) {
                // Parameter i is tied to primary parameter equality_mapping[i]
                // Add its gradient to the primary's gradient
                aggregated_grad[equality_mapping[i]] += full_grad[i];
            }
        }
    }

    // Now extract only free parameters
    free_grad.resize(nparam_free);
    int ifree = 0;
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) {
            free_grad[ifree++] = aggregated_grad[i];
        }
    }
}

void FactorModel::ExtractFreeHessian(const std::vector<double>& full_hess,
                                     std::vector<double>& free_hess)
{
    // Extract submatrix corresponding to free parameters
    // Input: full_hess is stored as nparam x nparam matrix (row-major, full matrix)
    // Output: free_hess is upper triangle of nparam_free x nparam_free matrix

    // If equality constraints exist, first aggregate tied Hessian elements
    // For tied param i -> primary p and tied param j -> primary q:
    //   H[p,q] += H[i,j]
    std::vector<double> aggregated_hess = full_hess;

    if (!equality_mapping.empty()) {
        // Create mapping from full param index to effective index (considering ties)
        std::vector<int> effective_idx(nparam);
        for (int i = 0; i < nparam; i++) {
            effective_idx[i] = (equality_mapping[i] >= 0) ? equality_mapping[i] : i;
        }

        // Aggregate Hessian elements to their effective (primary) positions
        // Need a clean aggregation matrix first
        std::vector<double> temp_hess(nparam * nparam, 0.0);

        for (int i = 0; i < nparam; i++) {
            int eff_i = effective_idx[i];
            for (int j = 0; j < nparam; j++) {
                int eff_j = effective_idx[j];
                // Original Hessian value at (i,j)
                double h_ij = full_hess[i * nparam + j];
                // Add to effective position
                temp_hess[eff_i * nparam + eff_j] += h_ij;
            }
        }
        aggregated_hess = temp_hess;
    }

    // Now extract only free parameters
    free_hess.clear();
    free_hess.reserve(nparam_free * (nparam_free + 1) / 2);

    std::vector<int> free_indices;
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) free_indices.push_back(i);
    }

    for (int i = 0; i < nparam_free; i++) {
        for (int j = i; j < nparam_free; j++) {
            int full_i = free_indices[i];
            int full_j = free_indices[j];
            int full_idx = full_i * nparam + full_j;
            free_hess.push_back(aggregated_hess[full_idx]);
        }
    }
}

double FactorModel::CalcLogLikelihood(const std::vector<double>& free_params)
{
    double logLkhd;
    std::vector<double> gradL, hessL;
    CalcLkhd(free_params, logLkhd, gradL, hessL, 1);  // flag=1: likelihood only
    return logLkhd;
}

void FactorModel::CalcLkhd(const std::vector<double>& free_params,
                           double& logLkhd,
                           std::vector<double>& gradL,
                           std::vector<double>& hessL,
                           int iflag)
{
    // ===== STEP 1: Map free parameters to full parameter vector =====
    MapFreeToFull(free_params);

    // ===== STEP 2: Extract factor parameters from param vector =====
    // Parameter organization depends on factor_structure and nmix:
    //   For nmix > 1: [variances (nmix*nfac_mix) | means ((nmix-1)*nfac_mix) | logweights (nmix-1) | structure_params | model_params]
    //   INDEPENDENT: nfac_mix = nfac
    //   SE_LINEAR/SE_QUADRATIC: nfac_mix = nfac-1 (input factors only), + SE params
    //   CORRELATION: nfac_mix = nfac (not supported with mixtures)

    // Factor correlation (for 2-factor correlated models)
    double rho = 0.0;
    double sqrt_1_minus_rho2 = 1.0;

    // SE parameters (for structural equation models)
    double se_intercept = 0.0;
    std::vector<double> se_linear_coef(n_input_factors, 0.0);
    std::vector<double> se_quadratic_coef(n_input_factors, 0.0);  // For SE_QUADRATIC
    std::vector<double> se_covariate_coef;  // Coefficients for SE covariates
    // Type-specific SE intercepts: se_type_offset[ityp] is the additional intercept
    // applied to the outcome factor when evaluating under type ityp. Index 0 (reference
    // type) is always 0; indices 1..ntyp-1 are free parameters when ntyp > 1.
    std::vector<double> se_type_offset(std::max(ntyp, 1), 0.0);
    double se_residual_var = 1.0;
    double sigma_eps = 1.0;  // sqrt of residual variance

    // Number of integration dimensions (for SE models: input factors + residuals)
    int n_integ_dims = nfac;  // Default: one dimension per factor

    // ===== STEP 2a: Compute mixture weights and means =====
    std::vector<double> mix_weights(nmix);
    std::vector<std::vector<double>> mix_means(nmix);
    std::vector<std::vector<double>> mix_sigmas(nmix);  // sqrt of variances per mixture

    ComputeMixtureWeights(mix_weights);
    ComputeMixtureMeans(mix_weights, mix_means);

    // ===== STEP 2b: Extract variances for each mixture =====
    for (int imix = 0; imix < nmix; imix++) {
        mix_sigmas[imix].resize(n_factors_for_mixture);
        for (int k = 0; k < n_factors_for_mixture; k++) {
            double var_k = std::fabs(param[GetFactorVarianceIndex(imix, k)]);
            mix_sigmas[imix][k] = std::sqrt(var_k);
        }
    }

    // For single mixture (backward compatibility), use first mixture's parameters
    std::vector<double> factor_var(nfac, 1.0);  // Default variances for all factors
    std::vector<double> factor_mean(nfac, 0.0);  // Default means for all factors

    if (factor_structure == FactorStructure::SE_LINEAR ||
        factor_structure == FactorStructure::SE_QUADRATIC) {
        // SE models: only input factors have variance parameters
        // For nmix=1, use GetFactorVarianceIndex; for nmix>1, handled in mixture loop
        for (int ifac = 0; ifac < n_input_factors; ifac++) {
            factor_var[ifac] = std::fabs(param[GetFactorVarianceIndex(0, ifac)]);
        }

        // SE parameters: intercept, linear coefficients, [quadratic coefficients], residual variance
        se_intercept = param[GetSEInterceptIndex()];
        for (int j = 0; j < n_input_factors; j++) {
            se_linear_coef[j] = param[GetSELinearIndex(j)];
        }

        // SE_QUADRATIC: also extract quadratic coefficients
        if (factor_structure == FactorStructure::SE_QUADRATIC) {
            for (int j = 0; j < n_input_factors; j++) {
                se_quadratic_coef[j] = param[GetSEQuadraticIndex(j)];
            }
        }

        // Type-specific SE intercepts (when ntyp > 1)
        // se_type_offset[0] stays 0 (reference type); types 2..ntyp read from param vector.
        if (ntyp > 1 && nse_type_intercepts > 0) {
            for (int t = 1; t < ntyp; t++) {
                se_type_offset[t] = param[GetSETypeInterceptIndex(t - 1)];
            }
        }

        se_residual_var = std::fabs(param[GetSEResidualVarIndex()]);  // Enforce positivity
        sigma_eps = std::sqrt(se_residual_var);

        // SE covariate coefficients (if enabled)
        if (use_se_covariates) {
            se_covariate_coef.resize(n_se_covariates);
            for (int j = 0; j < n_se_covariates; j++) {
                se_covariate_coef[j] = param[se_covariate_param_start + j];
            }
        }

        // Outcome factor variance is derived from SE equation (not a free parameter)
        factor_var[nfac - 1] = se_residual_var;  // For integration weighting

    } else if (factor_structure == FactorStructure::CORRELATION && nfac == 2) {
        // Correlated factors: all nfac variances + 1 correlation
        // Note: mixtures not supported with correlation structure
        for (int ifac = 0; ifac < nfac; ifac++) {
            factor_var[ifac] = std::fabs(param[ifac]);
        }
        rho = param[nfac];  // Correlation parameter
        // Clamp to valid range to prevent numerical issues
        if (rho > 0.999) rho = 0.999;
        if (rho < -0.999) rho = -0.999;
        sqrt_1_minus_rho2 = std::sqrt(1.0 - rho * rho);

    } else {
        // INDEPENDENT: all nfac variances (use first mixture for single-mixture compatibility)
        for (int ifac = 0; ifac < nfac; ifac++) {
            factor_var[ifac] = std::fabs(param[GetFactorVarianceIndex(0, ifac)]);
        }
    }

    // ===== STEP 3: Initialize outputs =====
    logLkhd = 0.0;
    std::vector<double> full_gradL(nparam, 0.0);
    std::vector<double> full_hessL;
    if (iflag == 3) {
        full_hessL.resize(nparam * nparam, 0.0);
    }

    // ===== STEP 4: Compute default number of integration points =====
    // (may be overridden per-observation in adaptive mode)
    int nint_points_default = 1;
    for (int i = 0; i < nfac; i++) {
        nint_points_default *= nquad_points;
    }

    // Temporary storage
    // OPTIMIZATION: Pre-allocate modEval and modHess to maximum needed size
    // This avoids repeated resize operations inside the hot loop
    int max_model_params = 0;
    for (size_t imod = 0; imod < models.size(); imod++) {
        if (param_model_count[imod] > max_model_params) {
            max_model_params = param_model_count[imod];
        }
    }
    int max_ngrad = 1 + nfac + max_model_params;  // Maximum gradient size for any model
    int max_hess_dim = nfac + max_model_params;   // Maximum Hessian dimension
    std::vector<double> modEval(max_ngrad, 0.0);
    std::vector<double> modHess;
    if (iflag == 3) {
        modHess.resize(max_hess_dim * max_hess_dim, 0.0);
    }
    std::vector<double> totalgrad(nparam, 0.0);
    std::vector<double> totalhess;
    if (iflag == 3) totalhess.resize(nparam * nparam, 0.0);

    // Pre-allocate working vectors outside the loops for performance
    // These are reused across observations and integration points
    std::vector<double> fac_val(nfac);
    // Base outcome factor value (SE models, without any type-specific SE intercept shift)
    // Used to restore fac_val[nfac-1] at the top of each iteration of the type loop.
    double base_f_outcome = 0.0;
    // Enable shifting fac_val[nfac-1] by the type-specific SE intercept inside the type loop.
    // Only active for SE_LINEAR / SE_QUADRATIC with ntyp > 1 and non-adaptive integration.
    // Adaptive integration uses Stage 1 factor scores directly for the outcome factor, so
    // the SE equation (and thus type intercepts) does not apply there.
    const bool se_type_shift_enabled =
        (nse_type_intercepts > 0) && !use_adaptive &&
        (factor_structure == FactorStructure::SE_LINEAR ||
         factor_structure == FactorStructure::SE_QUADRATIC);
    std::vector<double> type_weighted_grad(nparam, 0.0);
    std::vector<double> type_weighted_hess;
    if (iflag == 3) type_weighted_hess.resize(nparam * nparam, 0.0);
    std::vector<double> grad_this_type(nparam, 0.0);
    std::vector<double> hess_this_type;
    if (iflag == 3) hess_this_type.resize(nparam * nparam, 0.0);

    // OPTIMIZATION: Pre-allocate per-observation vectors outside the loop
    // These were previously allocated inside the observation loop, causing
    // nobs vector allocations per call to CalcLkhd
    std::vector<int> obs_nq(nfac);           // Number of quad points per factor
    std::vector<double> fac_center(nfac);    // Center for factor integration
    std::vector<double> fac_spread(nfac);    // Spread for factor integration (SE in adaptive mode, SD otherwise)
    std::vector<int> facint(nfac);           // Integration point indices (odometer)
    std::vector<double> type_probs(ntyp);    // Type probabilities (avoid allocation in ComputeTypeProbabilities)

    // OPTIMIZATION: Precompute sigma_fac ONCE (sqrt of factor variances)
    // These are constant for all observations and all integration points
    std::vector<double> sigma_fac(nfac);
    for (int ifac = 0; ifac < nfac; ifac++) {
        sigma_fac[ifac] = std::sqrt(factor_var[ifac]);
    }

    // OPTIMIZATION: Precompute integration point tables for likelihood/gradient only
    // For Hessian (iflag==3), inline computation is faster due to better cache locality
    // This gives us: faster likelihood/gradient from precomputation, fast Hessian from inline
    // NOTE: Cannot use precomputed tables when:
    // - factor_mean_covariates is enabled (factor_mean varies per observation)
    // - use_se_covariates is enabled
    // - nmix > 1 (each mixture has different sigma/mean, would need tables per mixture)
    bool use_precomputed_tables = !use_adaptive && (iflag < 3) && !use_factor_mean_covariates && !use_se_covariates && (nmix == 1);

    std::vector<std::vector<double>> fac_val_table;
    std::vector<double> probilk_table;
    std::vector<std::vector<double>> df_dsigma2_table;

    if (use_precomputed_tables) {
        fac_val_table.resize(nint_points_default, std::vector<double>(nfac));
        probilk_table.resize(nint_points_default);
        if (iflag >= 2) {
            df_dsigma2_table.resize(nint_points_default, std::vector<double>(nfac));
        }

        // Use odometer method to iterate through all integration points
        std::vector<int> precomp_facint(nfac, 0);
        for (int intpt = 0; intpt < nint_points_default; intpt++) {
            // Compute weight (product over dimensions)
            double probilk = 1.0;
            for (int ifac = 0; ifac < nfac; ifac++) {
                probilk *= quad_weights[precomp_facint[ifac]];
            }
            probilk_table[intpt] = probilk;

            // Compute factor values at this integration point
            if (factor_structure == FactorStructure::SE_LINEAR ||
                factor_structure == FactorStructure::SE_QUADRATIC) {
                // SE models: integrate over (f_1, ..., f_{k-1}, epsilon)
                // First (nfac-1) dimensions are input factors
                for (int j = 0; j < n_input_factors; j++) {
                    double x_node = quad_nodes[precomp_facint[j]];
                    fac_val_table[intpt][j] = sigma_fac[j] * x_node + factor_mean[j];
                }
                // Last dimension is residual epsilon (mean 0)
                double x_eps = quad_nodes[precomp_facint[nfac - 1]];
                double eps = sigma_eps * x_eps;

                // Compute outcome factor: f_k = se_intercept + sum(se_linear_j * f_j) [+ sum(se_quadratic_j * f_j^2)] + eps
                double f_outcome = se_intercept + eps;
                for (int j = 0; j < n_input_factors; j++) {
                    f_outcome += se_linear_coef[j] * fac_val_table[intpt][j];
                    // SE_QUADRATIC: add quadratic terms
                    if (factor_structure == FactorStructure::SE_QUADRATIC) {
                        f_outcome += se_quadratic_coef[j] * fac_val_table[intpt][j] * fac_val_table[intpt][j];
                    }
                }
                fac_val_table[intpt][nfac - 1] = f_outcome;

            } else if (fac_corr && nfac == 2) {
                // Correlated 2-factor model: use Cholesky transformation
                double x0 = quad_nodes[precomp_facint[0]];
                double x1 = quad_nodes[precomp_facint[1]];
                fac_val_table[intpt][0] = sigma_fac[0] * x0 + factor_mean[0];
                fac_val_table[intpt][1] = rho * sigma_fac[1] * x0 + sigma_fac[1] * sqrt_1_minus_rho2 * x1 + factor_mean[1];
            } else {
                // Independent factors
                for (int ifac = 0; ifac < nfac; ifac++) {
                    double x_node = quad_nodes[precomp_facint[ifac]];
                    fac_val_table[intpt][ifac] = sigma_fac[ifac] * x_node + factor_mean[ifac];
                }
            }

            // Precompute chain rule factors for gradients
            if (iflag >= 2) {
                for (int ifac = 0; ifac < nfac; ifac++) {
                    double x_node = quad_nodes[precomp_facint[ifac]];
                    df_dsigma2_table[intpt][ifac] = x_node / (2.0 * sigma_fac[ifac]);
                }
            }

            // Odometer increment
            for (int ifac = 0; ifac < nfac; ifac++) {
                precomp_facint[ifac]++;
                if (precomp_facint[ifac] < nquad_points) break;
                precomp_facint[ifac] = 0;
            }
        }
    }

    // Working vectors for adaptive mode (per-integration-point values)
    std::vector<double> x_node_fac(nfac);
    std::vector<double> df_dsigma2(nfac);
    std::vector<double> d2f_dsigma2_sq(nfac);

    // ===== STEP 5: Loop over observations =====
    for (int iobs = 0; iobs < nobs; iobs++) {
        int iobs_offset = iobs * nvar;

        // Compute observation-specific factor means if covariates are set
        // factor_mean[k] = sum_j(gamma_{k,j} * X_i_j)
        if (use_factor_mean_covariates) {
            for (int ifac = 0; ifac < n_factors_with_mean; ifac++) {
                double mean_k = 0.0;
                for (int icov = 0; icov < n_factor_mean_covariates; icov++) {
                    // Parameter index: factor_mean_param_start + ifac * n_covariates + icov
                    int param_idx = factor_mean_param_start + ifac * n_factor_mean_covariates + icov;
                    double gamma_k_j = param[param_idx];
                    double x_ij = factor_mean_covariate_data[iobs][icov];
                    mean_k += gamma_k_j * x_ij;
                }
                factor_mean[ifac] = mean_k;
            }
        }

        double totalprob = 0.0;
        // Zero gradient entries for params that need gradients (free + tied)
        if (iflag >= 2) {
            for (size_t gi = 0; gi < gradparlist.size(); gi++) {
                totalgrad[gradparlist[gi]] = 0.0;
            }
        }
        // Zero Hessian entries for every parameter that needs a Hessian
        // contribution. With equality constraints this includes tied
        // (derived) parameters because ExtractFreeHessian later aggregates
        // their cross-derivatives into their primaries; the tied positions
        // must therefore be populated on the upper triangle.
        if (iflag == 3) {
            for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                int i = gradparlist[fi];
                for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                    int j = gradparlist[fj];
                    totalhess[i * nparam + j] = 0.0;
                }
            }
        }

        // ===== STEP 6: Multi-dimensional integration over factors =====
        // For mixtures (nmix > 1), loop over mixture components
        // L = sum_m w_m * L_m where L_m is likelihood for mixture m

        // Model parameter indexing (needed for mixture Hessians)
        int first_model_param = (models.size() > 0) ? param_model_start[0] : nparam;
        int n_model_params = nparam - first_model_param;

        // Storage for per-mixture likelihoods and gradients (needed for mixture param derivatives)
        std::vector<double> mix_L(nmix, 0.0);  // L_m for each mixture
        std::vector<std::vector<double>> mix_dL_dmean;  // dL_m/d(factor_mean_k) for each mixture
        std::vector<std::vector<double>> mix_dL_dsigma;  // dL_m/d(sigma^2_k) for each mixture (for variance-weight cross)
        std::vector<std::vector<double>> mix_dL_dbeta;  // dL_m/d(beta_j) for each mixture (for weight-beta cross Hessian)
        if (nmix > 1 && iflag >= 2) {
            mix_dL_dmean.resize(nmix, std::vector<double>(n_factors_for_mixture, 0.0));
            mix_dL_dsigma.resize(nmix, std::vector<double>(n_factors_for_mixture, 0.0));
        }
        if (nmix > 1 && iflag == 3) {
            mix_dL_dbeta.resize(nmix, std::vector<double>(n_model_params, 0.0));
        }

        // For mixture Hessians: accumulate d²L_m/d(f_k1)d(f_k2) for each mixture
        // This is the factor-factor Hessian needed for mean-mean and cross-Hessian terms
        // Storage: mix_d2L_dfdf[imix][k1 * n_factors_for_mixture + k2]
        std::vector<std::vector<double>> mix_d2L_dfdf;
        // For variance-mean/variance-weight cross Hessians: accumulate d²L_m/d(σ²_k1)d(f_k2)
        // = d²L_m/d(f_k1)d(f_k2) * df_k1/d(σ²_k1) = d2L_dfdf * x_k1/(2*sigma_k1)
        // Storage: mix_d2L_dsigma_df[imix][k1 * n_factors_for_mixture + k2]
        std::vector<std::vector<double>> mix_d2L_dsigma_df;
        // For mean × model_param and weight × model_param cross-Hessians:
        // accumulate d²L_m/d(f_k)d(beta_j) for each mixture
        // Storage: mix_d2L_df_dbeta[imix][k * n_model_params + (param_idx - first_model_param)]
        std::vector<std::vector<double>> mix_d2L_df_dbeta;

        if (nmix > 1 && iflag == 3 && factor_structure == FactorStructure::INDEPENDENT) {
            int nfac_sq = n_factors_for_mixture * n_factors_for_mixture;
            int n_df_dbeta = n_factors_for_mixture * n_model_params;
            // IMPORTANT: Use assign() to force zero-initialization even when size unchanged
            // resize() only zeros new elements, not existing ones!
            mix_d2L_dfdf.resize(nmix);
            mix_d2L_dsigma_df.resize(nmix);
            mix_d2L_df_dbeta.resize(nmix);
            for (int m = 0; m < nmix; m++) {
                mix_d2L_dfdf[m].assign(nfac_sq, 0.0);
                mix_d2L_dsigma_df[m].assign(nfac_sq, 0.0);
                mix_d2L_df_dbeta[m].assign(n_df_dbeta, 0.0);
            }
        }

        for (int imix = 0; imix < nmix; imix++) {
            double mix_weight = mix_weights[imix];

            // Set sigma_fac and factor_mean for this mixture component
            // For SE models: only input factors have mixture distribution
            if (nmix > 1) {
                for (int k = 0; k < n_factors_for_mixture; k++) {
                    sigma_fac[k] = mix_sigmas[imix][k];
                    // Factor means from mixture (add to any covariate-based means)
                    if (!use_factor_mean_covariates) {
                        factor_mean[k] = mix_means[imix][k];
                    } else {
                        // With covariates: factor_mean already has covariate contribution
                        // Add mixture mean on top (mixture mean is relative to overall mean=0)
                        // Actually, with covariates, we typically don't use mixture means
                        // For now, just use mixture mean directly (covariates shift mean, mixture adds offset)
                        factor_mean[k] = mix_means[imix][k];  // Covariate effect applied separately
                    }
                }
            }

        // Determine per-observation integration settings
        // OPTIMIZATION: obs_nq and fac_center are pre-allocated outside the observation loop
        int nint_points;

        if (use_adaptive) {
            // Adaptive mode: use per-observation quadrature settings
            // Uses factor score SE as the spread (from importance sampling perspective)
            nint_points = 1;
            for (int ifac = 0; ifac < nfac; ifac++) {
                obs_nq[ifac] = obs_nquad[iobs][ifac];
                fac_center[ifac] = obs_fac_center[iobs][ifac];
                fac_spread[ifac] = obs_fac_se[iobs][ifac];  // SE for importance sampling
                nint_points *= obs_nq[ifac];
            }
        } else {
            // Standard mode: use global quadrature settings
            nint_points = nint_points_default;
            for (int ifac = 0; ifac < nfac; ifac++) {
                obs_nq[ifac] = nquad_points;
                fac_center[ifac] = factor_mean[ifac];
                fac_spread[ifac] = sigma_fac[ifac];  // Full SD for standard integration
            }
        }

        // Accumulators for this mixture component
        double mixprob = 0.0;
        std::vector<double> mixgrad(nparam, 0.0);
        std::vector<double> mixhess;
        if (iflag == 3) mixhess.resize(nparam * nparam, 0.0);

        // For mixture mean gradients: accumulate dL_m/d(mean_k) = Σ w_q * L_q * d(log L_q)/df_k
        // This is the "raw" gradient w.r.t. factor mean (df/d(mean) = 1)
        std::vector<double> mix_dL_df(n_factors_for_mixture, 0.0);

        // Use "odometer" method: facint tracks indices for each dimension
        // OPTIMIZATION: facint is pre-allocated outside the observation loop, just reset to 0
        std::fill(facint.begin(), facint.end(), 0);

        for (int intpt = 0; intpt < nint_points; intpt++) {
            // OPTIMIZATION: Three code paths for different scenarios:
            // 1. use_precomputed_tables: likelihood/gradient in non-adaptive mode (use precomputed tables)
            // 2. !use_adaptive && !use_precomputed_tables: Hessian in non-adaptive mode (inline computation)
            // 3. use_adaptive: any computation in adaptive mode (per-observation settings)
            double probilk;

            if (use_precomputed_tables) {
                // Non-adaptive likelihood/gradient: use precomputed tables
                probilk = probilk_table[intpt];
                for (int ifac = 0; ifac < nfac; ifac++) {
                    fac_val[ifac] = fac_val_table[intpt][ifac];
                    x_node_fac[ifac] = quad_nodes[facint[ifac]];
                }
                if (iflag >= 2) {
                    for (int ifac = 0; ifac < nfac; ifac++) {
                        df_dsigma2[ifac] = df_dsigma2_table[intpt][ifac];
                    }
                }
            } else if (!use_adaptive) {
                // Non-adaptive Hessian: compute inline for better cache locality
                // This avoids table allocation overhead that hurts Hessian performance
                probilk = 1.0;
                for (int ifac = 0; ifac < nfac; ifac++) {
                    probilk *= quad_weights[facint[ifac]];
                    x_node_fac[ifac] = quad_nodes[facint[ifac]];
                }

                // Compute factor values
                if (factor_structure == FactorStructure::SE_LINEAR ||
                    factor_structure == FactorStructure::SE_QUADRATIC) {
                    // SE models: integrate over (f_1, ..., f_{k-1}, epsilon)
                    for (int j = 0; j < n_input_factors; j++) {
                        fac_val[j] = sigma_fac[j] * x_node_fac[j] + factor_mean[j];
                    }
                    // Last dimension is residual epsilon
                    double x_eps = x_node_fac[nfac - 1];
                    double eps = sigma_eps * x_eps;

                    // Compute outcome factor
                    double f_outcome = se_intercept + eps;
                    for (int j = 0; j < n_input_factors; j++) {
                        f_outcome += se_linear_coef[j] * fac_val[j];
                        // SE_QUADRATIC: add quadratic terms
                        if (factor_structure == FactorStructure::SE_QUADRATIC) {
                            f_outcome += se_quadratic_coef[j] * fac_val[j] * fac_val[j];
                        }
                    }
                    // Add SE covariate contributions
                    if (use_se_covariates) {
                        for (int m = 0; m < n_se_covariates; m++) {
                            f_outcome += se_covariate_coef[m] * se_covariate_data[iobs][m];
                        }
                    }
                    fac_val[nfac - 1] = f_outcome;

                } else if (fac_corr && nfac == 2) {
                    double x0 = x_node_fac[0];
                    double x1 = x_node_fac[1];
                    fac_val[0] = sigma_fac[0] * x0 + factor_mean[0];
                    fac_val[1] = rho * sigma_fac[1] * x0 + sigma_fac[1] * sqrt_1_minus_rho2 * x1 + factor_mean[1];
                } else {
                    for (int ifac = 0; ifac < nfac; ifac++) {
                        fac_val[ifac] = sigma_fac[ifac] * x_node_fac[ifac] + factor_mean[ifac];
                    }
                }

                // Compute chain rule factors for gradient/Hessian
                for (int ifac = 0; ifac < nfac; ifac++) {
                    df_dsigma2[ifac] = x_node_fac[ifac] / (2.0 * sigma_fac[ifac]);
                }
                if (iflag == 3) {
                    for (int ifac = 0; ifac < nfac; ifac++) {
                        double sigma_cubed = sigma_fac[ifac] * sigma_fac[ifac] * sigma_fac[ifac];
                        d2f_dsigma2_sq[ifac] = -x_node_fac[ifac] / (4.0 * sigma_cubed);
                    }
                }
            } else {
                // Adaptive mode: compute per-observation
                // Get quadrature nodes for this observation
                for (int ifac = 0; ifac < nfac; ifac++) {
                    int nq = obs_nq[ifac];
                    x_node_fac[ifac] = adapt_nodes.at(nq)[facint[ifac]];
                }

                // Compute weight (product over dimensions)
                probilk = 1.0;
                for (int ifac = 0; ifac < nfac; ifac++) {
                    int nq = obs_nq[ifac];
                    probilk *= adapt_weights.at(nq)[facint[ifac]];
                }

                // Compute factor values using SE as spread (matching legacy code)
                // f = fac_center + fac_spread * x_node
                // where fac_center = factor score and fac_spread = SE (from Stage 1)
                for (int ifac = 0; ifac < nfac; ifac++) {
                    fac_val[ifac] = fac_spread[ifac] * x_node_fac[ifac] + fac_center[ifac];
                }

                // Apply importance sampling correction
                // Mathematical derivation:
                //   We want: ∫ L(f) N(f|0,σ²) df
                //   Using proposal N(f|fscore, SE²):
                //   = ∫ L(f) [N(f|0,σ²)/N(f|fscore,SE²)] N(f|fscore,SE²) df
                //   = Σ wᵢ L(fscore + SE·zᵢ) × IS(zᵢ)
                // where IS = N(f|0,σ²) / N(f|fscore,SE²) = (SE/σ) × exp(z²/2 - f²/(2σ²))
                //
                // The jacobian SE/σ arises from the ratio of the normal density normalizing constants
                for (int ifac = 0; ifac < nfac; ifac++) {
                    double x_node = x_node_fac[ifac];
                    double f = fac_val[ifac];
                    double se = fac_spread[ifac];  // Factor score SE from Stage 1
                    double sd = sigma_fac[ifac];   // sqrt(factor variance parameter)
                    double var = sd * sd;

                    // IS correction = (SE/σ) × exp(z²/2 - f²/(2σ²))
                    double jacobian = se / sd;
                    double is_correction = jacobian * std::exp(x_node * x_node / 2.0 - f * f / (2.0 * var));

                    probilk *= is_correction;
                }

                // Compute chain rule factors for gradients/Hessians in adaptive mode
                if (iflag >= 2) {
                    for (int ifac = 0; ifac < nfac; ifac++) {
                        df_dsigma2[ifac] = x_node_fac[ifac] / (2.0 * sigma_fac[ifac]);
                    }
                }
                if (iflag == 3) {
                    for (int ifac = 0; ifac < nfac; ifac++) {
                        double sigma_cubed = sigma_fac[ifac] * sigma_fac[ifac] * sigma_fac[ifac];
                        d2f_dsigma2_sq[ifac] = -x_node_fac[ifac] / (4.0 * sigma_cubed);
                    }
                }
            }

            // ===== STEP 8: Evaluate all models with type mixture =====
            // For ntyp > 1: L = Σ_t π_t(f) × L_t(y|f) where L_t = Π_m L_m(y_m|f, intercept_t)
            // For ntyp == 1: Standard case (no type mixture)

            // Save the base (reference-type) outcome factor value so we can apply
            // type-specific SE intercept shifts inside the type loop below without losing
            // the reference value. Type probabilities are computed using the base outcome
            // (the user-facing semantic: type is drawn before the residual epsilon and
            // does not depend on the realized type-specific outcome).
            if (se_type_shift_enabled) {
                base_f_outcome = fac_val[nfac - 1];
            }

            // Compute type probabilities (fills pre-allocated type_probs vector)
            // OPTIMIZATION: type_probs is pre-allocated outside the observation loop
            ComputeTypeProbabilities(fac_val, type_probs);

            // Accumulators for type-weighted contributions (pre-allocated, zero here)
            // OPTIMIZATION: For ntyp=1, skip zeroing type_weighted_* - we accumulate directly into total*
            double type_weighted_prob = 0.0;
            if (ntyp > 1) {
                if (iflag >= 2) std::fill(type_weighted_grad.begin(), type_weighted_grad.end(), 0.0);
                // Include tied params under equality constraints (see comment
                // at the outer totalhess zeroing loop).
                if (iflag == 3) {
                    for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                        int i = gradparlist[fi];
                        for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                            int j = gradparlist[fj];
                            type_weighted_hess[i * nparam + j] = 0.0;
                        }
                    }
                }
            }

            // Likelihood product for this type (declared outside loop for ntyp=1 optimization)
            double prob_this_type = 1.0;

            // For mixture Hessians: accumulate factor-factor log-likelihood Hessian and gradient
            // directly from modHess/modEval. This avoids numerical issues when extracting from
            // variance Hessian (which fails when df/dsigma² is near zero).
            // d2logL_dfdf[k1*n_factors_for_mixture+k2] = d²(log L)/d(f_{k1})d(f_{k2})
            // dlogL_df[k] = d(log L)/d(f_k)
            std::vector<double> d2logL_dfdf;
            std::vector<double> dlogL_df;
            // For factor × model_param LOG-Hessian at this integration point
            // d2logL_df_dbeta_thisintpt[k * n_model_params + (param - first_model_param)]
            std::vector<double> d2logL_df_dbeta_thisintpt;
            if (nmix > 1 && iflag == 3 && factor_structure == FactorStructure::INDEPENDENT && ntyp == 1) {
                int nfac_sq = n_factors_for_mixture * n_factors_for_mixture;
                d2logL_dfdf.resize(nfac_sq, 0.0);
                dlogL_df.resize(n_factors_for_mixture, 0.0);
                d2logL_df_dbeta_thisintpt.resize(n_factors_for_mixture * n_model_params, 0.0);
            }

            // Loop over types
            for (int ityp = 0; ityp < ntyp; ityp++) {
                double type_prob = type_probs[ityp];

                // Apply type-specific SE intercept shift to the outcome factor.
                // f_outcome_t = base_f_outcome + se_type_offset[ityp]
                // where se_type_offset[0] = 0 (reference type). All subsequent
                // model evaluations inside this type iteration see the shifted outcome.
                if (se_type_shift_enabled) {
                    fac_val[nfac - 1] = base_f_outcome + se_type_offset[ityp];
                }

                // Reset likelihood product for this type
                prob_this_type = 1.0;

                // Gradient accumulators for this type (pre-allocated, zero here)
                if (iflag >= 2) std::fill(grad_this_type.begin(), grad_this_type.end(), 0.0);
                // Include tied params under equality constraints (see comment
                // at the outer totalhess zeroing loop).
                if (iflag == 3) {
                    for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                        int i = gradparlist[fi];
                        for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                            int j = gradparlist[fj];
                            hess_this_type[i * nparam + j] = 0.0;
                        }
                    }
                }

                // Evaluate all models for this type
                for (size_t imod = 0; imod < models.size(); imod++) {
                    int firstpar = param_model_start[imod];

                    // Get type-specific intercept for types > 1 (only if model uses types)
                    // Type 1 (ityp=0) is reference with intercept = 0
                    double type_intercept = 0.0;
                    if (ityp > 0 && ntyp > 1 && models[imod]->GetUseTypes()) {
                        // Type intercepts are stored after base model parameters
                        // ityp-1 because type 1 is reference (no intercept)
                        int intercept_offset = GetTypeInterceptIndex(ityp - 1, imod);
                        type_intercept = param[firstpar + intercept_offset];
                    }

                    // If all parameters are fixed for this model, only compute likelihood (flag=1)
                    // EXCEPTION: For SE models, we still need factor gradients (dL/df_i) even when
                    // model params are fixed, because SE parameters depend on factor values.
                    // So for SE models, use iflag (to get factor gradients) but skip param gradients
                    // in the model (handled by setting all_params_fixed flag).
                    int model_flag;
                    if (models[imod]->GetAllParamsFixed()) {
                        if (factor_structure == FactorStructure::SE_LINEAR ||
                            factor_structure == FactorStructure::SE_QUADRATIC) {
                            // SE models: need factor gradients AND factor-factor Hessian
                            // for computing SE parameter derivatives.
                            // The model params are fixed, but we still need modHess for
                            // the factor-factor block (d2L/df_i*df_j) used by SE Hessian.
                            model_flag = iflag;
                        } else {
                            // Non-SE: all params fixed means no gradients needed
                            model_flag = 1;
                        }
                    } else {
                        model_flag = iflag;
                    }

                    // Pass model-local free indices for Hessian optimization
                    // Only pass when SOME (not all) parameters are fixed - avoids overhead when all are free
                    const std::vector<int>* free_indices_ptr = nullptr;
                    if (model_flag == 3 && imod < model_free_indices.size() &&
                        !model_free_indices[imod].empty() &&
                        model_free_indices[imod].size() < static_cast<size_t>(param_model_count[imod])) {
                        free_indices_ptr = &model_free_indices[imod];
                    }

                    models[imod]->Eval(iobs_offset, data, param, firstpar, fac_val,
                                      modEval, modHess, model_flag, type_intercept, free_indices_ptr);

                    // Multiply likelihood for this type
                    prob_this_type *= modEval[0];

                    // ===== STEP 9: Accumulate gradients for this type =====
                    if (iflag >= 2) {
                        // Gradients w.r.t. factor variances and structure-specific parameters
                        if (factor_structure == FactorStructure::SE_LINEAR ||
                            factor_structure == FactorStructure::SE_QUADRATIC) {
                            // SE models: f_k = se_intercept + Σ_j se_linear_j * f_j [+ Σ_j se_quadratic_j * f_j²] + eps
                            // modEval contains [dens, dL/df_1, ..., dL/df_k, dL/dparams...]
                            double dL_dfk = modEval[nfac];  // Gradient w.r.t. outcome factor

                            // Gradients w.r.t. input factor variances
                            // ∂f_k/∂f_j = α_j + 2*α_qj*f_j  (for SE_QUADRATIC)
                            // ∂L/∂σ²_j = [∂L/∂f_j + (∂L/∂f_k) * ∂f_k/∂f_j] * (x_j / (2*σ_j))
                            for (int j = 0; j < n_input_factors; j++) {
                                double dL_dfj = modEval[j + 1];
                                double dfk_dfj = se_linear_coef[j];
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    dfk_dfj += 2.0 * se_quadratic_coef[j] * fac_val[j];
                                }
                                double total_deriv = dL_dfj + dL_dfk * dfk_dfj;
                                grad_this_type[j] += total_deriv * df_dsigma2[j];
                            }

                            // Gradients w.r.t. SE parameters
                            // ∂L/∂(se_intercept) = ∂L/∂f_k * 1
                            grad_this_type[GetSEInterceptIndex()] += dL_dfk;

                            // ∂L/∂(se_intercept_type_{ityp}) = ∂L/∂f_k * 1
                            // Only non-reference types have free intercepts. The parameter
                            // only appears in this type's likelihood, so contributions from
                            // other types are zero and we only accumulate when ityp > 0.
                            if (ntyp > 1 && ityp > 0 && nse_type_intercepts > 0) {
                                grad_this_type[GetSETypeInterceptIndex(ityp - 1)] += dL_dfk;
                            }

                            // ∂L/∂(se_linear_j) = ∂L/∂f_k * f_j
                            for (int j = 0; j < n_input_factors; j++) {
                                grad_this_type[GetSELinearIndex(j)] += dL_dfk * fac_val[j];
                            }

                            // SE_QUADRATIC: ∂L/∂(se_quadratic_j) = ∂L/∂f_k * f_j²
                            if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                for (int j = 0; j < n_input_factors; j++) {
                                    grad_this_type[GetSEQuadraticIndex(j)] += dL_dfk * fac_val[j] * fac_val[j];
                                }
                            }

                            // ∂L/∂(se_residual_var) = ∂L/∂f_k * (x_eps / (2*sigma_eps))
                            double x_eps = x_node_fac[nfac - 1];
                            grad_this_type[GetSEResidualVarIndex()] += dL_dfk * x_eps / (2.0 * sigma_eps);

                            // SE covariate gradients: ∂L/∂β_m = ∂L/∂f_k * X_m
                            if (use_se_covariates) {
                                for (int m = 0; m < n_se_covariates; m++) {
                                    int param_idx = se_covariate_param_start + m;
                                    grad_this_type[param_idx] += dL_dfk * se_covariate_data[iobs][m];
                                }
                            }

                        } else if (fac_corr && nfac == 2) {
                            // Correlated 2-factor model: use Cholesky derivatives
                            // Use pre-computed sigma_fac and x_node_fac
                            double dL_df1 = modEval[1];
                            double dL_df2 = modEval[2];

                            grad_this_type[0] += dL_df1 * x_node_fac[0] / (2.0 * sigma_fac[0]);
                            double df2_dsigma2sq = (rho * x_node_fac[0] + sqrt_1_minus_rho2 * x_node_fac[1]) / (2.0 * sigma_fac[1]);
                            grad_this_type[1] += dL_df2 * df2_dsigma2sq;
                            double df2_drho = sigma_fac[1] * (x_node_fac[0] - rho * x_node_fac[1] / sqrt_1_minus_rho2);
                            grad_this_type[2] += dL_df2 * df2_drho;
                        } else {
                            // Independent factors: use pre-computed df_dsigma2
                            // For mixtures: accumulate to the correct variance parameter index
                            for (int ifac = 0; ifac < nfac; ifac++) {
                                // For factors with mixture distribution, use mixture-specific index
                                int var_idx = (nmix > 1 && ifac < n_factors_for_mixture) ?
                                    GetFactorVarianceIndex(imix, ifac) : ifac;
                                grad_this_type[var_idx] += modEval[ifac + 1] * df_dsigma2[ifac];
                            }
                            // Note: Mixture mean and log-weight gradients are computed separately in STEP 11d
                            // DO NOT add them to grad_this_type here to avoid double-counting
                        }

                        // Gradients w.r.t. factor mean covariate parameters
                        // d(log L)/d(gamma_{k,j}) = d(log L)/d(f_k) * X_i_j
                        if (use_factor_mean_covariates) {
                            for (int ifac = 0; ifac < n_factors_with_mean; ifac++) {
                                double dL_dfac = modEval[ifac + 1];
                                for (int icov = 0; icov < n_factor_mean_covariates; icov++) {
                                    int param_idx = factor_mean_param_start + ifac * n_factor_mean_covariates + icov;
                                    double x_ij = factor_mean_covariate_data[iobs][icov];
                                    grad_this_type[param_idx] += dL_dfac * x_ij;
                                }
                            }
                        }

                        // Gradients w.r.t. model parameters (base params only, not type intercepts)
                        // param_model_count includes type intercepts (if model uses types), but modEval only has base params
                        bool model_uses_types = models[imod]->GetUseTypes();
                        int n_type_intercepts = (ntyp > 1 && model_uses_types) ? (ntyp - 1) : 0;
                        int base_param_count = param_model_count[imod] - n_type_intercepts;
                        for (int iparam = 0; iparam < base_param_count; iparam++) {
                            grad_this_type[firstpar + iparam] += modEval[1 + nfac + iparam];
                        }

                        // Gradient w.r.t. this type's intercept (only for types > 1 and if model uses types)
                        // The gradient equals dL/d(linear_predictor) * 1
                        // For models with an intercept covariate (=1), this equals the base intercept gradient
                        if (ntyp > 1 && ityp > 0 && model_uses_types) {
                            int type_intercept_idx = firstpar + base_param_count + (ityp - 1);
                            // modEval[1 + nfac + 0] is gradient w.r.t. first covariate coefficient
                            // When first covariate is intercept (=1), this equals dL/d(linear_predictor)
                            grad_this_type[type_intercept_idx] += modEval[1 + nfac + 0];
                        }
                    }

                    // ===== STEP 9: Accumulate Hessians for this type =====
                    if (iflag == 3 && modHess.size() > 0) {
                        // Note: modHess from Model::Eval() only includes base parameters,
                        // not type-specific intercepts
                        bool model_uses_types_hess = models[imod]->GetUseTypes();
                        int n_type_intercepts_hess = (ntyp > 1 && model_uses_types_hess) ? (ntyp - 1) : 0;
                        int base_param_count_hess = param_model_count[imod] - n_type_intercepts_hess;
                        int nDimModHess = nfac + base_param_count_hess;

                        if (factor_structure == FactorStructure::SE_LINEAR ||
                            factor_structure == FactorStructure::SE_QUADRATIC) {
                            // ===== SE model Hessian =====
                            // SE_LINEAR: f_k = α + Σ_j α_j * f_j + ε
                            // SE_QUADRATIC: f_k = α + Σ_j α_j * f_j + Σ_j α_qj * f_j² + ε
                            // Need chain rule: ∂²L/∂θ₁∂θ₂ = Σ_i Σ_l (∂²L/∂f_i∂f_l)(∂f_i/∂θ₁)(∂f_l/∂θ₂) + Σ_i (∂L/∂f_i)(∂²f_i/∂θ₁∂θ₂)

                            double dL_dfk = modEval[nfac];  // Gradient w.r.t. outcome factor
                            double x_eps = x_node_fac[nfac - 1];

                            // Pre-compute derivatives
                            // For SE_QUADRATIC: ∂f_k/∂f_j = α_j + 2*α_qj*f_j, ∂²f_k/∂f_j² = 2*α_qj
                            // For SE_LINEAR: ∂f_k/∂f_j = α_j, ∂²f_k/∂f_j² = 0
                            std::vector<double> dfk_dfj(n_input_factors);
                            std::vector<double> d2fk_dfj2(n_input_factors, 0.0);
                            for (int j = 0; j < n_input_factors; j++) {
                                dfk_dfj[j] = se_linear_coef[j];
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    dfk_dfj[j] += 2.0 * se_quadratic_coef[j] * fac_val[j];
                                    d2fk_dfj2[j] = 2.0 * se_quadratic_coef[j];
                                }
                            }

                            double dfk_dres_var = x_eps / (2.0 * sigma_eps);
                            double sigma_eps_cubed = sigma_eps * sigma_eps * sigma_eps;
                            double d2fk_dres_var2 = -x_eps / (4.0 * sigma_eps_cubed);

                            // ===== Input factor variance terms (param indices 0..n_input_factors-1) =====
                            for (int i = 0; i < n_input_factors; i++) {
                                double dL_dfi = modEval[i + 1];
                                double dfk_dvar_i = dfk_dfj[i] * df_dsigma2[i];  // ∂f_k/∂σ²_i

                                // Diagonal: ∂²L/∂σ²_i²
                                // = (∂²L/∂f_i²)(∂f_i/∂σ²_i)² + (∂L/∂f_i)(∂²f_i/∂σ²_i²)
                                //   + 2(∂²L/∂f_i∂f_k)(∂f_i/∂σ²_i)(∂f_k/∂σ²_i) + (∂L/∂f_k)(∂f_k/∂f_i)(∂²f_i/∂σ²_i²)
                                //   + (∂²L/∂f_k²)(∂f_k/∂σ²_i)² + (∂L/∂f_k)(∂²f_k/∂f_i²)(∂f_i/∂σ²_i)² (for SE_QUADRATIC)
                                double d2L_dfidfi = modHess[i * nDimModHess + i];
                                double d2L_dfidfk = modHess[i * nDimModHess + (nfac - 1)];
                                double d2L_dfkdfk = modHess[(nfac - 1) * nDimModHess + (nfac - 1)];

                                double contrib = d2L_dfidfi * df_dsigma2[i] * df_dsigma2[i]
                                    + dL_dfi * d2f_dsigma2_sq[i]
                                    + 2.0 * d2L_dfidfk * df_dsigma2[i] * dfk_dvar_i
                                    + dL_dfk * dfk_dfj[i] * d2f_dsigma2_sq[i]
                                    + d2L_dfkdfk * dfk_dvar_i * dfk_dvar_i;
                                // SE_QUADRATIC: add term for ∂²f_k/∂f_i²
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    contrib += dL_dfk * d2fk_dfj2[i] * df_dsigma2[i] * df_dsigma2[i];
                                }
                                hess_this_type[i * nparam + i] += contrib;

                                // Cross-terms with other input factor variances
                                for (int j = i + 1; j < n_input_factors; j++) {
                                    double d2L_dfidfj = modHess[i * nDimModHess + j];
                                    double d2L_dfjdfk = modHess[j * nDimModHess + (nfac - 1)];
                                    double dfk_dvar_j = dfk_dfj[j] * df_dsigma2[j];

                                    hess_this_type[i * nparam + j] +=
                                        d2L_dfidfj * df_dsigma2[i] * df_dsigma2[j]
                                        + d2L_dfidfk * df_dsigma2[i] * dfk_dvar_j
                                        + d2L_dfjdfk * df_dsigma2[j] * dfk_dvar_i
                                        + d2L_dfkdfk * dfk_dvar_i * dfk_dvar_j;
                                }

                                // Cross-terms with SE parameters
                                int se_int_idx = GetSEInterceptIndex();
                                // ∂²L/∂σ²_i∂(se_intercept) = (∂²L/∂f_i∂f_k)(∂f_i/∂σ²_i)(1) + (∂²L/∂f_k²)(∂f_k/∂σ²_i)(1)
                                hess_this_type[i * nparam + se_int_idx] +=
                                    d2L_dfidfk * df_dsigma2[i] + d2L_dfkdfk * dfk_dvar_i;

                                for (int j = 0; j < n_input_factors; j++) {
                                    int se_lin_idx = GetSELinearIndex(j);
                                    // ∂²L/∂σ²_i∂(se_linear_j) = (∂²L/∂f_i∂f_k)(∂f_i/∂σ²_i)(f_j) + (∂²L/∂f_k²)(∂f_k/∂σ²_i)(f_j)
                                    //   + (∂L/∂f_k)(∂f_j/∂σ²_i) if i==j (second derivative of f_k w.r.t. σ²_i and α_i)
                                    double lin_contrib = d2L_dfidfk * df_dsigma2[i] * fac_val[j]
                                                   + d2L_dfkdfk * dfk_dvar_i * fac_val[j];
                                    if (i == j) {
                                        lin_contrib += dL_dfk * df_dsigma2[i];
                                    }
                                    hess_this_type[i * nparam + se_lin_idx] += lin_contrib;
                                }

                                // SE_QUADRATIC: cross-terms with quadratic coefficients
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    for (int j = 0; j < n_input_factors; j++) {
                                        int se_quad_idx = GetSEQuadraticIndex(j);
                                        // ∂²L/∂σ²_i∂(se_quadratic_j) = (∂²L/∂f_i∂f_k)(∂f_i/∂σ²_i)(f_j²) + (∂²L/∂f_k²)(∂f_k/∂σ²_i)(f_j²)
                                        //   + (∂L/∂f_k)(2*f_j)(∂f_j/∂σ²_i) if i==j
                                        double quad_contrib = d2L_dfidfk * df_dsigma2[i] * fac_val[j] * fac_val[j]
                                                       + d2L_dfkdfk * dfk_dvar_i * fac_val[j] * fac_val[j];
                                        if (i == j) {
                                            quad_contrib += dL_dfk * 2.0 * fac_val[j] * df_dsigma2[i];
                                        }
                                        hess_this_type[i * nparam + se_quad_idx] += quad_contrib;
                                    }
                                }

                                int se_res_idx = GetSEResidualVarIndex();
                                // ∂²L/∂σ²_i∂(se_residual_var) = (∂²L/∂f_i∂f_k)(∂f_i/∂σ²_i)(∂f_k/∂σ²_ε) + (∂²L/∂f_k²)(∂f_k/∂σ²_i)(∂f_k/∂σ²_ε)
                                hess_this_type[i * nparam + se_res_idx] +=
                                    d2L_dfidfk * df_dsigma2[i] * dfk_dres_var + d2L_dfkdfk * dfk_dvar_i * dfk_dres_var;

                                // ∂²L/∂σ²_i∂β_m = (∂²L/∂f_i∂f_k)(∂f_i/∂σ²_i) * X_m + (∂²L/∂f_k²)(∂f_k/∂σ²_i) * X_m
                                if (use_se_covariates) {
                                    for (int m = 0; m < n_se_covariates; m++) {
                                        int se_cov_idx = se_covariate_param_start + m;
                                        double X_m = se_covariate_data[iobs][m];
                                        hess_this_type[i * nparam + se_cov_idx] +=
                                            d2L_dfidfk * df_dsigma2[i] * X_m + d2L_dfkdfk * dfk_dvar_i * X_m;
                                    }
                                }
                            }

                            // ===== SE parameter terms =====
                            int se_int_idx = GetSEInterceptIndex();
                            int se_res_idx = GetSEResidualVarIndex();
                            double d2L_dfkdfk = modHess[(nfac - 1) * nDimModHess + (nfac - 1)];

                            // SE intercept diagonal: (∂²L/∂f_k²) * 1 * 1
                            hess_this_type[se_int_idx * nparam + se_int_idx] += d2L_dfkdfk;

                            // SE intercept × SE linear cross-terms
                            for (int j = 0; j < n_input_factors; j++) {
                                int se_lin_idx = GetSELinearIndex(j);
                                // ∂²L/∂(se_intercept)∂(se_linear_j) = (∂²L/∂f_k²) * 1 * f_j
                                hess_this_type[se_int_idx * nparam + se_lin_idx] += d2L_dfkdfk * fac_val[j];
                            }

                            // SE_QUADRATIC: SE intercept × SE quadratic cross-terms
                            if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                for (int j = 0; j < n_input_factors; j++) {
                                    int se_quad_idx = GetSEQuadraticIndex(j);
                                    hess_this_type[se_int_idx * nparam + se_quad_idx] += d2L_dfkdfk * fac_val[j] * fac_val[j];
                                }
                            }

                            // SE intercept × SE residual var
                            hess_this_type[se_int_idx * nparam + se_res_idx] += d2L_dfkdfk * dfk_dres_var;

                            // SE linear terms
                            for (int i = 0; i < n_input_factors; i++) {
                                int se_lin_i = GetSELinearIndex(i);
                                // Diagonal: (∂²L/∂f_k²) * f_i * f_i
                                hess_this_type[se_lin_i * nparam + se_lin_i] += d2L_dfkdfk * fac_val[i] * fac_val[i];

                                // Cross-terms with other SE linear
                                for (int j = i + 1; j < n_input_factors; j++) {
                                    int se_lin_j = GetSELinearIndex(j);
                                    hess_this_type[se_lin_i * nparam + se_lin_j] += d2L_dfkdfk * fac_val[i] * fac_val[j];
                                }

                                // SE_QUADRATIC: SE linear × SE quadratic cross-terms
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    for (int j = 0; j < n_input_factors; j++) {
                                        int se_quad_idx = GetSEQuadraticIndex(j);
                                        // ∂²L/∂(se_linear_i)∂(se_quadratic_j) = (∂²L/∂f_k²) * f_i * f_j²
                                        hess_this_type[se_lin_i * nparam + se_quad_idx] += d2L_dfkdfk * fac_val[i] * fac_val[j] * fac_val[j];
                                    }
                                }

                                // SE linear × SE residual var
                                hess_this_type[se_lin_i * nparam + se_res_idx] += d2L_dfkdfk * fac_val[i] * dfk_dres_var;
                            }

                            // SE_QUADRATIC: SE quadratic terms
                            if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                for (int i = 0; i < n_input_factors; i++) {
                                    int se_quad_i = GetSEQuadraticIndex(i);
                                    // Diagonal: (∂²L/∂f_k²) * f_i² * f_i²
                                    hess_this_type[se_quad_i * nparam + se_quad_i] += d2L_dfkdfk * fac_val[i] * fac_val[i] * fac_val[i] * fac_val[i];

                                    // Cross-terms with other SE quadratic
                                    for (int j = i + 1; j < n_input_factors; j++) {
                                        int se_quad_j = GetSEQuadraticIndex(j);
                                        hess_this_type[se_quad_i * nparam + se_quad_j] += d2L_dfkdfk * fac_val[i] * fac_val[i] * fac_val[j] * fac_val[j];
                                    }

                                    // SE quadratic × SE residual var
                                    hess_this_type[se_quad_i * nparam + se_res_idx] += d2L_dfkdfk * fac_val[i] * fac_val[i] * dfk_dres_var;
                                }
                            }

                            // SE residual var diagonal
                            hess_this_type[se_res_idx * nparam + se_res_idx] +=
                                d2L_dfkdfk * dfk_dres_var * dfk_dres_var + dL_dfk * d2fk_dres_var2;

                            // ===== Type-specific SE intercept Hessian contributions =====
                            // se_intercept_type_t enters the likelihood only for its own type
                            // (ityp == t), so we accumulate its Hessian terms only during
                            // that type's pass through the loop. Cross-terms between different
                            // se_intercept_type_t and se_intercept_type_{t'} (t ≠ t') are zero
                            // because the two parameters never appear in the same likelihood
                            // term — each type's contribution is additive and separable in its
                            // own intercept. ∂f_k/∂(se_intercept_type_t) = 1 exactly as for the
                            // base se_intercept, so the formulas mirror the se_intercept ones.
                            if (ntyp > 1 && ityp > 0 && nse_type_intercepts > 0) {
                                int typ_idx = GetSETypeInterceptIndex(ityp - 1);
                                // Layout invariants (Option A):
                                //   se_int_idx < se_lin_idx[j] < (se_quad_idx[j]) < typ_idx < se_res_idx
                                //   factor_var_i < typ_idx for all input-factor variances
                                //   typ_idx < se_cov params, factor_mean params, and model params

                                // Diagonal: (∂²L/∂f_k²) * 1 * 1
                                hess_this_type[typ_idx * nparam + typ_idx] += d2L_dfkdfk;

                                // typ_idx × se_intercept (se_int_idx < typ_idx)
                                hess_this_type[se_int_idx * nparam + typ_idx] += d2L_dfkdfk;

                                // typ_idx × se_linear_j (se_lin_idx < typ_idx)
                                for (int j = 0; j < n_input_factors; j++) {
                                    int se_lin_idx_j = GetSELinearIndex(j);
                                    hess_this_type[se_lin_idx_j * nparam + typ_idx] += d2L_dfkdfk * fac_val[j];
                                }

                                // typ_idx × se_quadratic_j (SE_QUADRATIC only; se_quad_idx < typ_idx)
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    for (int j = 0; j < n_input_factors; j++) {
                                        int se_quad_idx_j = GetSEQuadraticIndex(j);
                                        hess_this_type[se_quad_idx_j * nparam + typ_idx] += d2L_dfkdfk * fac_val[j] * fac_val[j];
                                    }
                                }

                                // typ_idx × se_residual_var (typ_idx < se_res_idx)
                                hess_this_type[typ_idx * nparam + se_res_idx] += d2L_dfkdfk * dfk_dres_var;

                                // typ_idx × factor_var_i (factor_var_i < typ_idx)
                                // Mirrors the factor_var × se_intercept cross-term (∂f_k/∂se_int_type = 1).
                                for (int i = 0; i < n_input_factors; i++) {
                                    double d2L_dfi_dfk = modHess[i * nDimModHess + (nfac - 1)];
                                    double dfk_dvar_i = dfk_dfj[i] * df_dsigma2[i];
                                    hess_this_type[i * nparam + typ_idx] +=
                                        d2L_dfi_dfk * df_dsigma2[i] + d2L_dfkdfk * dfk_dvar_i;
                                }

                                // typ_idx × SE covariate params (typ_idx < se_cov_idx)
                                if (use_se_covariates) {
                                    for (int m = 0; m < n_se_covariates; m++) {
                                        int se_cov_idx = se_covariate_param_start + m;
                                        double X_m = se_covariate_data[iobs][m];
                                        hess_this_type[typ_idx * nparam + se_cov_idx] += d2L_dfkdfk * X_m;
                                    }
                                }

                                // typ_idx × model params (typ_idx < firstpar + iparam)
                                for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                    int param_idx = firstpar + iparam;
                                    double d2L_dfk_dparam = modHess[(nfac - 1) * nDimModHess + (nfac + iparam)];
                                    hess_this_type[typ_idx * nparam + param_idx] += d2L_dfk_dparam;
                                }
                            }

                            // ===== SE covariate terms =====
                            if (use_se_covariates) {
                                // SE intercept × SE covariate cross-terms
                                for (int m = 0; m < n_se_covariates; m++) {
                                    int se_cov_idx = se_covariate_param_start + m;
                                    double X_m = se_covariate_data[iobs][m];
                                    // ∂²L/∂(se_intercept)∂β_m = (∂²L/∂f_k²) * 1 * X_m
                                    hess_this_type[se_int_idx * nparam + se_cov_idx] += d2L_dfkdfk * X_m;
                                }

                                // SE linear × SE covariate cross-terms
                                for (int j = 0; j < n_input_factors; j++) {
                                    int se_lin_idx = GetSELinearIndex(j);
                                    for (int m = 0; m < n_se_covariates; m++) {
                                        int se_cov_idx = se_covariate_param_start + m;
                                        double X_m = se_covariate_data[iobs][m];
                                        // ∂²L/∂(se_linear_j)∂β_m = (∂²L/∂f_k²) * f_j * X_m
                                        hess_this_type[se_lin_idx * nparam + se_cov_idx] += d2L_dfkdfk * fac_val[j] * X_m;
                                    }
                                }

                                // SE_QUADRATIC: SE quadratic × SE covariate cross-terms
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    for (int j = 0; j < n_input_factors; j++) {
                                        int se_quad_idx = GetSEQuadraticIndex(j);
                                        for (int m = 0; m < n_se_covariates; m++) {
                                            int se_cov_idx = se_covariate_param_start + m;
                                            double X_m = se_covariate_data[iobs][m];
                                            // ∂²L/∂(se_quadratic_j)∂β_m = (∂²L/∂f_k²) * f_j² * X_m
                                            hess_this_type[se_quad_idx * nparam + se_cov_idx] += d2L_dfkdfk * fac_val[j] * fac_val[j] * X_m;
                                        }
                                    }
                                }

                                // SE residual var × SE covariate cross-terms
                                for (int m = 0; m < n_se_covariates; m++) {
                                    int se_cov_idx = se_covariate_param_start + m;
                                    double X_m = se_covariate_data[iobs][m];
                                    // ∂²L/∂(se_residual_var)∂β_m = (∂²L/∂f_k²) * (∂f_k/∂σ²_ε) * X_m
                                    hess_this_type[se_res_idx * nparam + se_cov_idx] += d2L_dfkdfk * dfk_dres_var * X_m;
                                }

                                // SE covariate × SE covariate terms
                                for (int m = 0; m < n_se_covariates; m++) {
                                    int se_cov_m = se_covariate_param_start + m;
                                    double X_m = se_covariate_data[iobs][m];

                                    // Diagonal: (∂²L/∂f_k²) * X_m * X_m
                                    hess_this_type[se_cov_m * nparam + se_cov_m] += d2L_dfkdfk * X_m * X_m;

                                    // Cross-terms with other SE covariates
                                    for (int n = m + 1; n < n_se_covariates; n++) {
                                        int se_cov_n = se_covariate_param_start + n;
                                        double X_n = se_covariate_data[iobs][n];
                                        hess_this_type[se_cov_m * nparam + se_cov_n] += d2L_dfkdfk * X_m * X_n;
                                    }
                                }
                            }

                            // ===== Cross-terms with model parameters =====
                            for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                int param_idx = firstpar + iparam;

                                // Factor variance × model param
                                for (int i = 0; i < n_input_factors; i++) {
                                    double d2L_dfi_dparam = modHess[i * nDimModHess + (nfac + iparam)];
                                    double d2L_dfk_dparam = modHess[(nfac - 1) * nDimModHess + (nfac + iparam)];
                                    double dfk_dvar_i = dfk_dfj[i] * df_dsigma2[i];

                                    hess_this_type[i * nparam + param_idx] +=
                                        d2L_dfi_dparam * df_dsigma2[i] + d2L_dfk_dparam * dfk_dvar_i;
                                }

                                // SE params × model param
                                double d2L_dfk_dparam = modHess[(nfac - 1) * nDimModHess + (nfac + iparam)];

                                hess_this_type[se_int_idx * nparam + param_idx] += d2L_dfk_dparam;

                                for (int j = 0; j < n_input_factors; j++) {
                                    int se_lin_idx = GetSELinearIndex(j);
                                    hess_this_type[se_lin_idx * nparam + param_idx] += d2L_dfk_dparam * fac_val[j];
                                }

                                // SE_QUADRATIC: SE quadratic × model param
                                if (factor_structure == FactorStructure::SE_QUADRATIC) {
                                    for (int j = 0; j < n_input_factors; j++) {
                                        int se_quad_idx = GetSEQuadraticIndex(j);
                                        hess_this_type[se_quad_idx * nparam + param_idx] += d2L_dfk_dparam * fac_val[j] * fac_val[j];
                                    }
                                }

                                hess_this_type[se_res_idx * nparam + param_idx] += d2L_dfk_dparam * dfk_dres_var;

                                // SE covariate × model param
                                if (use_se_covariates) {
                                    for (int m = 0; m < n_se_covariates; m++) {
                                        int se_cov_idx = se_covariate_param_start + m;
                                        double X_m = se_covariate_data[iobs][m];
                                        // ∂²L/∂β_m∂(model_param) = (∂²L/∂f_k∂param) * X_m
                                        hess_this_type[se_cov_idx * nparam + param_idx] += d2L_dfk_dparam * X_m;
                                    }
                                }
                            }

                            // ===== Model param × model param =====
                            for (int i = 0; i < base_param_count_hess; i++) {
                                for (int j = i; j < base_param_count_hess; j++) {
                                    int modhess_idx = (nfac + i) * nDimModHess + (nfac + j);
                                    int full_i = firstpar + i;
                                    int full_j = firstpar + j;
                                    hess_this_type[full_i * nparam + full_j] += modHess[modhess_idx];
                                }
                            }

                        } else if (fac_corr && nfac == 2) {
                            // ===== Correlated 2-factor Hessian =====
                            // Use pre-computed sigma_fac and x_node_fac (z1, z2)
                            double df1_dsigma1sq = x_node_fac[0] / (2.0 * sigma_fac[0]);
                            double df2_dsigma2sq = (rho * x_node_fac[0] + sqrt_1_minus_rho2 * x_node_fac[1]) / (2.0 * sigma_fac[1]);
                            double df2_drho = sigma_fac[1] * (x_node_fac[0] - rho * x_node_fac[1] / sqrt_1_minus_rho2);
                            double d2f2_drho2 = -sigma_fac[1] * x_node_fac[1] / (sqrt_1_minus_rho2 * sqrt_1_minus_rho2 * sqrt_1_minus_rho2);
                            double d2f2_drho_dsigma2sq = (x_node_fac[0] - rho * x_node_fac[1] / sqrt_1_minus_rho2) / (2.0 * sigma_fac[1]);

                            double d2L_df1df1 = modHess[0 * nDimModHess + 0];
                            double d2L_df1df2 = modHess[0 * nDimModHess + 1];
                            double d2L_df2df2 = modHess[1 * nDimModHess + 1];
                            double dL_df1 = modEval[1];
                            double dL_df2 = modEval[2];

                            double sigma1_cubed = sigma_fac[0] * sigma_fac[0] * sigma_fac[0];
                            double sigma2_cubed = sigma_fac[1] * sigma_fac[1] * sigma_fac[1];
                            double d2f1_dsigma1sq_sq = -x_node_fac[0] / (4.0 * sigma1_cubed);
                            hess_this_type[0 * nparam + 0] += d2L_df1df1 * df1_dsigma1sq * df1_dsigma1sq
                                                              + dL_df1 * d2f1_dsigma1sq_sq;
                            hess_this_type[0 * nparam + 1] += d2L_df1df2 * df1_dsigma1sq * df2_dsigma2sq;
                            hess_this_type[0 * nparam + 2] += d2L_df1df2 * df1_dsigma1sq * df2_drho;

                            double d2f2_dsigma2sq_sq = -(rho * x_node_fac[0] + sqrt_1_minus_rho2 * x_node_fac[1]) / (4.0 * sigma2_cubed);
                            hess_this_type[1 * nparam + 1] += d2L_df2df2 * df2_dsigma2sq * df2_dsigma2sq
                                                              + dL_df2 * d2f2_dsigma2sq_sq;
                            hess_this_type[1 * nparam + 2] += d2L_df2df2 * df2_dsigma2sq * df2_drho
                                                              + dL_df2 * d2f2_drho_dsigma2sq;
                            hess_this_type[2 * nparam + 2] += d2L_df2df2 * df2_drho * df2_drho
                                                              + dL_df2 * d2f2_drho2;

                            // Use base_param_count_hess for loop bounds since modHess
                            // only contains base parameters, not type-specific intercepts
                            for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                int param_idx = firstpar + iparam;
                                double d2L_df1_dparam = modHess[0 * nDimModHess + (nfac + iparam)];
                                double d2L_df2_dparam = modHess[1 * nDimModHess + (nfac + iparam)];
                                hess_this_type[0 * nparam + param_idx] += d2L_df1_dparam * df1_dsigma1sq;
                                hess_this_type[1 * nparam + param_idx] += d2L_df2_dparam * df2_dsigma2sq;
                                hess_this_type[2 * nparam + param_idx] += d2L_df2_dparam * df2_drho;
                            }

                            for (int i = 0; i < base_param_count_hess; i++) {
                                for (int j = i; j < base_param_count_hess; j++) {
                                    int modhess_idx = (nfac + i) * nDimModHess + (nfac + j);
                                    int full_i = firstpar + i;
                                    int full_j = firstpar + j;
                                    hess_this_type[full_i * nparam + full_j] += modHess[modhess_idx];
                                }
                            }
                        } else {
                            // ===== Independent factors Hessian =====
                            // OPTIMIZATION: Separate loops to avoid conditionals in inner loops
                            // This matches the legacy TMinLkhd.cc structure for better performance

                            // 1. Factor × Factor terms (with chain rule and diagonal second derivative)
                            // For mixture models: variance indices are mixture-specific
                            for (int i = 0; i < nfac; i++) {
                                // Get the correct parameter index for this factor's variance
                                // For nmix=1: var_idx_i = i
                                // For nmix>1: var_idx_i = GetFactorVarianceIndex(imix, i)
                                int var_idx_i = (nmix > 1 && i < n_factors_for_mixture)
                                              ? GetFactorVarianceIndex(imix, i) : i;

                                for (int j = i; j < nfac; j++) {
                                    int var_idx_j = (nmix > 1 && j < n_factors_for_mixture)
                                                  ? GetFactorVarianceIndex(imix, j) : j;

                                    int modhess_idx = i * nDimModHess + j;
                                    int full_idx = var_idx_i * nparam + var_idx_j;
                                    double chain_factor = df_dsigma2[i] * df_dsigma2[j];
                                    hess_this_type[full_idx] += modHess[modhess_idx] * chain_factor;
                                    // Diagonal second derivative term
                                    if (i == j) {
                                        hess_this_type[full_idx] += modEval[i + 1] * d2f_dsigma2_sq[i];
                                    }

                                    // For mixture Hessians: accumulate factor-factor log-likelihood Hessian
                                    // directly from modHess (without df/dsigma² factors)
                                    if (nmix > 1 && ntyp == 1 && i < n_factors_for_mixture && j < n_factors_for_mixture) {
                                        int ff_idx = i * n_factors_for_mixture + j;
                                        d2logL_dfdf[ff_idx] += modHess[modhess_idx];
                                        if (i != j) {
                                            // Symmetric entry
                                            int ff_idx_sym = j * n_factors_for_mixture + i;
                                            d2logL_dfdf[ff_idx_sym] += modHess[modhess_idx];
                                        }
                                    }
                                }
                            }

                            // For mixture Hessians: accumulate factor gradient dlogL/df_k from modEval
                            if (nmix > 1 && ntyp == 1) {
                                for (int k = 0; k < n_factors_for_mixture; k++) {
                                    dlogL_df[k] += modEval[k + 1];  // modEval[k+1] = d(log L_model)/df_k
                                }
                            }

                            // 2. Factor × Model param terms (with chain rule on factor side only)
                            for (int i = 0; i < nfac; i++) {
                                // Use mixture-specific variance index
                                int var_idx_i = (nmix > 1 && i < n_factors_for_mixture)
                                              ? GetFactorVarianceIndex(imix, i) : i;

                                for (int jparam = 0; jparam < base_param_count_hess; jparam++) {
                                    int j = nfac + jparam;
                                    int modhess_idx = i * nDimModHess + j;
                                    int full_idx = var_idx_i * nparam + (firstpar + jparam);
                                    hess_this_type[full_idx] += modHess[modhess_idx] * df_dsigma2[i];
                                }
                            }

                            // For mixture Hessians: accumulate d²logL/df_k d(beta_j) for RAW conversion later
                            // Used for mean × beta and weight × beta cross-Hessians in STEP 11e
                            if (nmix > 1 && ntyp == 1 && factor_structure == FactorStructure::INDEPENDENT) {
                                for (int k = 0; k < n_factors_for_mixture; k++) {
                                    for (int jparam = 0; jparam < base_param_count_hess; jparam++) {
                                        int modhess_idx = k * nDimModHess + (nfac + jparam);
                                        int storage_idx = k * n_model_params + (firstpar + jparam - first_model_param);
                                        d2logL_df_dbeta_thisintpt[storage_idx] += modHess[modhess_idx];
                                    }
                                }
                            }

                            // 3. Model param × Model param terms (no chain rule - direct copy)
                            for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                for (int jparam = iparam; jparam < base_param_count_hess; jparam++) {
                                    int modhess_idx = (nfac + iparam) * nDimModHess + (nfac + jparam);
                                    int full_idx = (firstpar + iparam) * nparam + (firstpar + jparam);
                                    hess_this_type[full_idx] += modHess[modhess_idx];
                                }
                            }

                            // 4. Factor mean covariate Hessian terms
                            if (use_factor_mean_covariates) {
                                // Factor mean × Factor mean terms
                                // d²L/d(gamma_{k1,j1})d(gamma_{k2,j2}) = d²L/df_{k1}df_{k2} * X_{j1} * X_{j2}
                                for (int k1 = 0; k1 < n_factors_with_mean; k1++) {
                                    for (int j1 = 0; j1 < n_factor_mean_covariates; j1++) {
                                        int param_idx1 = factor_mean_param_start + k1 * n_factor_mean_covariates + j1;
                                        double x_j1 = factor_mean_covariate_data[iobs][j1];

                                        for (int k2 = k1; k2 < n_factors_with_mean; k2++) {
                                            int j2_start = (k1 == k2) ? j1 : 0;
                                            for (int j2 = j2_start; j2 < n_factor_mean_covariates; j2++) {
                                                int param_idx2 = factor_mean_param_start + k2 * n_factor_mean_covariates + j2;
                                                double x_j2 = factor_mean_covariate_data[iobs][j2];

                                                int modhess_idx = k1 * nDimModHess + k2;
                                                int full_idx = param_idx1 * nparam + param_idx2;
                                                hess_this_type[full_idx] += modHess[modhess_idx] * x_j1 * x_j2;
                                            }
                                        }
                                    }
                                }

                                // Factor variance × Factor mean cross terms
                                // d²L/d(sigma²_k1)d(gamma_{k2,j}) = d²L/df_{k1}df_{k2} * df/dsigma² * X_j
                                for (int k1 = 0; k1 < nfac; k1++) {
                                    for (int k2 = 0; k2 < n_factors_with_mean; k2++) {
                                        for (int j = 0; j < n_factor_mean_covariates; j++) {
                                            int param_idx_mean = factor_mean_param_start + k2 * n_factor_mean_covariates + j;
                                            double x_j = factor_mean_covariate_data[iobs][j];

                                            int modhess_idx = k1 * nDimModHess + k2;
                                            int full_idx = k1 * nparam + param_idx_mean;
                                            hess_this_type[full_idx] += modHess[modhess_idx] * df_dsigma2[k1] * x_j;
                                        }
                                    }
                                }

                                // Factor mean × Model param cross terms
                                // d²L/d(gamma_{k,j})d(theta_m) = d²L/df_k d(theta_m) * X_j
                                for (int k = 0; k < n_factors_with_mean; k++) {
                                    for (int j = 0; j < n_factor_mean_covariates; j++) {
                                        int param_idx_mean = factor_mean_param_start + k * n_factor_mean_covariates + j;
                                        double x_j = factor_mean_covariate_data[iobs][j];

                                        for (int mparam = 0; mparam < base_param_count_hess; mparam++) {
                                            int modhess_idx = k * nDimModHess + (nfac + mparam);
                                            int full_idx = param_idx_mean * nparam + (firstpar + mparam);
                                            hess_this_type[full_idx] += modHess[modhess_idx] * x_j;
                                        }
                                    }
                                }
                            }
                        }

                        // ===== STEP 9b: Add Hessian terms for type-specific intercepts =====
                        // Type-specific intercepts have the same derivative structure as base intercept.
                        // The Hessian contribution mirrors that of the base intercept (index 0 in model params).
                        if (ntyp > 1 && ityp > 0 && model_uses_types_hess) {
                            int type_intercept_idx = firstpar + base_param_count_hess + (ityp - 1);

                            // Diagonal term: d²L/d(type_intercept)² = d²L/d(intercept)²
                            // In modHess, intercept is at position (nfac, nfac) -> index nfac * nDimModHess + nfac
                            double hess_intercept_diag = modHess[nfac * nDimModHess + nfac];
                            hess_this_type[type_intercept_idx * nparam + type_intercept_idx] += hess_intercept_diag;

                            // Cross-terms with factor variances: d²L/d(type_intercept)d(factor_var_k)
                            for (int k = 0; k < nfac; k++) {
                                // modHess[k * nDimModHess + nfac] = d²L/df_k d(intercept)
                                double hess_fk_intercept = modHess[k * nDimModHess + nfac];

                                // Use pre-computed df_dsigma2[k] = x_node / (2*sigma)
                                // Hessian cross-term: d²L/d(factor_var_k)d(type_intercept)
                                int full_idx = k * nparam + type_intercept_idx;
                                hess_this_type[full_idx] += hess_fk_intercept * df_dsigma2[k];
                            }

                            // Cross-terms with base model parameters
                            for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                int param_idx = firstpar + iparam;
                                // modHess[(nfac + iparam) * nDimModHess + nfac] = d²L/d(param_iparam)d(intercept)
                                // Note: Hessian is symmetric, upper triangular stored
                                int mh_row = nfac;  // intercept row
                                int mh_col = nfac + iparam;  // param column
                                if (mh_row > mh_col) std::swap(mh_row, mh_col);
                                double hess_param_intercept = modHess[mh_row * nDimModHess + mh_col];

                                // Add to cross-term (type_intercept, param)
                                int full_i = std::min(type_intercept_idx, param_idx);
                                int full_j = std::max(type_intercept_idx, param_idx);
                                hess_this_type[full_i * nparam + full_j] += hess_param_intercept;
                            }
                        }
                    }
                } // End of model loop for this type

                // ===== Weight by type probability and accumulate =====
                // OPTIMIZATION: For ntyp=1, skip type_weighted_* and accumulate directly into total* later
                if (ntyp > 1) {
                    type_weighted_prob += type_prob * prob_this_type;
                }

                if (iflag >= 2) {
                    // Gradient of type-weighted likelihood:
                    // d(π_t * L_t)/dθ = (dπ_t/dθ) * L_t + π_t * (dL_t/dθ)
                    // For model parameters: dπ_t/dθ = 0, so just π_t * (dL_t/dθ)
                    // For ntyp=1, skip - we use grad_this_type directly later.
                    // Include tied params under equality constraints: the
                    // downstream aggregation (mixgrad += type_weighted_grad)
                    // uses gradparlist, so tied positions of type_weighted_grad
                    // must be populated here or their contribution is lost.
                    if (ntyp > 1) {
                        for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                            int ipar = gradparlist[fi];
                            type_weighted_grad[ipar] += type_prob * prob_this_type * grad_this_type[ipar];
                        }
                    }

                    // Gradient w.r.t. type model parameters (for types > 1)
                    // η_t = typeprob_intercept_t + Σ_k λ_{t,k} * f_k
                    // π_t = exp(η_t) / (1 + Σ_s exp(η_s))
                    // dπ_ityp/dη_{t+1} = π_ityp * (δ_{ityp,t+1} - π_{t+1})
                    if (ntyp > 1) {
                        // Gradient w.r.t. typeprob intercepts
                        // dπ_ityp/d(intercept_t) = π_ityp * (δ_{ityp,t+1} - π_{t+1})
                        for (int t = 0; t < ntyp - 1; t++) {
                            int intercept_idx = GetTypeProbInterceptIndex(t);
                            double dpi_dintercept = type_probs[ityp] *
                                ((ityp == t + 1 ? 1.0 : 0.0) - type_probs[t + 1]);
                            type_weighted_grad[intercept_idx] += dpi_dintercept * prob_this_type;
                        }

                        // Gradient w.r.t. type loadings
                        // dπ_ityp/dλ_{s,k} = π_ityp * (δ_{ityp,s+1} - π_{s+1}) * f_k
                        for (int t = 0; t < ntyp - 1; t++) {  // Types 2, 3, ... (t+1 in 1-indexed)
                            for (int k = 0; k < nfac; k++) {
                                int param_idx = GetTypeLoadingIndex(t, k);
                                // ityp is 0-based (0 = type 1, 1 = type 2, etc.)
                                // t is 0-based index for types with loadings (0 = type 2, etc.)
                                double dpi_dlambda = type_probs[ityp] *
                                    ((ityp == t + 1 ? 1.0 : 0.0) - type_probs[t + 1]) * fac_val[k];
                                type_weighted_grad[param_idx] += dpi_dlambda * prob_this_type;
                            }
                        }

                        // Additional gradient for factor variance through type probabilities
                        // dπ_t/dσ²_k = dπ_t/df_k * df_k/dσ²_k
                        // where df_k/dσ²_k = x_q / (2 * σ_k) since f = σ * x and dσ/dσ² = 1/(2σ)
                        // and dπ_t/df_k = Σ_s π_t * (δ_{t,s} - π_s) * λ_{s,k}
                        for (int k = 0; k < nfac; k++) {
                            // Compute dπ_t/df_k = Σ_s π_t * (δ_{t,s} - π_s) * λ_{s,k}
                            double dpi_df_k = 0.0;
                            for (int s = 0; s < ntyp - 1; s++) {
                                int lambda_idx = GetTypeLoadingIndex(s, k);
                                double lambda_sk = param[lambda_idx];
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                dpi_df_k += type_probs[ityp] * (delta_ts - type_probs[s + 1]) * lambda_sk;
                            }

                            // Add (dπ_t/dσ²_k) * L_t to factor variance gradient
                            // Use pre-computed df_dsigma2[k]
                            type_weighted_grad[k] += dpi_df_k * df_dsigma2[k] * prob_this_type;
                        }
                    }
                }

                if (iflag == 3) {
                    // ===== STEP 9c: Hessian for model parameters (original formula) =====
                    // d²L_t/dθdφ = L_t * [d²(log L_t)/dθdφ + d(log L_t)/dθ * d(log L_t)/dφ]
                    // Weighted by π_t: contribution to d²L_mix/dθdφ is π_t * d²L_t/dθdφ
                    // For ntyp=1, skip - we use hess_this_type directly later.
                    // Include tied params under equality constraints (see comment
                    // at the outer totalhess zeroing loop).
                    if (ntyp > 1) {
                        for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                            int i = gradparlist[fi];
                            for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                                int j = gradparlist[fj];
                                int idx = i * nparam + j;
                                // Accumulate Hessian weighted by type probability
                                type_weighted_hess[idx] += type_prob * prob_this_type *
                                    (hess_this_type[idx] + grad_this_type[i] * grad_this_type[j]);
                            }
                        }
                    }

                    // ===== STEP 9d: Additional Hessian terms for type model parameters =====
                    // The mixture L_mix = Σ_t π_t * L_t depends on type parameters through π_t
                    // d²L_mix/dθdφ = Σ_t [d²π_t/dθdφ * L_t + dπ_t/dθ * dL_t/dφ + dπ_t/dφ * dL_t/dθ + π_t * d²L_t/dθdφ]
                    // The last term is already handled above.
                    // For typeprob intercept α_s: dπ_t/dα_s = π_t * (δ_{t,s+1} - π_{s+1})
                    // For type loading λ_{s,k}: dπ_t/dλ = π_t * (δ_{t,s+1} - π_{s+1}) * f_k
                    // For model parameters: dπ_t/dθ = 0
                    if (ntyp > 1) {
                        // 0. Second derivative terms for typeprob intercepts: d²π_t/dα_s dα_r * L_t
                        // d²π_t/dη_s dη_r = π_t * (δ_{t,s} - π_s) * (δ_{t,r} - π_r) - π_t * π_s * (δ_{s,r} - π_r)
                        for (int s = 0; s < ntyp - 1; s++) {
                            for (int r = s; r < ntyp - 1; r++) {  // r >= s for upper triangle
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                double delta_tr = (ityp == r + 1) ? 1.0 : 0.0;
                                double delta_sr = (s == r) ? 1.0 : 0.0;

                                double d2pi_deta = type_probs[ityp] * (delta_ts - type_probs[s + 1])
                                                                    * (delta_tr - type_probs[r + 1])
                                                 - type_probs[ityp] * type_probs[s + 1]
                                                                    * (delta_sr - type_probs[r + 1]);

                                int intercept_idx_s = GetTypeProbInterceptIndex(s);
                                int intercept_idx_r = GetTypeProbInterceptIndex(r);
                                int pi = std::min(intercept_idx_s, intercept_idx_r);
                                int pj = std::max(intercept_idx_s, intercept_idx_r);
                                int idx = pi * nparam + pj;

                                // d²π_t/dα_s dα_r = d2pi_deta (no f_k factor for intercepts)
                                type_weighted_hess[idx] += d2pi_deta * prob_this_type;
                            }
                        }

                        // 0b. Cross terms: d²π_t/dα_s dλ_{r,l} * L_t (intercept × loading)
                        for (int s = 0; s < ntyp - 1; s++) {
                            int intercept_idx_s = GetTypeProbInterceptIndex(s);
                            for (int r = 0; r < ntyp - 1; r++) {
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                double delta_tr = (ityp == r + 1) ? 1.0 : 0.0;
                                double delta_sr = (s == r) ? 1.0 : 0.0;

                                double d2pi_deta = type_probs[ityp] * (delta_ts - type_probs[s + 1])
                                                                    * (delta_tr - type_probs[r + 1])
                                                 - type_probs[ityp] * type_probs[s + 1]
                                                                    * (delta_sr - type_probs[r + 1]);

                                for (int l = 0; l < nfac; l++) {
                                    int loading_idx_rl = GetTypeLoadingIndex(r, l);
                                    int pi = std::min(intercept_idx_s, loading_idx_rl);
                                    int pj = std::max(intercept_idx_s, loading_idx_rl);
                                    int idx = pi * nparam + pj;

                                    // d²π_t/dα_s dλ_{r,l} = d2pi_deta * f_l
                                    type_weighted_hess[idx] += d2pi_deta * fac_val[l] * prob_this_type;
                                }
                            }
                        }

                        // 0c. Cross terms: dπ_t/dα_s * dL_t/dθ (intercept × model_param)
                        for (int s = 0; s < ntyp - 1; s++) {
                            int intercept_idx = GetTypeProbInterceptIndex(s);
                            double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                            double dpi_dintercept = type_probs[ityp] * (delta_ts - type_probs[s + 1]);

                            for (int ipar = 0; ipar < nparam; ipar++) {
                                // Skip if ipar is a type parameter (intercept or loading)
                                if (ntyp_param > 0 && ipar >= type_param_start &&
                                    ipar < type_param_start + ntyp_param) {
                                    continue;
                                }

                                int pi = std::min(intercept_idx, ipar);
                                int pj = std::max(intercept_idx, ipar);
                                int idx = pi * nparam + pj;

                                // dπ_t/dα_s * dL_t/dθ = dπ_t/dα_s * L_t * d(log L_t)/dθ
                                type_weighted_hess[idx] += dpi_dintercept * prob_this_type * grad_this_type[ipar];
                            }
                        }

                        // 1. Second derivative terms: d²π_t/dλdλ * L_t
                        // d²π_t/dη_s dη_r = π_t * (δ_{t,s} - π_s) * (δ_{t,r} - π_r) - π_t * π_s * (δ_{s,r} - π_r)
                        for (int s = 0; s < ntyp - 1; s++) {
                            for (int r = s; r < ntyp - 1; r++) {  // r >= s for upper triangle
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                double delta_tr = (ityp == r + 1) ? 1.0 : 0.0;
                                double delta_sr = (s == r) ? 1.0 : 0.0;

                                double d2pi_deta = type_probs[ityp] * (delta_ts - type_probs[s + 1])
                                                                    * (delta_tr - type_probs[r + 1])
                                                 - type_probs[ityp] * type_probs[s + 1]
                                                                    * (delta_sr - type_probs[r + 1]);

                                for (int k = 0; k < nfac; k++) {
                                    for (int l = (s == r ? k : 0); l < nfac; l++) {
                                        int param_idx_sk = GetTypeLoadingIndex(s, k);
                                        int param_idx_rl = GetTypeLoadingIndex(r, l);

                                        int pi = std::min(param_idx_sk, param_idx_rl);
                                        int pj = std::max(param_idx_sk, param_idx_rl);
                                        int idx = pi * nparam + pj;

                                        double d2pi_dlambda = d2pi_deta * fac_val[k] * fac_val[l];
                                        type_weighted_hess[idx] += d2pi_dlambda * prob_this_type;
                                    }
                                }
                            }
                        }

                        // 2. Cross terms: dπ_t/dλ * dL_t/dθ (for λ × model_param)
                        // Note: dL_t/dθ = L_t * d(log L_t)/dθ = prob_this_type * grad_this_type
                        for (int s = 0; s < ntyp - 1; s++) {
                            for (int k = 0; k < nfac; k++) {
                                int type_loading_idx = GetTypeLoadingIndex(s, k);
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                double dpi_dlambda = type_probs[ityp] * (delta_ts - type_probs[s + 1]) * fac_val[k];

                                // Cross with model parameters (not type loadings)
                                for (int ipar = 0; ipar < nparam; ipar++) {
                                    // Skip if ipar is a type loading (handled in section 1)
                                    if (ntyp_param > 0 && ipar >= type_param_start &&
                                        ipar < type_param_start + ntyp_param) {
                                        continue;
                                    }

                                    int pi = std::min(type_loading_idx, ipar);
                                    int pj = std::max(type_loading_idx, ipar);
                                    int idx = pi * nparam + pj;

                                    // dπ_t/dλ * dL_t/dθ = dπ_t/dλ * L_t * d(log L_t)/dθ
                                    type_weighted_hess[idx] += dpi_dlambda * prob_this_type * grad_this_type[ipar];
                                }
                            }
                        }

                        // 3. Hessian terms for factor variance through type probabilities
                        // dπ_t/dσ²_k = dπ_t/df_k * df_k/dσ²_k
                        // where df/dσ² = x_q / (2σ) since f = σ * x and dσ/dσ² = 1/(2σ)
                        // Use pre-computed df_dsigma2[k]
                        for (int k = 0; k < nfac; k++) {
                            // Compute dπ_t/df_k = Σ_s π_t * (δ_{t,s} - π_s) * λ_{s,k}
                            double dpi_df_k = 0.0;
                            for (int s = 0; s < ntyp - 1; s++) {
                                int lambda_idx = GetTypeLoadingIndex(s, k);
                                double lambda_sk = param[lambda_idx];
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                dpi_df_k += type_probs[ityp] * (delta_ts - type_probs[s + 1]) * lambda_sk;
                            }
                            double dpi_dsigma2_k = dpi_df_k * df_dsigma2[k];

                            // 3a. Cross terms: dπ_t/dσ²_k * dL_t/dθ (for σ² × model_param)
                            for (int ipar = 0; ipar < nparam; ipar++) {
                                // Skip type loadings (handled separately below)
                                if (ntyp_param > 0 && ipar >= type_param_start &&
                                    ipar < type_param_start + ntyp_param) {
                                    continue;
                                }
                                // Skip factor variances (handled in 3b)
                                if (ipar < nfac) continue;

                                int pi = std::min(k, ipar);
                                int pj = std::max(k, ipar);
                                int idx = pi * nparam + pj;

                                // dπ_t/dσ² * dL_t/dθ
                                type_weighted_hess[idx] += dpi_dsigma2_k * prob_this_type * grad_this_type[ipar];
                            }

                            // 3b. Cross terms: dπ_t/dσ²_k * dL_t/dσ²_l (for σ²_k × σ²_l)
                            // Plus second derivative d²π_t/dσ²_k dσ²_l * L_t
                            // Use pre-computed df_dsigma2[l]
                            for (int l = k; l < nfac; l++) {
                                // Compute dπ_t/df_l
                                double dpi_df_l = 0.0;
                                for (int s = 0; s < ntyp - 1; s++) {
                                    int lambda_idx = GetTypeLoadingIndex(s, l);
                                    double lambda_sl = param[lambda_idx];
                                    double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                    dpi_df_l += type_probs[ityp] * (delta_ts - type_probs[s + 1]) * lambda_sl;
                                }
                                double dpi_dsigma2_l = dpi_df_l * df_dsigma2[l];

                                int idx = k * nparam + l;

                                // Cross term: dπ_t/dσ²_k * dL_t/dσ²_l + dπ_t/dσ²_l * dL_t/dσ²_k
                                // Note: dL_t/dσ² = L_t * d(log L_t)/dσ² = prob_this_type * grad_this_type[factor_idx]
                                if (k == l) {
                                    // Diagonal: both cross terms coincide to 2 * dπ_t/dσ² * dL_t/dσ²
                                    type_weighted_hess[idx] += 2.0 * dpi_dsigma2_k * prob_this_type * grad_this_type[l];
                                } else {
                                    // Off-diagonal: add both cross terms
                                    type_weighted_hess[idx] += dpi_dsigma2_k * prob_this_type * grad_this_type[l];
                                    type_weighted_hess[idx] += dpi_dsigma2_l * prob_this_type * grad_this_type[k];
                                }

                                // Second derivative: d²π_t/dσ²_k dσ²_l * L_t
                                // d²π_t/df_k df_l = Σ_s,r (d²π_t/dη_s dη_r) * λ_{s,k} * λ_{r,l}
                                double d2pi_df_k_df_l = 0.0;
                                for (int s = 0; s < ntyp - 1; s++) {
                                    for (int r = 0; r < ntyp - 1; r++) {
                                        int lambda_sk_idx = GetTypeLoadingIndex(s, k);
                                        int lambda_rl_idx = GetTypeLoadingIndex(r, l);
                                        double lambda_sk = param[lambda_sk_idx];
                                        double lambda_rl = param[lambda_rl_idx];

                                        double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                        double delta_tr = (ityp == r + 1) ? 1.0 : 0.0;
                                        double delta_sr = (s == r) ? 1.0 : 0.0;

                                        double d2pi_deta_s_r = type_probs[ityp] * (delta_ts - type_probs[s + 1])
                                                                                * (delta_tr - type_probs[r + 1])
                                                             - type_probs[ityp] * type_probs[s + 1]
                                                                                * (delta_sr - type_probs[r + 1]);
                                        d2pi_df_k_df_l += d2pi_deta_s_r * lambda_sk * lambda_rl;
                                    }
                                }
                                double d2pi_dsigma2_k_l = d2pi_df_k_df_l * df_dsigma2[k] * df_dsigma2[l];
                                // Diagonal (k == l): also add dπ_t/df_k * d²f_k/d(σ²_k)²
                                // The chain rule for d²π_t/d(σ²_k)² is
                                //   d²π/df_k² * (df_k/dσ²_k)² + dπ/df_k * d²f_k/d(σ²_k)²
                                // The first term is already in d2pi_dsigma2_k_l; the second is only
                                // non-zero when k == l because f_k's GH scaling depends on σ²_k alone.
                                if (k == l) {
                                    // dpi_df_k was computed above as Σ_s π_t (δ_{t,s+1} - π_{s+1}) λ_{s,k}
                                    d2pi_dsigma2_k_l += dpi_df_k * d2f_dsigma2_sq[k];
                                }
                                type_weighted_hess[idx] += d2pi_dsigma2_k_l * prob_this_type;
                            }

                            // 3c. Cross terms: dπ_t/dσ²_k * dL_t/dλ_{s,l} + d²π_t/dσ²_k dλ_{s,l} * L_t
                            // (for σ²_k × λ_{s,l})
                            for (int s = 0; s < ntyp - 1; s++) {
                                for (int l = 0; l < nfac; l++) {
                                    int type_loading_idx = GetTypeLoadingIndex(s, l);

                                    // Note: dL_t/dλ = 0 (type-specific likelihoods don't depend on type loadings)
                                    // So we only need d²π_t/dσ²_k dλ_{s,l} * L_t

                                    // d²π_t/df_k dλ_{s,l} = d²π_t/df_k dη_s * f_l + dπ_t/dη_s * δ_{k,l} * df/dσ² * dσ²/dσ²
                                    //                     = (d²π_t/dη dη) * λ_{?,k} * f_l + ... (complex)
                                    // Simpler: d(dπ_t/dλ_{s,l})/dσ²_k = d(π_t * (δ_{t,s} - π_s) * f_l)/dσ²_k
                                    //        = (dπ_t/dσ²_k) * (δ_{t,s} - π_s) * f_l
                                    //          + π_t * (-dπ_s/dσ²_k) * f_l
                                    //          + π_t * (δ_{t,s} - π_s) * df_l/dσ²_k

                                    double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                    double delta_kl = (k == l) ? 1.0 : 0.0;

                                    // df_l/dσ²_k = 0 if k != l, else x_node_k / (2 * σ_k)
                                    double df_l_dsigma2_k = delta_kl * df_dsigma2[k];

                                    // dπ_s/dσ²_k (for the type s+1)
                                    double dpi_s_df_k = 0.0;
                                    for (int r = 0; r < ntyp - 1; r++) {
                                        int lambda_rk_idx = GetTypeLoadingIndex(r, k);
                                        double lambda_rk = param[lambda_rk_idx];
                                        double delta_s1_r1 = (s == r) ? 1.0 : 0.0;
                                        dpi_s_df_k += type_probs[s + 1] * (delta_s1_r1 - type_probs[r + 1]) * lambda_rk;
                                    }
                                    double dpi_s_dsigma2_k = dpi_s_df_k * df_dsigma2[k];

                                    // d²(π_t*(δ_{t,s}-π_s)*f_l)/dσ²_k
                                    double d2_term = dpi_dsigma2_k * (delta_ts - type_probs[s + 1]) * fac_val[l]
                                                   - type_probs[ityp] * dpi_s_dsigma2_k * fac_val[l]
                                                   + type_probs[ityp] * (delta_ts - type_probs[s + 1]) * df_l_dsigma2_k;

                                    int pi = std::min(k, type_loading_idx);
                                    int pj = std::max(k, type_loading_idx);
                                    int idx = pi * nparam + pj;

                                    type_weighted_hess[idx] += d2_term * prob_this_type;
                                }
                            }

                            // 3d. Cross terms: dπ_t/dσ²_k * dL_t/dα_s + d²π_t/dσ²_k dα_s * L_t
                            // (for σ²_k × typeprob_intercept_s)
                            // Note: dL_t/dα = 0 (type-specific likelihoods don't depend on typeprob intercepts)
                            // So we only need d²π_t/dσ²_k dα_s * L_t
                            for (int s = 0; s < ntyp - 1; s++) {
                                int intercept_idx = GetTypeProbInterceptIndex(s);

                                // d(dπ_t/dα_s)/dσ²_k = d(π_t * (δ_{t,s+1} - π_{s+1}))/dσ²_k
                                //                    = (dπ_t/dσ²_k) * (δ_{t,s+1} - π_{s+1})
                                //                      - π_t * (dπ_{s+1}/dσ²_k)

                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;

                                // dπ_{s+1}/dσ²_k (probability of type s+1, not ityp)
                                double dpi_s_df_k = 0.0;
                                for (int r = 0; r < ntyp - 1; r++) {
                                    int lambda_rk_idx = GetTypeLoadingIndex(r, k);
                                    double lambda_rk = param[lambda_rk_idx];
                                    double delta_s1_r1 = (s == r) ? 1.0 : 0.0;
                                    dpi_s_df_k += type_probs[s + 1] * (delta_s1_r1 - type_probs[r + 1]) * lambda_rk;
                                }
                                double dpi_s_dsigma2_k = dpi_s_df_k * df_dsigma2[k];

                                // d²(π_t*(δ_{t,s+1}-π_{s+1}))/dσ²_k
                                double d2_term = dpi_dsigma2_k * (delta_ts - type_probs[s + 1])
                                               - type_probs[ityp] * dpi_s_dsigma2_k;

                                int pi = std::min(k, intercept_idx);
                                int pj = std::max(k, intercept_idx);
                                int idx = pi * nparam + pj;

                                type_weighted_hess[idx] += d2_term * prob_this_type;
                            }
                        }
                    }
                }
            } // End of type loop

            // ===== STEP 10: Accumulate weighted contributions =====
            // OPTIMIZATION: Accumulate raw likelihood derivatives directly, convert once at end.
            // This eliminates the O(nparam²) hessilk conversion loop per integration point.
            //
            // Mathematical justification:
            // We need: totalhess = Σ w_q * d²L/dθ²  (raw Hessian, not log-likelihood)
            //          totalgrad = Σ w_q * dL/dθ    (raw gradient)
            //          totalprob = Σ w_q * L        (total probability)
            //
            // type_weighted_hess = L * (d²(log L)/dθ² + d(log L)/dθ * d(log L)/dθ')
            //                    = d²L/dθ²  (this IS the raw Hessian, not log!)
            // type_weighted_grad = L * d(log L)/dθ = dL/dθ  (raw gradient)
            // type_weighted_prob = L
            //
            // So: totalhess += quad_weight * type_weighted_hess = w_q * d²L/dθ²
            //     totalgrad += quad_weight * type_weighted_grad = w_q * dL/dθ
            //     totalprob += quad_weight * type_weighted_prob = w_q * L
            //
            // At the end per observation, convert to log-likelihood:
            //     d²(log L)/dθ² = totalhess/totalprob - totalgrad*totalgrad'/(totalprob²)

            // probilk is the quadrature weight at this point
            double quad_weight = probilk;

            // OPTIMIZATION: For ntyp=1, use grad_this_type/hess_this_type directly
            // instead of going through type_weighted_* intermediates. This combines
            // two O(nparam²) loops into one for the Hessian.
            //
            // For mixtures: accumulate into mixprob/mixgrad/mixhess first
            // then weight by mix_weight when accumulating to total
            if (ntyp == 1) {
                // For ntyp=1, prob_this_type is the likelihood at this quad point
                // (type_prob = 1.0 for single type)
                mixprob += quad_weight * prob_this_type;

                if (iflag >= 2) {
                    // grad_this_type = d(log L)/dθ, so raw gradient = L * grad = prob_this_type * grad_this_type
                    // Accumulate for all params that need gradients (free + tied)
                    double combined_weight = quad_weight * prob_this_type;
                    for (size_t gi = 0; gi < gradparlist.size(); gi++) {
                        int ipar = gradparlist[gi];
                        mixgrad[ipar] += combined_weight * grad_this_type[ipar];
                    }

                    // For mixture mean gradients: accumulate dL_m/d(f_k)
                    // d(log L)/df_k = grad_this_type[var_idx] / df_dsigma2[k] for independent factors
                    // where var_idx = GetFactorVarianceIndex(imix, k) is the variance param index
                    // dL/df_k = L * d(log L)/df_k
                    // Only needed when nmix > 1
                    if (nmix > 1 && factor_structure == FactorStructure::INDEPENDENT) {
                        for (int k = 0; k < n_factors_for_mixture; k++) {
                            // df_dsigma2[k] = x_node / (2*sigma), non-zero when x_node != 0
                            // Avoid division by zero: if x_node = 0, then df/dsigma2 = 0
                            // and grad_this_type[k] contribution from df/dsigma2 is also 0
                            if (std::fabs(df_dsigma2[k]) > 1e-12) {
                                // Use mixture-specific variance index to get the correct gradient
                                int var_idx = GetFactorVarianceIndex(imix, k);
                                double dlogL_dfk = grad_this_type[var_idx] / df_dsigma2[k];
                                mix_dL_df[k] += combined_weight * dlogL_dfk;
                            }
                        }
                    }
                }

                if (iflag == 3) {
                    // Raw Hessian = L * (d²(log L)/dθ² + d(log L)/dθ * d(log L)/dθ')
                    //             = prob_this_type * (hess_this_type + grad_this_type * grad_this_type')
                    double combined_weight = quad_weight * prob_this_type;
                    // Include tied params under equality constraints (see comment
                    // at the outer totalhess zeroing loop).
                    for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                        int i = gradparlist[fi];
                        for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                            int j = gradparlist[fj];
                            int idx = i * nparam + j;
                            mixhess[idx] += combined_weight *
                                (hess_this_type[idx] + grad_this_type[i] * grad_this_type[j]);
                        }
                    }

                    // For mixture Hessians: convert factor-factor log-likelihood Hessian to raw Hessian
                    // and accumulate into mix_d2L_dfdf and mix_d2L_dsigma_df.
                    // d2logL_dfdf was accumulated directly from modHess during the model loop.
                    // Raw Hessian: d²L/d(f_k1)d(f_k2) = L * (d²(log L)/df df + dlogL/df * dlogL/df')
                    if (nmix > 1 && factor_structure == FactorStructure::INDEPENDENT) {
                        for (int k1 = 0; k1 < n_factors_for_mixture; k1++) {
                            for (int k2 = k1; k2 < n_factors_for_mixture; k2++) {
                                int ff_idx = k1 * n_factors_for_mixture + k2;

                                // Convert log-likelihood Hessian to raw Hessian
                                // d²L/df df = L * (d²logL/df df + dlogL/df_k1 * dlogL/df_k2)
                                double d2L_df_df = prob_this_type *
                                    (d2logL_dfdf[ff_idx] + dlogL_df[k1] * dlogL_df[k2]);

                                // Accumulate into mix_d2L_dfdf (weighted by quadrature weight)
                                mix_d2L_dfdf[imix][ff_idx] += quad_weight * d2L_df_df;

                                // Accumulate into mix_d2L_dsigma_df for variance-mean cross Hessian
                                // d²L/(d(σ²_k1) d(f_k2)) = d²L/(df_k1 df_k2) * df_k1/d(σ²_k1)
                                //                       = d2L_df_df * df_dsigma2[k1]
                                mix_d2L_dsigma_df[imix][ff_idx] += quad_weight * d2L_df_df * df_dsigma2[k1];

                                // Symmetric entry
                                if (k1 != k2) {
                                    int ff_idx_sym = k2 * n_factors_for_mixture + k1;
                                    mix_d2L_dfdf[imix][ff_idx_sym] += quad_weight * d2L_df_df;
                                    // For variance-mean cross: d²L/(d(σ²_k2) d(f_k1)) = d2L_df_df * df_dsigma2[k2]
                                    mix_d2L_dsigma_df[imix][ff_idx_sym] += quad_weight * d2L_df_df * df_dsigma2[k2];
                                }
                            }
                        }

                        // Convert factor × model_param LOG-Hessian to RAW Hessian and accumulate
                        // d²L/df_k d(param_j) = L * (d²logL/df_k d(param_j) + dlogL/df_k * dlogL/d(param_j))
                        // Where dlogL/df_k and dlogL/d(param_j) are TOTALS from all models
                        for (int k = 0; k < n_factors_for_mixture; k++) {
                            for (int j = 0; j < n_model_params; j++) {
                                int storage_idx = k * n_model_params + j;
                                int param_idx = first_model_param + j;
                                // Convert log-Hessian to raw Hessian
                                double d2L_df_dbeta = prob_this_type *
                                    (d2logL_df_dbeta_thisintpt[storage_idx] +
                                     dlogL_df[k] * grad_this_type[param_idx]);
                                // Accumulate weighted by quadrature weight
                                mix_d2L_df_dbeta[imix][storage_idx] += quad_weight * d2L_df_dbeta;
                            }
                        }

                        // Reset accumulators for next integration point
                        std::fill(d2logL_dfdf.begin(), d2logL_dfdf.end(), 0.0);
                        std::fill(dlogL_df.begin(), dlogL_df.end(), 0.0);
                        std::fill(d2logL_df_dbeta_thisintpt.begin(), d2logL_df_dbeta_thisintpt.end(), 0.0);
                    }
                }
            } else {
                // ntyp > 1: Use type_weighted_* which already accumulated across types
                mixprob += quad_weight * type_weighted_prob;

                if (iflag >= 2) {
                    // Accumulate for all params that need gradients (free + tied)
                    for (size_t gi = 0; gi < gradparlist.size(); gi++) {
                        int ipar = gradparlist[gi];
                        mixgrad[ipar] += quad_weight * type_weighted_grad[ipar];
                    }
                }

                if (iflag == 3) {
                    // Include tied params under equality constraints (see comment
                    // at the outer totalhess zeroing loop).
                    for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                        int i = gradparlist[fi];
                        for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                            int j = gradparlist[fj];
                            int idx = i * nparam + j;
                            mixhess[idx] += quad_weight * type_weighted_hess[idx];
                        }
                    }
                }
            }

            // ===== STEP 11: Increment integration point indices (odometer style) =====
            if (intpt < nint_points - 1) {
                for (int ifac = 0; ifac < nfac; ifac++) {
                    facint[ifac]++;
                    if (facint[ifac] < obs_nq[ifac]) break;
                    facint[ifac] = 0;
                }
            }
        }  // End integration point loop

        // ===== STEP 11b: Accumulate mixture contribution =====
        // Weight this mixture's contribution by mix_weight
        // Note: For nmix=1, mix_weight=1 so this is equivalent to old code
        totalprob += mix_weight * mixprob;

        if (iflag >= 2) {
            // Mixture gradient contribution
            // d/dθ[Σ_m w_m L_m] = Σ_m [dw_m/dθ * L_m + w_m * dL_m/dθ]
            // For now, just accumulate w_m * dL_m/dθ
            // Gradient w.r.t. mixture weights/means handled separately below
            for (size_t gi = 0; gi < gradparlist.size(); gi++) {
                int ipar = gradparlist[gi];
                totalgrad[ipar] += mix_weight * mixgrad[ipar];
            }
        }

        if (iflag == 3) {
            // Mixture Hessian contribution. Include tied params under equality
            // constraints (see comment at the outer totalhess zeroing loop).
            for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                int i = gradparlist[fi];
                for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                    int j = gradparlist[fj];
                    int idx = i * nparam + j;
                    totalhess[idx] += mix_weight * mixhess[idx];
                }
            }
        }

        // ===== STEP 11c: Store per-mixture values for mixture parameter gradients =====
        // Store L_m for this mixture (needed for weight parameter gradients)
        mix_L[imix] = mixprob;

        // Store dL_m/d(mean_k) for this mixture (needed for mean parameter gradients)
        // dL_m/d(mean_k) = Σ_intpt w_q * L_intpt * d(log L_intpt)/d(f_k) * d(f_k)/d(mean_k)
        //                = Σ_intpt w_q * L_intpt * d(log L_intpt)/d(f_k) * 1
        // Store mix_dL_df for this mixture (for mean gradient computation after loop)
        if (nmix > 1 && iflag >= 2) {
            for (int k = 0; k < n_factors_for_mixture; k++) {
                mix_dL_dmean[imix][k] = mix_dL_df[k];
                // Store dL_m/d(σ²_k) for variance-weight cross Hessian
                // This is the unweighted variance gradient accumulated in mixgrad
                int var_idx = GetFactorVarianceIndex(imix, k);
                mix_dL_dsigma[imix][k] = mixgrad[var_idx];
            }
        }
        // Store dL_m/d(beta_j) for this mixture (for weight-beta cross Hessian)
        if (nmix > 1 && iflag == 3) {
            for (int j = 0; j < n_model_params; j++) {
                int beta_idx = first_model_param + j;
                mix_dL_dbeta[imix][j] = mixgrad[beta_idx];
            }
        }

        }  // End mixture loop

        // ===== STEP 11d: Compute mixture parameter gradients after all mixtures =====
        // L = Σ_m w_m * L_m
        //
        // Gradient w.r.t. log-weight parameters (fw_m for m < nmix-1):
        // dL/d(fw_m) = Σ_n (dw_n/d(fw_m)) * L_n
        // Using softmax: dw_m/d(fw_m) = w_m(1 - w_m), dw_n/d(fw_m) = -w_m * w_n for n ≠ m
        // So dL/d(fw_m) = w_m(1 - w_m) * L_m + Σ_{n≠m} (-w_m * w_n) * L_n
        //              = w_m * L_m - w_m * (w_m * L_m + Σ_{n≠m} w_n * L_n)
        //              = w_m * (L_m - L)  where L = totalprob
        //
        // Gradient w.r.t. mean parameters (mu_{m,k} for m < nmix-1, k < n_factors_for_mixture):
        // L_m depends on mean through the factor distribution shift
        // dL/d(mu_{m,k}) = w_m * dL_m/d(mu_{m,k}) + w_last * dL_last/d(mu_{last,k}) * d(mu_{last,k})/d(mu_{m,k})
        // where d(mu_{last,k})/d(mu_{m,k}) = -w_m / w_last (from constraint Σ w_m * mu_m = 0)
        // = w_m * dL_m/d(mu_{m,k}) - w_m * dL_last/d(mu_{last,k})
        // where dL_m/d(mu_{m,k}) is the derivative of L_m w.r.t. the mixture-specific mean
        //
        // Note: We don't have dL_m/d(mu_{m,k}) stored directly. We would need to modify
        // the integration loop to accumulate this. For now, we'll use the fact that
        // the mean gradient equals the sum over integration points of w_q * dL_intpt/df_k
        // This requires re-computing, which is expensive. A proper implementation would
        // accumulate during the main loop.
        //
        // SIMPLIFIED APPROACH: For the initial implementation, we note that:
        // If we shift the mean for mixture m, factor values shift by the same amount.
        // dL_m/d(mu_k) = Σ_intpt w_q * L_intpt * d(log L_intpt)/d(f_k)
        // This is proportional to mixgrad[k] / df_dsigma2, but df_dsigma2 varies.
        // As a first approximation, we can use the average.
        if (nmix > 1 && iflag >= 2) {
            // Log-weight gradients
            // dL/d(fw_m) = w_m * (L_m - L) + dL/d(mean_last) * d(mean_last)/d(fw_m)
            //
            // The second term accounts for the fact that changing mixture weights
            // also changes the derived mean for the last mixture via constraint:
            // mean_last = -Σ_j w_j * mean_j / w_last
            //
            // d(mean_last_k)/d(fw_m) = -mean_m_k * w_m / w_last  (for independent factors)
            //
            // So the additional contribution is:
            // Σ_k dL/d(mean_last_k) * (-mean_m_k * w_m / w_last)
            // = -w_m / w_last * Σ_k dL_last/d(f_k) * mean_m_k
            // = -w_m / w_last * Σ_k mix_dL_dmean[last][k] * mix_means[m][k]
            //
            int last_mix = nmix - 1;
            for (int m = 0; m < nmix - 1; m++) {
                int logwt_idx = GetMixtureLogWeightIndex(m);
                // Direct effect: changing weights
                double dL_dfw = mix_weights[m] * (mix_L[m] - totalprob);

                // Indirect effect: changing mean_last through constraint
                // Only applies for independent factors
                //
                // For 2-mixture softmax with constraint mean_last = -w_m*mean_m/w_last:
                // d(mean_last)/d(fw_m) = -mean_m * d(w_m/w_last)/d(fw_m)
                //                      = -mean_m * w_m/w_last  (since d(w_m/w_last)/d(fw_m) = w_m/w_last for 2-mix)
                //
                // The indirect effect on L is:
                // w_last * dL_last/d(mean_last) * d(mean_last)/d(fw_m)
                // = w_last * mix_dL_dmean[last][k] * (-mean_m * w_m/w_last)
                // = -w_m * mean_m * mix_dL_dmean[last][k]
                //
                if (factor_structure == FactorStructure::INDEPENDENT) {
                    for (int k = 0; k < n_factors_for_mixture; k++) {
                        // mix_dL_dmean[last_mix][k] is dL_last/d(mean_last_k)
                        // mix_means[m][k] is the mean for mixture m, factor k
                        double indirect = -mix_weights[m] * mix_means[m][k] * mix_dL_dmean[last_mix][k];
                        dL_dfw += indirect;
                    }
                }

                totalgrad[logwt_idx] += dL_dfw;
            }

            // Mean gradients (for independent factors only for now)
            // dL/d(mu_{m,k}) = w_m * dL_m/d(mu_k) + w_last * dL_last/d(mu_k) * d(mu_last_k)/d(mu_{m,k})
            // where d(mu_last_k)/d(mu_{m,k}) = -w_m / w_last (from constraint)
            // = w_m * (dL_m/d(mu_k) - dL_last/d(mu_k))
            if (factor_structure == FactorStructure::INDEPENDENT) {
                int last_mix = nmix - 1;
                for (int m = 0; m < nmix - 1; m++) {
                    for (int k = 0; k < n_factors_for_mixture; k++) {
                        int mean_idx = GetFactorMeanIndex(m, k);
                        // dL_m/d(mu_k) and dL_last/d(mu_k) are stored in mix_dL_dmean
                        double dL_dmu = mix_weights[m] * (mix_dL_dmean[m][k] - mix_dL_dmean[last_mix][k]);
                        totalgrad[mean_idx] += dL_dmu;
                    }
                }
            }

            // ===== STEP 11e: Compute mixture parameter Hessians after all mixtures =====
            // These are the second derivatives of L = Σ_m w_m * L_m w.r.t. mixture parameters
            // Using mix_d2L_dfdf (factor-factor Hessian per mixture) accumulated during integration
            // Only access mix_d2L_df_dbeta when it's been initialized (iflag == 3)
            // (removing debug output for now)
            if (iflag == 3 && factor_structure == FactorStructure::INDEPENDENT && ntyp == 1) {
                int last_mix = nmix - 1;

                // (1) Weight-weight Hessian: d²L/d(fw_m1)d(fw_m2)
                // From L = Σ_n w_n * L_n, the raw Hessian is:
                // d²L/d(fw_m)² = Σ_n [d²w_n/d(fw_m)² * L_n + 2 * dw_n/d(fw_m) * dL_n/d(fw_m) + w_n * d²L_n/d(fw_m)²]
                //
                // Term 1: Softmax second derivative contribution
                // For 2-mixture: d²w_0/d(fw_0)² = w_0*(1-w_0)*(1-2*w_0)
                //                d²w_1/d(fw_0)² = -w_0*(1-w_0)*(1-2*w_0)
                // So: Σ_n d²w_n/d(fw_0)² * L_n = w_0*(1-w_0)*(1-2*w_0) * (L_0 - L_1)
                //
                // Term 2: Cross effect (only last mixture contributes)
                // 2 * dw_last/d(fw_m) * dL_last/d(fw_m)
                // = 2 * (-w_last * w_m) * (-exp(fw_m) * Σ_k mix_dL_dmean[last][k] * mean_m_k)
                // = 2 * w_last * w_m * exp(fw_m) * Σ_k mix_dL_dmean[last][k] * mean_m_k
                //
                // Term 3: w_last * d²L_last/d(fw_m)² (through mean_last dependence)
                for (int m1 = 0; m1 < nmix - 1; m1++) {
                    int fw_idx1 = GetMixtureLogWeightIndex(m1);
                    double exp_fw_m1 = std::exp(param[fw_idx1]);
                    double w_m1 = mix_weights[m1];

                    for (int m2 = m1; m2 < nmix - 1; m2++) {
                        int fw_idx2 = GetMixtureLogWeightIndex(m2);
                        double exp_fw_m2 = std::exp(param[fw_idx2]);
                        double w_m2 = mix_weights[m2];
                        int hess_idx = fw_idx1 * nparam + fw_idx2;

                        double contrib = 0.0;

                        // Term 1: Softmax second derivative contribution
                        // For diagonal (m1 == m2):
                        // d²w_m/d(fw_m)² = w_m * (1-w_m) * (1 - 2*w_m)
                        // d²w_n/d(fw_m)² = w_n * w_m * (2*w_m - 1) for n != m
                        // Σ_n d²w_n/d(fw_m)² * L_n computed below
                        if (m1 == m2) {
                            // Diagonal: compute Σ_n d²w_n/d(fw_m)² * L_n
                            double d2w_m_dfw2 = w_m1 * (1.0 - w_m1) * (1.0 - 2.0 * w_m1);
                            double d2w_last_dfw2 = mix_weights[last_mix] * w_m1 * (2.0 * w_m1 - 1.0);

                            // Contribution from non-last mixtures
                            contrib += d2w_m_dfw2 * mix_L[m1];

                            // Contribution from other non-last mixtures (n != m, n < last)
                            for (int n = 0; n < nmix - 1; n++) {
                                if (n != m1) {
                                    // d²w_n/d(fw_m)² = w_n * w_m * (2*w_m - 1) for n != m
                                    double d2w_n_dfw2 = mix_weights[n] * w_m1 * (2.0 * w_m1 - 1.0);
                                    contrib += d2w_n_dfw2 * mix_L[n];
                                }
                            }

                            // Contribution from last mixture
                            contrib += d2w_last_dfw2 * mix_L[last_mix];
                        } else {
                            // Off-diagonal: d²w_m/d(fw_m1)d(fw_m2) terms
                            // d²w_m1/d(fw_m1)d(fw_m2) = w_m1 * (δ_{m1} - w_m1) * (-w_m2) = -w_m1 * (1-w_m1) * w_m2
                            // d²w_m2/d(fw_m1)d(fw_m2) = w_m2 * (-w_m1) * (δ_{m2} - w_m2) = -w_m2 * w_m1 * (1-w_m2)
                            // d²w_last/d(fw_m1)d(fw_m2) = -(-w_last*w_m1) * w_m2 - w_last * (-w_m1*w_m2)
                            //                          = w_last*w_m1*w_m2 + w_last*w_m1*w_m2 = 2*w_last*w_m1*w_m2
                            double d2w_m1_dfw12 = -w_m1 * (1.0 - w_m1) * w_m2;
                            double d2w_m2_dfw12 = -w_m2 * w_m1 * (1.0 - w_m2);
                            double d2w_last_dfw12 = 2.0 * mix_weights[last_mix] * w_m1 * w_m2;

                            contrib += d2w_m1_dfw12 * mix_L[m1] + d2w_m2_dfw12 * mix_L[m2];
                            contrib += d2w_last_dfw12 * mix_L[last_mix];
                        }

                        // Term 2: Cross effect (2 * dw_last/d(fw_m) * dL_last/d(fw_m))
                        // dw_last/d(fw_m) = -w_last * w_m
                        // dL_last/d(fw_m) = Σ_k mix_dL_dmean[last][k] * (-exp(fw_m)*mean_m_k)
                        if (m1 == m2) {
                            double dw_last_dfw = -mix_weights[last_mix] * w_m1;
                            double dL_last_dfw = 0.0;
                            for (int k = 0; k < n_factors_for_mixture; k++) {
                                dL_last_dfw += mix_dL_dmean[last_mix][k] * (-exp_fw_m1 * mix_means[m1][k]);
                            }
                            contrib += 2.0 * dw_last_dfw * dL_last_dfw;
                        }

                        // Term 3: w_last * d²L_last/d(fw_m1)d(fw_m2) through mean_last
                        // d²L_last/d(fw_m1)d(fw_m2) = Σ_k1,k2 d²L_last/d(f_k1)d(f_k2) * df_k1/d(fw_m1) * df_k2/d(fw_m2)
                        //                          + δ_{m1,m2} * Σ_k dL_last/d(f_k) * d²(mean_last_k)/d(fw_m)²
                        // where df_k/d(fw_m) = d(mean_last_k)/d(fw_m) = -exp(fw_m)*mean_m_k
                        for (int k1 = 0; k1 < n_factors_for_mixture; k1++) {
                            double df_k1_dfw_m1 = -exp_fw_m1 * mix_means[m1][k1];
                            for (int k2 = 0; k2 < n_factors_for_mixture; k2++) {
                                double df_k2_dfw_m2 = -exp_fw_m2 * mix_means[m2][k2];
                                int ff_idx = k1 * n_factors_for_mixture + k2;
                                contrib += mix_weights[last_mix] * mix_d2L_dfdf[last_mix][ff_idx] * df_k1_dfw_m1 * df_k2_dfw_m2;
                            }
                        }

                        // Second derivative of mean_last w.r.t. fw_m (only for diagonal)
                        // d²(mean_last_k)/d(fw_m)² = d(-exp(fw_m)*mean_m_k)/d(fw_m) = -exp(fw_m)*mean_m_k
                        if (m1 == m2) {
                            for (int k = 0; k < n_factors_for_mixture; k++) {
                                double d2_mean_last_dfw2 = -exp_fw_m1 * mix_means[m1][k];
                                contrib += mix_weights[last_mix] * mix_dL_dmean[last_mix][k] * d2_mean_last_dfw2;
                            }
                        }

                        totalhess[hess_idx] += contrib;
                    }
                }

                // (2) Mean-mean Hessian: d²L/d(mu_m1_k1)d(mu_m2_k2)
                // dL/d(mu_m_k) = w_m * (dL_m/d(f_k) - dL_last/d(f_k))
                // d²L/d(mu_m1_k1)d(mu_m2_k2) =
                //   δ_{m1,m2} * w_m1 * d²L_m1/d(f_k1)d(f_k2)   [direct: same mixture]
                // + (w_m1 * w_m2 / w_last) * d²L_last/d(f_k1)d(f_k2)   [indirect: via mean_last]
                // - δ_{m1,m2} * w_m1 * d²L_last/d(f_k1)d(f_k2) * (-w_m1/w_last)   [cross term]
                //
                // Simplifying:
                // = δ_{m1,m2} * w_m1 * d²L_m1/d(f_k1)d(f_k2)
                // + (w_m1 * w_m2 / w_last) * d²L_last/d(f_k1)d(f_k2)
                for (int m1 = 0; m1 < nmix - 1; m1++) {
                    for (int k1 = 0; k1 < n_factors_for_mixture; k1++) {
                        int mean_idx1 = GetFactorMeanIndex(m1, k1);

                        for (int m2 = m1; m2 < nmix - 1; m2++) {
                            for (int k2 = (m1 == m2 ? k1 : 0); k2 < n_factors_for_mixture; k2++) {
                                int mean_idx2 = GetFactorMeanIndex(m2, k2);
                                int hess_idx = mean_idx1 * nparam + mean_idx2;
                                if (mean_idx1 > mean_idx2) {
                                    hess_idx = mean_idx2 * nparam + mean_idx1;
                                }

                                int ff_idx = k1 * n_factors_for_mixture + k2;
                                double contrib = 0.0;

                                // Direct term: same mixture
                                if (m1 == m2) {
                                    contrib += mix_weights[m1] * mix_d2L_dfdf[m1][ff_idx];
                                }

                                // Indirect term: via last mixture's mean constraint
                                // d(mean_last_k1)/d(mu_m1_k1) = -w_m1/w_last
                                // d(mean_last_k2)/d(mu_m2_k2) = -w_m2/w_last
                                contrib += (mix_weights[m1] * mix_weights[m2] / mix_weights[last_mix])
                                         * mix_d2L_dfdf[last_mix][ff_idx];

                                totalhess[hess_idx] += contrib;
                            }
                        }
                    }
                }

                // (3) Weight-mean cross Hessian: d²L/d(fw_m)d(mu_n_k)
                // From dL/d(fw_m) = w_m*(L_m - L) + Σ_j dL_last/df_j * (-w_m * mean_m_j)
                // and dL/d(mu_n_k) = w_n*(dL_n/d(f_k) - dL_last/d(f_k))
                //
                // Taking derivative of dL/d(fw_m) w.r.t. mu_n_k:
                // Term A (direct): w_m * (dL_m/d(mu_n_k) - dL/d(mu_n_k))
                // Term B (indirect second-order): Σ_j d²L_last/df_j d(f_k) * (-w_m*mean_m_j) * d(f_k)/d(mu_n_k)
                // Term C (indirect first-order, m=n only): dL_last/df_k * d(-w_m*mean_m_k)/d(mu_n_k)
                for (int m = 0; m < nmix - 1; m++) {
                    int fw_idx = GetMixtureLogWeightIndex(m);

                    for (int n = 0; n < nmix - 1; n++) {
                        for (int k = 0; k < n_factors_for_mixture; k++) {
                            int mean_idx = GetFactorMeanIndex(n, k);
                            int hess_idx = std::min(fw_idx, mean_idx) * nparam + std::max(fw_idx, mean_idx);

                            // Term A: w_m * (dL_m/d(mu_n_k) - dL/d(mu_n_k))
                            // dL_m/d(mu_n_k) = dL_m/df_k if m=n, else 0
                            // dL/d(mu_n_k) = w_n * (dL_n/df_k - dL_last/df_k)
                            double dL_m_dmu_n_k = (m == n) ? mix_dL_dmean[m][k] : 0.0;
                            double dL_dmu_n_k = mix_weights[n] * (mix_dL_dmean[n][k] - mix_dL_dmean[last_mix][k]);
                            double direct_contrib = mix_weights[m] * (dL_m_dmu_n_k - dL_dmu_n_k);

                            // Term B: second-order contribution from last mixture
                            // d(f_k)/d(mu_n_k) affects last mixture: = -w_n/w_last
                            // d(f_j)/d(fw_m) affects last mixture: = -w_m*mean_m_j / w_last
                            // Combined with w_last: w_last * (-w_m*mean_m_j/w_last) * (-w_n/w_last)
                            //                     = w_m * w_n * mean_m_j / w_last
                            double second_order_contrib = 0.0;
                            for (int j = 0; j < n_factors_for_mixture; j++) {
                                int ff_idx = j * n_factors_for_mixture + k;
                                second_order_contrib += mix_d2L_dfdf[last_mix][ff_idx] *
                                    (mix_weights[m] * mix_weights[n] / mix_weights[last_mix]) * mix_means[m][j];
                            }

                            // Term C: first derivative term when m == n
                            // d(-w_m*mean_m_k)/d(mu_m_k) = -w_m (since mean_m_k is parameter)
                            // w_last cancels with d(mean_last)/d(mu_n) = -w_n/w_last, giving -w_n
                            // But we're differentiating w.r.t. mu_n_k, so this only contributes when m=n
                            double first_order_contrib = 0.0;
                            if (m == n) {
                                first_order_contrib = -mix_weights[m] * mix_dL_dmean[last_mix][k];
                            }

                            totalhess[hess_idx] += direct_contrib + second_order_contrib + first_order_contrib;
                        }
                    }
                }

                // (4) Variance-mean cross Hessian: d²L/d(var_imix_k1)d(mu_m_k2)
                // mix_d2L_dsigma_df[imix][k1,k2] = d²L_imix/(d(σ²_k1) df_k2) was accumulated during integration
                // For variance var_imix_k1: only mixture imix contributes (variance only affects its own mixture)
                // d²L/d(var_imix_k1)d(mu_m_k2) = w_imix * mix_d2L_dsigma_df[imix][k1,k2] * df_k2/d(mu_m_k2)
                //
                // df_k2/d(mu_m_k2) depends on which mixture:
                //   - For imix == m (non-last): df_k2/d(mu_m_k2) = 1
                //   - For imix == last: df_k2/d(mu_m_k2) = d(mu_last_k2)/d(mu_m_k2) = -w_m/w_last
                for (int imix = 0; imix < nmix; imix++) {
                    for (int k1 = 0; k1 < n_factors_for_mixture; k1++) {
                        int var_idx = GetFactorVarianceIndex(imix, k1);

                        for (int m = 0; m < nmix - 1; m++) {
                            for (int k2 = 0; k2 < n_factors_for_mixture; k2++) {
                                int mean_idx = GetFactorMeanIndex(m, k2);
                                int hess_idx = std::min(var_idx, mean_idx) * nparam + std::max(var_idx, mean_idx);
                                int ff_idx = k1 * n_factors_for_mixture + k2;

                                double contrib = 0.0;
                                if (imix == m) {
                                    // Same mixture: df_k2/d(mu_m_k2) = 1
                                    contrib = mix_weights[imix] * mix_d2L_dsigma_df[imix][ff_idx];
                                } else if (imix == last_mix) {
                                    // Last mixture: df_k2/d(mu_m_k2) = -w_m/w_last
                                    double df_dmu = -mix_weights[m] / mix_weights[last_mix];
                                    contrib = mix_weights[imix] * mix_d2L_dsigma_df[imix][ff_idx] * df_dmu;
                                }
                                // Other mixtures: variance doesn't affect this mixture, no contribution

                                totalhess[hess_idx] += contrib;
                            }
                        }
                    }
                }

                // (5) Variance-weight cross Hessian: d²L/d(var_imix_k)d(fw_m)
                // d²L/d(σ²_imix_k)d(fw_m) = Σ_n [d²w_n/d(fw_m) * dL_n/d(σ²_imix_k) + dw_n/d(fw_m) * d²L_n/(d(σ²_imix_k)d(fw_m))]
                //
                // But dL_n/d(σ²_imix_k) = 0 for n != imix (variance only affects its own mixture)
                // And d²L_n/(d(σ²_imix_k)d(fw_m)) = 0 for n != last (only last mixture has fw_m dependency through mean_last)
                //
                // So the terms are:
                // 1. dw_imix/d(fw_m) * dL_imix/d(σ²_imix_k)  [first-order cross term]
                // 2. For imix == last: w_last * d²L_last/(d(σ²_k)d(fw_m))  [second-order term through mean_last]
                for (int imix = 0; imix < nmix; imix++) {
                    for (int k = 0; k < n_factors_for_mixture; k++) {
                        int var_idx = GetFactorVarianceIndex(imix, k);

                        for (int m = 0; m < nmix - 1; m++) {
                            int fw_idx = GetMixtureLogWeightIndex(m);
                            int hess_idx = std::min(var_idx, fw_idx) * nparam + std::max(var_idx, fw_idx);

                            double contrib = 0.0;

                            // Term 1: First-order cross term: dw_imix/d(fw_m) * dL_imix/d(σ²_k)
                            // dw_imix/d(fw_m) via softmax:
                            //   If imix == m: w_m * (1 - w_m)
                            //   If imix != m and imix < nmix-1: -w_imix * w_m
                            //   If imix == last: -w_last * w_m
                            double dw_imix_dfw_m;
                            if (imix == m) {
                                dw_imix_dfw_m = mix_weights[m] * (1.0 - mix_weights[m]);
                            } else {
                                dw_imix_dfw_m = -mix_weights[imix] * mix_weights[m];
                            }
                            contrib += dw_imix_dfw_m * mix_dL_dsigma[imix][k];

                            // Term 2: Second-order term for last mixture through mean_last
                            // d²L_last/(d(σ²_k)d(fw_m)) = Σ_j d²L_last/(d(σ²_k)df_j) * df_j/d(fw_m)
                            //   where df_j/d(fw_m) = d(mean_last_j)/d(fw_m) = -exp(fw_m) * mean_m_j
                            if (imix == last_mix) {
                                double exp_fw_m = std::exp(param[fw_idx]);
                                double second_order_contrib = 0.0;
                                for (int j = 0; j < n_factors_for_mixture; j++) {
                                    double df_j_dfw_m = -exp_fw_m * mix_means[m][j];
                                    int sf_idx = k * n_factors_for_mixture + j;
                                    second_order_contrib += mix_d2L_dsigma_df[last_mix][sf_idx] * df_j_dfw_m;
                                }
                                contrib += mix_weights[last_mix] * second_order_contrib;
                            }

                            totalhess[hess_idx] += contrib;
                        }
                    }
                }

                // (6) Mean × Model param cross-Hessian: d²L/d(mean_m_k)d(beta_j)
                // Using mix_d2L_df_dbeta[imix][k,j] = d²L_imix/d(f_k)d(beta_j) accumulated during integration
                //
                // For mean_m_k (m < last):
                // dL/d(mean_m_k) = w_m * dL_m/df_k + w_last * dL_last/df_k * (-w_m/w_last)
                //                = w_m * (dL_m/df_k - dL_last/df_k)
                //
                // d²L/d(mean_m_k)d(beta_j) = w_m * (d²L_m/df_k d(beta_j) - d²L_last/df_k d(beta_j))
                //                          = w_m * (mix_d2L_df_dbeta[m][k,j] - mix_d2L_df_dbeta[last][k,j])

                for (int m = 0; m < nmix - 1; m++) {
                    for (int k = 0; k < n_factors_for_mixture; k++) {
                        int mean_idx = GetFactorMeanIndex(m, k);

                        for (int j = 0; j < n_model_params; j++) {
                            int beta_idx = first_model_param + j;
                            int hess_idx = std::min(mean_idx, beta_idx) * nparam + std::max(mean_idx, beta_idx);
                            int storage_idx = k * n_model_params + j;

                            double contrib = mix_weights[m] *
                                (mix_d2L_df_dbeta[m][storage_idx] - mix_d2L_df_dbeta[last_mix][storage_idx]);
                            totalhess[hess_idx] += contrib;
                        }
                    }
                }

                // (7) Weight × Model param cross-Hessian: d²L/d(fw_m)d(beta_j)
                // From gradient: dL/d(fw_m) = w_m * (L_m - L) + Σ_k dL_last/df_k * (-w_m * mean_m_k)
                //
                // Taking derivative w.r.t. beta_j:
                // d²L/d(fw_m)d(beta_j) = w_m * (dL_m/d(beta_j) - dL/d(beta_j))   [Term A: direct]
                //                      + Σ_k d²L_last/d(f_k)d(beta_j) * (-w_m * mean_m_k)  [Term B: indirect]
                //
                // Term A: direct effect from changing weights
                // Term B: indirect effect from mean_last changing with fw_m
                for (int m = 0; m < nmix - 1; m++) {
                    int fw_idx = GetMixtureLogWeightIndex(m);

                    for (int j = 0; j < n_model_params; j++) {
                        int beta_idx = first_model_param + j;
                        int hess_idx = std::min(fw_idx, beta_idx) * nparam + std::max(fw_idx, beta_idx);

                        // Term A: w_m * (dL_m/d(beta_j) - dL/d(beta_j))
                        // dL/d(beta_j) = Σ_n w_n * dL_n/d(beta_j) = totalgrad[beta_idx] (after weighting)
                        // But totalgrad hasn't been fully accumulated yet at this point in the obs loop
                        // So we compute: dL/d(beta_j) = Σ_n w_n * mix_dL_dbeta[n][j]
                        double dL_dbeta = 0.0;
                        for (int n = 0; n < nmix; n++) {
                            dL_dbeta += mix_weights[n] * mix_dL_dbeta[n][j];
                        }
                        double direct_contrib = mix_weights[m] * (mix_dL_dbeta[m][j] - dL_dbeta);

                        // Term B: Σ_k d²L_last/d(f_k)d(beta_j) * (-w_m * mean_m_k)
                        double indirect_contrib = 0.0;
                        for (int k = 0; k < n_factors_for_mixture; k++) {
                            int storage_idx = k * n_model_params + j;
                            indirect_contrib += mix_d2L_df_dbeta[last_mix][storage_idx] * (-mix_weights[m] * mix_means[m][k]);
                        }

                        totalhess[hess_idx] += direct_contrib + indirect_contrib;
                    }
                }
            }
        }

        // ===== STEP 12: Compute log-likelihood for this observation =====
        // Apply observation weight (default 1.0 if no weights set)
        double obs_weight = use_weights ? obs_weights[iobs] : 1.0;

        if (totalprob > 1e-100) {
            logLkhd += obs_weight * std::log(totalprob);

            // Gradient: d(log L)/dθ = (1/L) * dL/dθ (weighted by observation)
            // Accumulate for all params that need gradients (free + tied)
            if (iflag >= 2) {
                for (size_t gi = 0; gi < gradparlist.size(); gi++) {
                    int ipar = gradparlist[gi];
                    full_gradL[ipar] += obs_weight * totalgrad[ipar] / totalprob;
                }
            }

            // Hessian: d²(log L)/dθ² = (1/L)*d²L/dθ² - (1/L²)*(dL/dθ)² (weighted)
            // Include tied params under equality constraints: with ties the
            // analytical Hessian at a primary position must accumulate the
            // cross-derivatives ∂²L/∂primary∂derived, ∂²L/∂derived∂primary,
            // and ∂²L/∂derived∂derived (see ExtractFreeHessian for the
            // aggregation step). That requires populating the
            // full-parameter-space Hessian at tied positions too, which is
            // why we iterate gradparlist here rather than freeparlist.
            if (iflag == 3) {
                for (size_t fi = 0; fi < gradparlist.size(); fi++) {
                    int i = gradparlist[fi];
                    for (size_t fj = fi; fj < gradparlist.size(); fj++) {
                        int j = gradparlist[fj];
                        int idx = i * nparam + j;
                        double term1 = totalhess[idx] / totalprob;
                        double term2 = (totalgrad[i] * totalgrad[j]) / (totalprob * totalprob);
                        full_hessL[idx] += obs_weight * (term1 - term2);
                    }
                }
            }
        } else {
            // Numerical underflow
            logLkhd += obs_weight * (-1e10);
        }
    }

    // ===== STEP 13: Extract free parameter gradients/Hessians =====
    if (iflag >= 2) {
        ExtractFreeGradient(full_gradL, gradL);
    }
    if (iflag == 3) {
        // Under equality constraints, ExtractFreeHessian iterates all (i, j)
        // pairs and reads full_hess[i * nparam + j]. The per-observation
        // accumulation above only populated the upper triangle. For tied
        // parameters, the aggregation needs to find the symmetric entry on
        // the lower triangle too (e.g. when derived i > primary j, the
        // populated upper-triangle entry lives at (j, i), not (i, j)).
        // Symmetrising full_hessL here makes the aggregation lookup correct
        // in both triangles.
        if (!equality_mapping.empty()) {
            for (int i = 0; i < nparam; i++) {
                for (int j = i + 1; j < nparam; j++) {
                    full_hessL[j * nparam + i] = full_hessL[i * nparam + j];
                }
            }
        }
        ExtractFreeHessian(full_hessL, hessL);
    }
}

int FactorModel::GetFactorVarianceIndex(int imix, int ifac)
{
    // Factor variances come first in parameter vector
    // Layout: [var_m0_f0, var_m0_f1, ..., var_m0_fK-1, var_m1_f0, ..., var_mN-1_fK-1]
    // For SE models, only input factors (n_factors_for_mixture) have variances
    return imix * n_factors_for_mixture + ifac;
}

int FactorModel::GetFactorMeanIndex(int imix, int ifac)
{
    // Factor means come after all variances (if nmix > 1)
    // Layout: [variances...] [mean_m0_f0, mean_m0_f1, ..., mean_mN-2_fK-1]
    // Note: imix < nmix-1 (last mixture mean is derived from constraint)
    if (nmix == 1) return -1;  // No means for single mixture
    if (imix >= nmix - 1) return -1;  // Last mixture mean is derived
    int n_var_params = nmix * n_factors_for_mixture;
    return n_var_params + imix * n_factors_for_mixture + ifac;
}

int FactorModel::GetMixtureLogWeightIndex(int imix)
{
    // Log-weights come after variances and means (if nmix > 1)
    // Layout: [variances...] [means...] [logwt_m0, logwt_m1, ..., logwt_mN-2]
    // Note: imix < nmix-1 (last mixture weight is reference)
    if (nmix == 1) return -1;  // No weights for single mixture
    if (imix >= nmix - 1) return -1;  // Last mixture weight is reference
    int n_var_params = nmix * n_factors_for_mixture;
    int n_mean_params = (nmix - 1) * n_factors_for_mixture;
    return n_var_params + n_mean_params + imix;
}

int FactorModel::GetMixtureWeightIndex(int imix)
{
    // Legacy function - same as GetMixtureLogWeightIndex
    return GetMixtureLogWeightIndex(imix);
}

void FactorModel::ComputeMixtureWeights(std::vector<double>& weights)
{
    // Compute mixture weights from log-weight parameters using softmax
    // w_m = exp(fw_m) / (1 + sum_{j<nmix-1} exp(fw_j)) for m < nmix-1
    // w_{nmix-1} = 1 / (1 + sum_{j<nmix-1} exp(fw_j))
    weights.resize(nmix);

    if (nmix == 1) {
        weights[0] = 1.0;
        return;
    }

    double denom = 1.0;  // Start with 1 for reference mixture
    for (int m = 0; m < nmix - 1; m++) {
        double log_wt = param[GetMixtureLogWeightIndex(m)];
        denom += std::exp(log_wt);
    }

    for (int m = 0; m < nmix - 1; m++) {
        double log_wt = param[GetMixtureLogWeightIndex(m)];
        weights[m] = std::exp(log_wt) / denom;
    }
    weights[nmix - 1] = 1.0 / denom;  // Reference mixture
}

void FactorModel::ComputeMixtureMeans(const std::vector<double>& weights,
                                       std::vector<std::vector<double>>& means)
{
    // Compute mixture means from parameters and constraint E[f] = 0
    // For m < nmix-1: mean[m][k] = param[GetFactorMeanIndex(m, k)]
    // For m = nmix-1: mean[m][k] = -sum_{j<nmix-1}(w_j * mean[j][k]) / w_{nmix-1}
    means.resize(nmix);
    for (int m = 0; m < nmix; m++) {
        means[m].resize(n_factors_for_mixture);
    }

    if (nmix == 1) {
        // Single mixture: means are all 0
        for (int k = 0; k < n_factors_for_mixture; k++) {
            means[0][k] = 0.0;
        }
        return;
    }

    // Read means for non-reference mixtures
    for (int m = 0; m < nmix - 1; m++) {
        for (int k = 0; k < n_factors_for_mixture; k++) {
            means[m][k] = param[GetFactorMeanIndex(m, k)];
        }
    }

    // Derive last mixture's means from constraint: sum(w_m * mu_m) = 0
    // mu_last = -sum_{m<nmix-1}(w_m * mu_m) / w_last
    double w_last = weights[nmix - 1];
    for (int k = 0; k < n_factors_for_mixture; k++) {
        double weighted_sum = 0.0;
        for (int m = 0; m < nmix - 1; m++) {
            weighted_sum += weights[m] * means[m][k];
        }
        means[nmix - 1][k] = -weighted_sum / w_last;
    }
}

void FactorModel::SetModelParameters(const std::vector<double>& params)
{
    if (params.size() != nparam) {
        throw std::runtime_error("Parameter vector size mismatch in SetModelParameters");
    }
    param = params;
}

void FactorModel::CalcLkhdSingleObs(int iobs,
                                    const std::vector<double>& factor_values,
                                    const std::vector<double>& model_params,
                                    double& logLkhd,
                                    std::vector<double>& gradL,
                                    std::vector<double>& hessL,
                                    int iflag,
                                    bool includePrior)
{
    // ===== Factor score estimation mode =====
    // Evaluates likelihood for a single observation at given factor values.
    // No quadrature - factor values are the parameters we're optimizing.
    // Includes factor prior: L = p(y|f,θ) * φ(f|0,σ²)

    if (factor_values.size() != nfac) {
        throw std::runtime_error("Factor values size mismatch");
    }
    if (model_params.size() != nparam) {
        throw std::runtime_error("Model parameters size mismatch");
    }

    // Set the model parameters
    param = model_params;

    // Extract factor variances from param vector
    std::vector<double> factor_var(nfac);
    for (int ifac = 0; ifac < nfac; ifac++) {
        factor_var[ifac] = std::fabs(param[ifac]);  // Enforce positivity
    }

    // Initialize outputs
    logLkhd = 0.0;
    if (iflag >= 2) {
        gradL.resize(nfac, 0.0);
        std::fill(gradL.begin(), gradL.end(), 0.0);
    }
    if (iflag == 3) {
        // Upper triangle storage for nfac x nfac Hessian
        hessL.resize(nfac * (nfac + 1) / 2, 0.0);
        std::fill(hessL.begin(), hessL.end(), 0.0);
    }

    int iobs_offset = iobs * nvar;

    // ===== Part 1: Observation likelihood p(y|f,θ) =====
    // We compute log p(y|f) = sum_m log p_m(y_m|f) directly for numerical stability
    double logObs = 0.0;  // Sum of log-likelihoods
    std::vector<double> gradLogObs(nfac, 0.0);  // d(log L_obs)/df = sum_m d(log L_m)/df
    std::vector<double> hessLogObs;
    if (iflag == 3) {
        hessLogObs.resize(nfac * nfac, 0.0);
    }

    // Temporary storage for model evaluation
    std::vector<double> modEval;
    std::vector<double> modHess;

    // Evaluate all models at the given factor values
    for (size_t imod = 0; imod < models.size(); imod++) {
        int firstpar = param_model_start[imod];

        // For factor score estimation, we only need derivatives w.r.t. factor values
        // Pass iflag to get gradients/Hessians w.r.t. factors
        // No free indices needed here since we only use factor-factor part of Hessian
        models[imod]->Eval(iobs_offset, data, param, firstpar, factor_values,
                          modEval, modHess, iflag, 0.0, nullptr);

        // Get model likelihood
        double Lm = modEval[0];
        if (Lm < 1e-100) {
            Lm = 1e-100;  // Prevent log(0)
        }
        logObs += std::log(Lm);

        // Accumulate gradients: d(log L_obs)/df = sum_m d(log L_m)/df
        // NOTE: Model::Eval returns d(log L_m)/df in modEval[ifac+1], NOT dL_m/df
        // So we simply sum them directly without dividing by Lm
        if (iflag >= 2) {
            for (int ifac = 0; ifac < nfac; ifac++) {
                // modEval[ifac + 1] is d(log L_m)/df_ifac from this model
                gradLogObs[ifac] += modEval[ifac + 1];
            }
        }

        // Accumulate Hessians: d²(log L_obs)/df_i df_j = sum_m d²(log L_m)/df_i df_j
        // NOTE: Model::Eval returns d²(log L_m)/df² in modHess, NOT d²L_m/df²
        // So we simply sum them directly without applying the chain rule
        if (iflag == 3 && modHess.size() > 0) {
            int nDimModHess = nfac + param_model_count[imod];
            // Extract only the factor-factor part of the Hessian
            for (int i = 0; i < nfac; i++) {
                for (int j = i; j < nfac; j++) {
                    int modhess_idx = i * nDimModHess + j;
                    int hess_idx = i * nfac + j;

                    // d²(log L_m)/df_i df_j - directly from model
                    hessLogObs[hess_idx] += modHess[modhess_idx];
                }
            }
        }
    }

    // ===== Part 2: Factor prior φ(f|0,σ²) =====
    // log φ(f|0,σ²) = -0.5 * Σ [f²/σ² + log(2π*σ²)]
    // = -0.5 * Σ [f²/σ² + log(2π) + log(σ²)]
    // Note: Legacy code does NOT include the prior when computing factor score SEs
    // (see TMinLkhd.cc:5332-5336 where nfac_eff=0 when predicting)
    // Set includePrior=false to match legacy behavior.
    double logPrior = 0.0;
    std::vector<double> gradPrior(nfac, 0.0);
    std::vector<double> hessPrior(nfac, 0.0);  // Diagonal only for uncorrelated factors

    if (includePrior) {
        for (int ifac = 0; ifac < nfac; ifac++) {
            double f = factor_values[ifac];
            double sigma2 = factor_var[ifac];

            // Log-prior contribution
            logPrior -= 0.5 * (f * f / sigma2 + std::log(2.0 * M_PI * sigma2));

            // Gradient: d(logPrior)/df = -f/σ²
            if (iflag >= 2) {
                gradPrior[ifac] = -f / sigma2;
            }

            // Hessian: d²(logPrior)/df² = -1/σ²
            if (iflag == 3) {
                hessPrior[ifac] = -1.0 / sigma2;
            }
        }
    }

    // ===== Part 3: Combine log-likelihood and log-prior =====
    // log(L_total) = log(L_obs) + log(L_prior)
    logLkhd = logObs + logPrior;

    // Gradient: d(log L_total)/df = d(log L_obs)/df + d(log L_prior)/df
    if (iflag >= 2) {
        for (int ifac = 0; ifac < nfac; ifac++) {
            gradL[ifac] = gradLogObs[ifac] + gradPrior[ifac];
        }
    }

    // Hessian: d²(log L_total)/df² = d²(log L_obs)/df² + d²(log L_prior)/df²
    if (iflag == 3) {
        int hess_out_idx = 0;
        for (int i = 0; i < nfac; i++) {
            for (int j = i; j < nfac; j++) {
                int hess_idx = i * nfac + j;
                double obs_hess = hessLogObs[hess_idx];

                // Add prior Hessian (diagonal only for uncorrelated factors)
                double prior_hess = (i == j) ? hessPrior[i] : 0.0;

                hessL[hess_out_idx++] = obs_hess + prior_hess;
            }
        }
    }
}

int FactorModel::GetTypeProbInterceptIndex(int ityp)
{
    // Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // Returns: index in parameter vector for typeprob_intercept_{ityp+2}
    // Intercepts come first, then loadings
    if (ntyp <= 1 || type_param_start < 0) {
        return -1;  // No type parameters
    }
    return type_param_start + ityp;
}

int FactorModel::GetTypeLoadingIndex(int ityp, int ifac)
{
    // Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // ifac: 0-based factor index
    // Returns: index in parameter vector for lambda_{ityp+2, ifac+1}
    // Intercepts come first (n_typeprob_intercepts), then loadings
    if (ntyp <= 1 || type_param_start < 0) {
        return -1;  // No type parameters
    }
    return type_param_start + n_typeprob_intercepts + ityp * nfac + ifac;
}

int FactorModel::GetTypeInterceptIndex(int ityp, int model_idx)
{
    // Type-specific intercepts are stored in each model's parameter block
    // They come after the model's base parameters
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // model_idx: 0-based model index
    //
    // Note: The actual implementation depends on how model parameters are organized
    // This returns the index within the model's parameter block
    // The caller needs to add param_model_start[model_idx] to get the full index
    if (ntyp <= 1) {
        return -1;  // No type-specific intercepts
    }
    // Type intercepts come after base model parameters
    // Each model has (ntyp - 1) type-specific intercepts
    // param_model_count[model_idx] includes these, so we need to find the offset
    // Base parameters for model = param_model_count[model_idx] - (ntyp - 1)
    int base_params = param_model_count[model_idx] - (ntyp - 1);
    return base_params + ityp;
}

void FactorModel::ComputeTypeProbabilities(const std::vector<double>& fac, std::vector<double>& probs)
{
    // Compute type probabilities using multinomial logit
    // Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    //
    // P(type=1) = 1 / (1 + sum_{t=2}^T exp(eta_t))
    // P(type=t) = exp(eta_t) / (1 + sum_{t=2}^T exp(eta_t))
    // where eta_t = typeprob_t_intercept + sum_k lambda_t_k * f_k
    //
    // OPTIMIZATION: probs is pre-allocated by caller to avoid allocation per call

    if (ntyp == 1) {
        probs[0] = 1.0;
        return;
    }

    // Compute log-odds for types 2, 3, ... relative to type 1
    // OPTIMIZATION: Use stack-based array for small ntyp, heap otherwise
    double log_odds_stack[16];  // Stack allocation for common case
    double* log_odds = (ntyp <= 17) ? log_odds_stack : new double[ntyp - 1];
    double max_log_odds = 0.0;  // For numerical stability

    for (int t = 0; t < ntyp - 1; t++) {
        // Start with the intercept for this type
        int intercept_idx = GetTypeProbInterceptIndex(t);
        log_odds[t] = param[intercept_idx];

        // Add factor loadings
        for (int k = 0; k < nfac; k++) {
            int param_idx = GetTypeLoadingIndex(t, k);
            log_odds[t] += param[param_idx] * fac[k];
        }
        if (log_odds[t] > max_log_odds) {
            max_log_odds = log_odds[t];
        }
    }

    // Compute denominator with log-sum-exp trick for numerical stability
    double sum_exp = std::exp(-max_log_odds);  // This is exp(0 - max) for type 1
    for (int t = 0; t < ntyp - 1; t++) {
        sum_exp += std::exp(log_odds[t] - max_log_odds);
    }

    // Compute probabilities
    probs[0] = std::exp(-max_log_odds) / sum_exp;  // Type 1 (reference)
    for (int t = 0; t < ntyp - 1; t++) {
        probs[t + 1] = std::exp(log_odds[t] - max_log_odds) / sum_exp;
    }

    // Clean up if we used heap allocation
    if (ntyp > 17) {
        delete[] log_odds;
    }
}
