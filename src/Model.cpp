#include "Model.h"
#include "distributions.h"
#include <cmath>
#include <cstring>
#include <algorithm>

Model::Model(ModelType type, int outcome, int missing,
             const std::vector<int>& regs, int nfac, int ntyp,
             const std::vector<double>& fnorm,
             int nchoice, int nrank, bool params_fixed,
             FactorSpec fspec,
             bool dynamic, int outcome_fac_idx,
             const std::vector<int>& outcome_idxs,
             bool excl_chosen, int rankshare_idx,
             bool uses_types)
    : modtype(type), outcome_idx(outcome), missing_idx(missing),
      regressors(regs), facnorm(fnorm),
      numfac(nfac), numtyp(ntyp),
      numchoice(nchoice), numrank(nrank),
      exclude_chosen(excl_chosen), ranksharevar_idx(rankshare_idx),
      ignore(false), all_params_fixed(params_fixed),
      factor_spec(fspec),
      is_dynamic(dynamic), outcome_factor_idx(outcome_fac_idx),
      use_types(uses_types)
{
    nregressors = regressors.size();

    // Initialize outcome_indices for exploded logit
    if (outcome_idxs.empty()) {
        // Single outcome - create vector with just the outcome_idx
        outcome_indices.push_back(outcome_idx);
    } else {
        outcome_indices = outcome_idxs;
    }

    // Compute number of quadratic and interaction loadings based on factor_spec
    // For dynamic models, exclude the outcome factor from quadratic/interaction terms
    int effective_fac = is_dynamic ? (numfac - 1) : numfac;

    if (factor_spec == FactorSpec::QUADRATIC || factor_spec == FactorSpec::FULL) {
        n_quadratic_loadings = effective_fac;
    } else {
        n_quadratic_loadings = 0;
    }

    if ((factor_spec == FactorSpec::INTERACTIONS || factor_spec == FactorSpec::FULL) && effective_fac >= 2) {
        n_interaction_loadings = effective_fac * (effective_fac - 1) / 2;
    } else {
        n_interaction_loadings = 0;
    }
}

void Model::Eval(int iobs_offset, const std::vector<double>& data,
                 const std::vector<double>& param, int firstpar,
                 const std::vector<double>& fac,
                 std::vector<double>& modEval,
                 std::vector<double>& hess,
                 int flag,
                 double type_intercept,
                 const std::vector<int>* model_free_indices)
{
    // Count free factor loadings FIRST (needed for gradient vector sizing)
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    // Determine size of gradient vector
    // PERFORMANCE: Only resize if vector is too small, reuse existing capacity
    if (flag >= 2) {
        // Layout: 1 (likelihood) + numfac (d/dtheta for factor variances) +
        //         nregressors (d/dbeta) + ifreefac (d/dalpha for free linear loadings) +
        //         n_quadratic_loadings (d/dalpha_quad) + n_interaction_loadings (d/dalpha_inter) +
        //         model-specific parameters
        int ngrad = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

        // Model-specific parameters
        if (modtype == ModelType::LINEAR) ngrad += 1;  // sigma
        if (modtype == ModelType::LOGIT && numchoice > 2) {
            ngrad = 1 + numfac + (numchoice-1)*(nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings);
        }
        if (modtype == ModelType::OPROBIT) ngrad += (numchoice - 1);  // thresholds

        // OPTIMIZATION: Only resize if needed, use memset for faster zeroing
        if (modEval.size() < static_cast<size_t>(ngrad)) {
            modEval.resize(ngrad);
        }
        std::memset(modEval.data(), 0, ngrad * sizeof(double));

        // Clear hess at start of Eval (matches legacy TModel.cc line 354)
        // This sets size to 0 so that resize(n, 0.0) in Eval* functions
        // will allocate and zero all n elements properly
        if (flag == 3) {
            hess.clear();
        }
    } else {
        if (modEval.size() < 1) {
            modEval.resize(1);
        }
    }

    modEval[0] = 1.0;

    // Check missing indicator
    if (missing_idx > -1) {
        if (data[iobs_offset + missing_idx] == 0) {
            // Clear Hessian to avoid stale data contaminating cross-component derivatives
            if (flag == 3 && hess.size() > 0) {
                std::fill(hess.begin(), hess.end(), 0.0);
            }
            return;
        }
    }
    if (ignore) {
        // Clear Hessian to avoid stale data
        if (flag == 3 && hess.size() > 0) {
            std::fill(hess.begin(), hess.end(), 0.0);
        }
        return;
    }

    // Build linear predictor(s)
    // PERFORMANCE: Use thread_local static to avoid allocation on every call
    // This is safe because Model::Eval is not called recursively
    int numlogitchoice = 2;
    static thread_local std::vector<double> expres;
    expres.assign(1, 0.0);

    if (modtype == ModelType::LOGIT && numchoice > 2) {
        expres.resize(numchoice - 1);
        std::fill(expres.begin(), expres.end(), 0.0);
        numlogitchoice = numchoice;
    }

    for (int ichoice = 0; ichoice < numlogitchoice - 1; ichoice++) {
        int ifree_local = 0;
        int nparamchoice = nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

        // Add regressor terms
        for (int i = 0; i < nregressors; i++) {
            expres[ichoice] += param[i + firstpar + ichoice*nparamchoice] * getRegValue(i, data, iobs_offset);
        }

        // Add factor and type loadings
        // IMPORTANT: Loop bound must respect actual vector sizes to avoid out-of-bounds access
        // When facnorm is provided, use its size; otherwise use numfac
        // The numtyp term is legacy and not currently used in the R interface
        int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
        for (int i = 0; i < fac_loop_bound; i++) {
            if (facnorm.size() == 0) {
                // All loadings are free
                expres[ichoice] += param[ifree_local + firstpar + ichoice*nparamchoice + nregressors] * fac[i];
                ifree_local++;
            } else {
                // Check if this loading is normalized
                if (facnorm[i] > -9998) {
                    // Fixed loading
                    expres[ichoice] += facnorm[i] * fac[i];
                } else {
                    // Free loading
                    expres[ichoice] += param[ifree_local + firstpar + ichoice*nparamchoice + nregressors] * fac[i];
                    ifree_local++;
                }
            }
        }

        // Add quadratic factor terms: sum(lambda_quad_k * f_k^2)
        // For dynamic models, skip the outcome factor
        if (n_quadratic_loadings > 0) {
            int quad_start = nregressors + ifreefac;
            int quad_idx = 0;
            for (int k = 0; k < numfac; k++) {
                if (is_dynamic && k == outcome_factor_idx) continue;  // Skip outcome factor
                double fk_sq = fac[k] * fac[k];
                expres[ichoice] += param[firstpar + ichoice*nparamchoice + quad_start + quad_idx] * fk_sq;
                quad_idx++;
            }
        }

        // Add interaction factor terms: sum(lambda_inter_jk * f_j * f_k) for j < k
        // For dynamic models, skip pairs involving the outcome factor
        if (n_interaction_loadings > 0) {
            int inter_start = nregressors + ifreefac + n_quadratic_loadings;
            int inter_idx = 0;
            for (int j = 0; j < numfac - 1; j++) {
                if (is_dynamic && j == outcome_factor_idx) continue;  // Skip outcome factor
                for (int k = j + 1; k < numfac; k++) {
                    if (is_dynamic && k == outcome_factor_idx) continue;  // Skip outcome factor
                    double fj_fk = fac[j] * fac[k];
                    expres[ichoice] += param[firstpar + ichoice*nparamchoice + inter_start + inter_idx] * fj_fk;
                    inter_idx++;
                }
            }
        }
    }

    // Add type-specific intercept to linear predictor(s)
    // For models with multiple types, this shifts the index function
    for (int ichoice = 0; ichoice < numlogitchoice - 1; ichoice++) {
        expres[ichoice] += type_intercept;
    }

    // Dispatch to model-specific evaluation
    if (modtype == ModelType::LINEAR) {
        // Z is residual: outcome - predicted value
        // For dynamic models (outcome_idx == -2), outcome is 0
        double Z = -expres[0];  // Initialize assuming outcome = 0
        if (outcome_idx > -1) Z = data[iobs_offset + outcome_idx] - expres[0];
        double sigma = std::fabs(param[firstpar + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings]);
        EvalLinear(Z, sigma, fac, param, firstpar, modEval, hess, flag, data, iobs_offset, model_free_indices);
    }
    else if (modtype == ModelType::PROBIT) {
        double obsSign = 1.0;
        if (int(data[iobs_offset + outcome_idx]) == 0) obsSign = -1.0;
        EvalProbit(expres[0], obsSign, fac, param, firstpar, modEval, hess, flag, data, iobs_offset, model_free_indices);
    }
    else if (modtype == ModelType::LOGIT) {
        double outcome_val = data[iobs_offset + outcome_idx];
        EvalLogit(expres, outcome_val, fac, param, firstpar, modEval, hess, flag, data, iobs_offset, model_free_indices);
    }
    else if (modtype == ModelType::OPROBIT) {
        int outcome_val = int(data[iobs_offset + outcome_idx]);
        EvalOprobit(expres[0], outcome_val, fac, param, firstpar, modEval, hess, flag, data, iobs_offset, model_free_indices);
    }
}

void Model::EvalLinear(double Z, double sigma, const std::vector<double>& fac,
                       const std::vector<double>& param, int firstpar,
                       std::vector<double>& modEval, std::vector<double>& hess,
                       int flag, const std::vector<double>& data, int iobs_offset,
                       const std::vector<int>* model_free_indices)
{
    // Note: model_free_indices is accepted for API consistency but not currently used
    // for optimization in EvalLinear. The main performance benefit comes from EvalLogit.
    (void)model_free_indices;  // Suppress unused warning
    // Compute likelihood
    modEval[0] = normal_pdf(Z, 0.0, sigma);

    if (flag < 2) return;

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    int npar = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings + 1; // +1 for sigma

    // Initialize Hessian if needed
    // OPTIMIZATION: Only resize if too small, avoid repeated allocation
    if (flag == 3) {
        size_t hess_size = static_cast<size_t>(npar * npar);
        if (hess.size() < hess_size) {
            hess.resize(hess_size);
        }
        double neg_inv_sigma2 = -1.0 / (sigma*sigma);
        for (int i = 0; i < npar; i++) {
            for (int j = i; j < npar; j++) {
                hess[i*npar + j] = neg_inv_sigma2;
            }
        }
    }

    // Compute gradients
    double sigma2 = sigma * sigma;
    ifreefac = 0;

    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
    for (int i = 0; i < fac_loop_bound; i++) {
        if (facnorm.size() == 0) {
            // All free
            if (i < numfac) {
                // Gradient w.r.t. factor variance
                modEval[i+1] = Z * param[ifreefac + firstpar + nregressors] / sigma2;
            }
            // Gradient w.r.t. factor loading
            modEval[1 + numfac + nregressors + ifreefac] = Z * fac[i] / sigma2;

            if (flag == 3) {
                if (i < numfac) {
                    for (int j = i; j < npar; j++)
                        hess[i*npar + j] *= param[ifreefac + firstpar + nregressors];
                    for (int j = 0; j <= i; j++)
                        hess[j*npar + i] *= param[ifreefac + firstpar + nregressors];
                }
                int index = numfac + nregressors + ifreefac;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= fac[i];
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= fac[i];
                hess[i*npar + index] += Z / sigma2;
            }
            ifreefac++;
        } else {
            // Mixed fixed/free
            if (facnorm[i] > -9998.0) {
                // Fixed loading
                if (i < numfac) {
                    modEval[i+1] = Z * facnorm[i] / sigma2;
                    if (flag == 3) {
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= facnorm[i];
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= facnorm[i];
                    }
                }
            } else {
                // Free loading
                if (i < numfac) {
                    modEval[i+1] = Z * param[ifreefac + firstpar + nregressors] / sigma2;
                }
                modEval[1 + numfac + nregressors + ifreefac] = Z * fac[i] / sigma2;

                if (flag == 3) {
                    if (i < numfac) {
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= param[ifreefac + firstpar + nregressors];
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= param[ifreefac + firstpar + nregressors];
                    }
                    int index = numfac + nregressors + ifreefac;
                    for (int j = index; j < npar; j++)
                        hess[index*npar + j] *= fac[i];
                    for (int j = 0; j <= index; j++)
                        hess[j*npar + index] *= fac[i];
                    hess[i*npar + index] += Z / sigma2;
                }
                ifreefac++;
            }
        }
    }

    // Gradients w.r.t. regression coefficients
    for (int ireg = 0; ireg < nregressors; ireg++) {
        modEval[ireg + numfac + 1] = Z * getRegValue(ireg, data, iobs_offset) / sigma2;
        if (flag == 3) {
            int index = numfac + ireg;
            for (int j = index; j < npar; j++)
                hess[index*npar + j] *= getRegValue(ireg, data, iobs_offset);
            for (int j = 0; j <= index; j++)
                hess[j*npar + index] *= getRegValue(ireg, data, iobs_offset);
        }
    }

    // Gradients w.r.t. quadratic factor loadings (lambda_quad_k)
    // d(logL)/d(lambda_quad_k) = Z * f_k^2 / sigma^2
    // For dynamic models, skip the outcome factor
    if (n_quadratic_loadings > 0) {
        int quad_grad_start = 1 + numfac + nregressors + ifreefac;
        int quad_idx = 0;
        for (int k = 0; k < numfac; k++) {
            if (is_dynamic && k == outcome_factor_idx) continue;  // Skip outcome factor

            double fk_sq = fac[k] * fac[k];
            modEval[quad_grad_start + quad_idx] = Z * fk_sq / sigma2;

            // Add quadratic contribution to factor variance gradient
            // d(xb)/d(f_k) includes 2*lambda_quad_k*f_k term
            double lambda_quad_k = param[firstpar + nregressors + ifreefac + quad_idx];
            modEval[k+1] += Z * 2.0 * lambda_quad_k * fac[k] / sigma2;

            if (flag == 3) {
                int index = numfac + nregressors + ifreefac + quad_idx;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= fk_sq;
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= fk_sq;

                // Cross-derivative d^2(logL)/(d(theta_k) d(lambda_quad_k)) = 2 * f_k / sigma^2
                hess[k*npar + index] += 2.0 * fac[k] * Z / sigma2;

                // Second derivative of xb w.r.t. f_k: d^2(xb)/d(f_k)^2 = 2*lambda_quad_k
                // This contributes Z * 2*lambda_quad_k / sigma^2 to the factor-factor Hessian diagonal
                hess[k*npar + k] += 2.0 * lambda_quad_k * Z / sigma2;
            }
            quad_idx++;
        }
    }

    // Gradients w.r.t. interaction factor loadings (lambda_inter_jk) for j < k
    // d(logL)/d(lambda_inter_jk) = Z * f_j * f_k / sigma^2
    // For dynamic models, skip pairs involving the outcome factor
    if (n_interaction_loadings > 0) {
        int inter_grad_start = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings;
        int inter_idx = 0;
        for (int j = 0; j < numfac - 1; j++) {
            if (is_dynamic && j == outcome_factor_idx) continue;  // Skip outcome factor
            for (int k = j + 1; k < numfac; k++) {
                if (is_dynamic && k == outcome_factor_idx) continue;  // Skip outcome factor

                double fj_fk = fac[j] * fac[k];
                modEval[inter_grad_start + inter_idx] = Z * fj_fk / sigma2;

                // Add interaction contribution to factor variance gradients
                // d(xb)/d(f_j) includes lambda_inter * f_k, d(xb)/d(f_k) includes lambda_inter * f_j
                double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx];
                modEval[j+1] += Z * lambda_inter * fac[k] / sigma2;
                modEval[k+1] += Z * lambda_inter * fac[j] / sigma2;

                if (flag == 3) {
                    int index = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                    for (int jj = index; jj < npar; jj++)
                        hess[index*npar + jj] *= fj_fk;
                    for (int jj = 0; jj <= index; jj++)
                        hess[jj*npar + index] *= fj_fk;

                    // Cross-derivatives d^2(logL)/(d(theta) d(lambda_inter))
                    hess[j*npar + index] += fac[k] * Z / sigma2;  // d(theta_j)/d(lambda_inter_jk)
                    hess[k*npar + index] += fac[j] * Z / sigma2;  // d(theta_k)/d(lambda_inter_jk)

                    // Cross-derivative d^2(logL)/(d(f_j) d(f_k)) from second derivative of xb
                    // d²(xb)/d(f_j)d(f_k) = lambda_inter_jk (non-zero for interaction terms)
                    // Contributes: dL/d(xb) * lambda_inter = Z/sigma^2 * lambda_inter
                    hess[j*npar + k] += lambda_inter * Z / sigma2;
                }
                inter_idx++;
            }
        }
    }

    // ===== HESSIAN CORRECTION FOR FACTOR VARIANCE BLOCK =====
    // The factor variance Hessian needs full derivative d(xb)/d(f_k), not just λ_k.
    // d(xb)/d(f_k) = λ_k + 2*λ_quad_k*f_k + Σ λ_inter_jk*f_j
    // The multiplicative structure above only used λ_k, so we add corrections.
    if (flag == 3 && (n_quadratic_loadings > 0 || n_interaction_loadings > 0)) {
        // Compute additional derivative terms for each factor
        std::vector<double> additional(numfac, 0.0);
        std::vector<double> linear_loading(numfac, 0.0);

        // Get linear loading values
        int ifree = 0;
        for (int k = 0; k < numfac; k++) {
            if (facnorm.size() == 0 || facnorm[k] <= -9998) {
                linear_loading[k] = param[firstpar + nregressors + ifree];
                ifree++;
            } else {
                linear_loading[k] = facnorm[k];
            }
        }

        // Add quadratic contributions to additional
        // For dynamic models, skip the outcome factor
        if (n_quadratic_loadings > 0) {
            int quad_idx = 0;
            for (int k = 0; k < numfac; k++) {
                if (is_dynamic && k == outcome_factor_idx) continue;  // Skip outcome factor
                double lambda_quad_k = param[firstpar + nregressors + ifreefac + quad_idx];
                additional[k] += 2.0 * lambda_quad_k * fac[k];
                quad_idx++;
            }
        }

        // Add interaction contributions to additional
        // For dynamic models, skip pairs involving the outcome factor
        if (n_interaction_loadings > 0) {
            int inter_idx = 0;
            for (int j = 0; j < numfac - 1; j++) {
                if (is_dynamic && j == outcome_factor_idx) continue;  // Skip outcome factor
                for (int k = j + 1; k < numfac; k++) {
                    if (is_dynamic && k == outcome_factor_idx) continue;  // Skip outcome factor
                    double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx];
                    additional[j] += lambda_inter * fac[k];
                    additional[k] += lambda_inter * fac[j];
                    inter_idx++;
                }
            }
        }

        // Correction for factor variance block: d²(logL)/(dθ_j)(dθ_k)
        // Current: -1/σ² * λ_j * λ_k
        // Should be: -1/σ² * (λ_j + add_j) * (λ_k + add_k)
        // Correction: -1/σ² * (λ_j*add_k + λ_k*add_j + add_j*add_k)
        for (int j = 0; j < numfac; j++) {
            for (int k = j; k < numfac; k++) {
                double correction = linear_loading[j] * additional[k] +
                                   linear_loading[k] * additional[j] +
                                   additional[j] * additional[k];
                if (j == k) {
                    // Diagonal: only add once (not double counted)
                    correction = 2.0 * linear_loading[j] * additional[j] + additional[j] * additional[j];
                }
                hess[j * npar + k] += -correction / sigma2;
            }
        }

        // Correction for cross-derivatives: d²(logL)/(dθ_k)(d other_param)
        // Current: -1/σ² * λ_k * d(xb)/d(other_param)
        // Should be: -1/σ² * (λ_k + add_k) * d(xb)/d(other_param)
        // Correction: -1/σ² * add_k * d(xb)/d(other_param)
        for (int k = 0; k < numfac; k++) {
            if (std::abs(additional[k]) < 1e-15) continue;

            // Cross with regression coefficients
            for (int ireg = 0; ireg < nregressors; ireg++) {
                double dxb_dreg = getRegValue(ireg, data, iobs_offset);
                int reg_idx = numfac + ireg;
                hess[k * npar + reg_idx] += -additional[k] * dxb_dreg / sigma2;
            }

            // Cross with free linear loadings
            int ifree2 = 0;
            for (int m = 0; m < numfac; m++) {
                if (facnorm.size() == 0 || facnorm[m] <= -9998) {
                    int load_idx = numfac + nregressors + ifree2;
                    hess[k * npar + load_idx] += -additional[k] * fac[m] / sigma2;
                    ifree2++;
                }
            }

            // Cross with quadratic loadings
            // For dynamic models, skip the outcome factor
            if (n_quadratic_loadings > 0) {
                int quad_idx = 0;
                for (int m = 0; m < numfac; m++) {
                    if (is_dynamic && m == outcome_factor_idx) continue;  // Skip outcome factor
                    int q_idx = numfac + nregressors + ifreefac + quad_idx;
                    hess[k * npar + q_idx] += -additional[k] * fac[m] * fac[m] / sigma2;
                    quad_idx++;
                }
            }

            // Cross with interaction loadings
            // For dynamic models, skip pairs involving the outcome factor
            if (n_interaction_loadings > 0) {
                int inter_idx = 0;
                for (int j2 = 0; j2 < numfac - 1; j2++) {
                    if (is_dynamic && j2 == outcome_factor_idx) continue;  // Skip outcome factor
                    for (int k2 = j2 + 1; k2 < numfac; k2++) {
                        if (is_dynamic && k2 == outcome_factor_idx) continue;  // Skip outcome factor
                        int idx = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        hess[k * npar + idx] += -additional[k] * fac[j2] * fac[k2] / sigma2;
                        inter_idx++;
                    }
                }
            }

            // Cross with sigma (will be multiplied by 2*Z/sigma later)
            int sigma_idx = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;
            hess[k * npar + sigma_idx] += -additional[k] / sigma2;
        }
    }

    // Gradient w.r.t. sigma
    int sigma_grad_idx = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;
    modEval[sigma_grad_idx] = (Z*Z / sigma - sigma) / sigma2;

    // Hessian for sigma
    if (flag == 3) {
        int sigma_index = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

        // Multiply sigma row and column by 2*Z/sigma for cross-derivatives
        for (int j = sigma_index; j < npar; j++) {
            hess[sigma_index * npar + j] *= 2.0 * Z / sigma;
        }
        for (int j = 0; j < sigma_index; j++) {
            hess[j * npar + sigma_index] *= 2.0 * Z / sigma;
        }

        // Separately handle diagonal element with correct formula
        hess[sigma_index * npar + sigma_index] = -3.0 * Z * Z / (sigma*sigma*sigma*sigma) + 1.0 / (sigma*sigma);
    }
}

void Model::EvalProbit(double expres, double obsSign, const std::vector<double>& fac,
                       const std::vector<double>& param, int firstpar,
                       std::vector<double>& modEval, std::vector<double>& hess,
                       int flag, const std::vector<double>& data, int iobs_offset,
                       const std::vector<int>* model_free_indices)
{
    // Note: model_free_indices is accepted for API consistency but not currently used
    // for optimization in EvalProbit. The main performance benefit comes from EvalLogit.
    (void)model_free_indices;  // Suppress unused warning

    // Compute likelihood
    modEval[0] = normal_cdf(obsSign * expres);

    if (flag < 2) return;

    double Z = obsSign * expres;
    double pdf = normal_pdf(obsSign * expres);
    double cdf = modEval[0];
    if (obsSign * expres < -35.0) cdf = 1.0e-50;

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    int npar = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

    // OPTIMIZATION: Only resize if too small, avoid repeated allocation
    if (flag == 3) {
        size_t hess_size = static_cast<size_t>(npar * npar);
        if (hess.size() < hess_size) {
            hess.resize(hess_size);
        }
        for (int i = 0; i < npar; i++) {
            for (int j = i; j < npar; j++) {
                hess[i*npar + j] = 1.0;
            }
        }
    }

    // Common term for probit gradients: d(logL)/d(xb) = phi/Phi * obsSign
    double dlogL_dxb = pdf * obsSign / cdf;

    // Compute gradients
    ifreefac = 0;
    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
    for (int i = 0; i < fac_loop_bound; i++) {
        if (facnorm.size() == 0) {
            if (i < numfac) {
                modEval[i+1] = pdf * obsSign * param[ifreefac + firstpar + nregressors] / cdf;
            }
            modEval[1 + numfac + nregressors + ifreefac] = pdf * (obsSign * fac[i]) / cdf;

            if (flag == 3) {
                double lambda_val = pdf / cdf;  // λ = φ/Φ
                if (i < numfac) {
                    double dZ_dtheta = obsSign * param[ifreefac + firstpar + nregressors];
                    double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                    for (int j = i; j < npar; j++)
                        hess[i*npar + j] *= row_mult;
                    for (int j = 0; j <= i; j++)
                        hess[j*npar + i] *= dZ_dtheta;
                }
                int index = numfac + nregressors + ifreefac;
                double dZ_dtheta = obsSign * fac[i];
                double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= row_mult;
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= dZ_dtheta;
            }
            ifreefac++;
        } else {
            if (facnorm[i] > -9998.0) {
                if (i < numfac) {
                    modEval[i+1] = pdf * obsSign * facnorm[i] / cdf;
                    if (flag == 3) {
                        double lambda_val = pdf / cdf;  // λ = φ/Φ
                        double dZ_dtheta = obsSign * facnorm[i];
                        double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= row_mult;
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= dZ_dtheta;
                    }
                }
            } else {
                if (i < numfac) {
                    modEval[i+1] = pdf * obsSign * param[ifreefac + firstpar + nregressors] / cdf;
                }
                modEval[1 + numfac + nregressors + ifreefac] = pdf * (obsSign * fac[i]) / cdf;

                if (flag == 3) {
                    double lambda_val = pdf / cdf;  // λ = φ/Φ
                    if (i < numfac) {
                        double dZ_dtheta = obsSign * param[ifreefac + firstpar + nregressors];
                        double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= row_mult;
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= dZ_dtheta;
                    }
                    int index = numfac + nregressors + ifreefac;
                    double dZ_dtheta_loading = obsSign * fac[i];
                    double row_mult_loading = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta_loading;
                    for (int j = index; j < npar; j++)
                        hess[index*npar + j] *= row_mult_loading;
                    for (int j = 0; j <= index; j++)
                        hess[j*npar + index] *= dZ_dtheta_loading;
                }
                ifreefac++;
            }
        }
    }

    // Gradients w.r.t. regression coefficients
    for (int ireg = 0; ireg < nregressors; ireg++) {
        modEval[ireg + numfac + 1] = pdf * (obsSign * getRegValue(ireg, data, iobs_offset)) / cdf;
        if (flag == 3) {
            double lambda_val = pdf / cdf;  // λ = φ/Φ
            int index = numfac + ireg;
            double dZ_dtheta = obsSign * getRegValue(ireg, data, iobs_offset);
            double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
            for (int j = index; j < npar; j++)
                hess[index*npar + j] *= row_mult;
            for (int j = 0; j <= index; j++)
                hess[j*npar + index] *= dZ_dtheta;
        }
    }

    // Gradients w.r.t. quadratic factor loadings (lambda_quad_k)
    // d(logL)/d(lambda_quad_k) = dlogL_dxb * f_k^2
    if (n_quadratic_loadings > 0) {
        int quad_grad_start = 1 + numfac + nregressors + ifreefac;
        for (int k = 0; k < numfac; k++) {
            double fk_sq = fac[k] * fac[k];
            modEval[quad_grad_start + k] = dlogL_dxb * fk_sq;

            // Add quadratic contribution to factor variance gradient
            double lambda_quad_k = param[firstpar + nregressors + ifreefac + k];
            modEval[k+1] += dlogL_dxb * 2.0 * lambda_quad_k * fac[k];

            if (flag == 3) {
                double lambda_val = pdf / cdf;
                int index = numfac + nregressors + ifreefac + k;
                double dZ_dtheta = obsSign * fk_sq;
                double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= row_mult;
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= dZ_dtheta;

                // Cross-derivative d(theta_k)/d(lambda_quad_k) = 2 * f_k
                hess[k*npar + index] += 2.0 * fac[k] * lambda_val * obsSign;

                // Second derivative term for factor-factor diagonal: d²(xb)/df_k² = 2*λ_q_k
                // This comes from: d²L/df² = hess_factor * (dxb/df)² + λ_mills * obsSign * d²(xb)/df²
                // The first term is in the HESSIAN CORRECTION block; this is the second term.
                hess[k*npar + k] += 2.0 * lambda_quad_k * lambda_val * obsSign;
            }
        }
    }

    // Gradients w.r.t. interaction factor loadings (lambda_inter_jk) for j < k
    // d(logL)/d(lambda_inter_jk) = dlogL_dxb * f_j * f_k
    if (n_interaction_loadings > 0) {
        int inter_grad_start = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings;
        int inter_idx = 0;
        for (int j = 0; j < numfac - 1; j++) {
            for (int k = j + 1; k < numfac; k++) {
                double fj_fk = fac[j] * fac[k];
                modEval[inter_grad_start + inter_idx] = dlogL_dxb * fj_fk;

                // Add interaction contribution to factor variance gradients
                double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx];
                modEval[j+1] += dlogL_dxb * lambda_inter * fac[k];
                modEval[k+1] += dlogL_dxb * lambda_inter * fac[j];

                if (flag == 3) {
                    double lambda_val = pdf / cdf;
                    int index = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                    double dZ_dtheta = obsSign * fj_fk;
                    double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                    for (int jj = index; jj < npar; jj++)
                        hess[index*npar + jj] *= row_mult;
                    for (int jj = 0; jj <= index; jj++)
                        hess[jj*npar + index] *= dZ_dtheta;

                    // Cross-derivatives d(theta)/d(lambda_inter)
                    hess[j*npar + index] += fac[k] * lambda_val * obsSign;
                    hess[k*npar + index] += fac[j] * lambda_val * obsSign;

                    // Second derivative term for factor-factor cross: d²(xb)/df_j df_k = λ_inter
                    // This comes from: d²L/df_j df_k = hess_factor * (dxb/df_j)*(dxb/df_k) + λ_mills * obsSign * d²(xb)/df_j df_k
                    // The first term is in the HESSIAN CORRECTION block; this is the second term.
                    hess[j*npar + k] += lambda_inter * lambda_val * obsSign;
                }
                inter_idx++;
            }
        }
    }

    // ===== HESSIAN CORRECTION FOR FACTOR VARIANCE BLOCK (PROBIT) =====
    // The factor variance Hessian needs full derivative d(xb)/d(f_k), not just λ_k.
    // d(xb)/d(f_k) = λ_k + 2*λ_quad_k*f_k + Σ λ_inter_jk*f_j
    // For probit: d²(logL)/(dθ_j)(dθ_k) = (-Z*mills - mills²) * d(xb)/d(f_j) * d(xb)/d(f_k)
    if (flag == 3 && (n_quadratic_loadings > 0 || n_interaction_loadings > 0)) {
        double lambda_val = pdf / cdf;  // Mills ratio
        double hess_factor = -Z * lambda_val - lambda_val * lambda_val;

        // Compute additional derivative terms for each factor
        std::vector<double> additional(numfac, 0.0);
        std::vector<double> linear_loading(numfac, 0.0);

        // Get linear loading values
        int ifree = 0;
        for (int k = 0; k < numfac; k++) {
            if (facnorm.size() == 0 || facnorm[k] <= -9998) {
                linear_loading[k] = param[firstpar + nregressors + ifree];
                ifree++;
            } else {
                linear_loading[k] = facnorm[k];
            }
        }

        // Add quadratic contributions to additional
        if (n_quadratic_loadings > 0) {
            for (int k = 0; k < numfac; k++) {
                double lambda_quad_k = param[firstpar + nregressors + ifreefac + k];
                additional[k] += 2.0 * lambda_quad_k * fac[k];
            }
        }

        // Add interaction contributions to additional
        if (n_interaction_loadings > 0) {
            int inter_idx = 0;
            for (int j = 0; j < numfac - 1; j++) {
                for (int k = j + 1; k < numfac; k++) {
                    double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx];
                    additional[j] += lambda_inter * fac[k];
                    additional[k] += lambda_inter * fac[j];
                    inter_idx++;
                }
            }
        }

        // Correction for factor variance block: d²(logL)/(dθ_j)(dθ_k)
        for (int j = 0; j < numfac; j++) {
            for (int k = j; k < numfac; k++) {
                double correction = linear_loading[j] * additional[k] +
                                   linear_loading[k] * additional[j] +
                                   additional[j] * additional[k];
                if (j == k) {
                    correction = 2.0 * linear_loading[j] * additional[j] + additional[j] * additional[j];
                }
                hess[j * npar + k] += hess_factor * correction;
            }
        }

        // Correction for cross-derivatives: d²(logL)/(dθ_k)(d other_param)
        // Current: hess[k,j] = hess_factor * λ_k * d(xb)/d(j)  (after obsSign² = 1 cancellation)
        // Target:  hess[k,j] = hess_factor * (λ_k + add_k) * d(xb)/d(j)
        // Correction: hess_factor * add_k * d(xb)/d(j)
        for (int k = 0; k < numfac; k++) {
            if (std::abs(additional[k]) < 1e-15) continue;

            // Cross with regression coefficients
            for (int ireg = 0; ireg < nregressors; ireg++) {
                double dxb_dreg = getRegValue(ireg, data, iobs_offset);
                int reg_idx = numfac + ireg;
                hess[k * npar + reg_idx] += hess_factor * additional[k] * dxb_dreg;
            }

            // Cross with free linear loadings
            int ifree2 = 0;
            for (int m = 0; m < numfac; m++) {
                if (facnorm.size() == 0 || facnorm[m] <= -9998) {
                    int load_idx = numfac + nregressors + ifree2;
                    hess[k * npar + load_idx] += hess_factor * additional[k] * fac[m];
                    ifree2++;
                }
            }

            // Cross with quadratic loadings
            if (n_quadratic_loadings > 0) {
                for (int m = 0; m < numfac; m++) {
                    int quad_idx = numfac + nregressors + ifreefac + m;
                    hess[k * npar + quad_idx] += hess_factor * additional[k] * fac[m] * fac[m];
                }
            }

            // Cross with interaction loadings
            if (n_interaction_loadings > 0) {
                int inter_idx = 0;
                for (int j2 = 0; j2 < numfac - 1; j2++) {
                    for (int k2 = j2 + 1; k2 < numfac; k2++) {
                        int idx = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        hess[k * npar + idx] += hess_factor * additional[k] * fac[j2] * fac[k2];
                        inter_idx++;
                    }
                }
            }
        }
    }

    // Add cross-derivative terms for Hessian
    // d²(log L)/df_k dλ_k = (-Z*λ - λ²) * λ_k * f_k + λ * s
    // where λ = φ/Φ (Mills ratio), s = obsSign
    // The first term comes from the multiplicative Hessian structure,
    // but we need to add the second term λ * s = (pdf/cdf) * obsSign
    if (flag == 3) {
        double lambda_val = pdf / cdf;  // Mills ratio
        ifreefac = 0;
        for (int i = 0; i < numfac; i++) {
            if (facnorm.size() == 0 || facnorm[i] <= -9998) {
                int index = numfac + nregressors + ifreefac;
                hess[i*npar + index] += lambda_val * obsSign;
                ifreefac++;
            }
        }
    }
}

void Model::EvalLogit(const std::vector<double>& expres, double outcome,
                      const std::vector<double>& fac,
                      const std::vector<double>& param, int firstpar,
                      std::vector<double>& modEval, std::vector<double>& hess,
                      int flag, const std::vector<double>& data, int iobs_offset,
                      const std::vector<int>* model_free_indices)
{
    // Multinomial logit with K choices (numchoice)
    // Choice 0 is reference category with Z_0 = 0
    // For choices 1 to K-1, we have separate parameters
    //
    // For exploded logit (numrank > 1), we have multiple ranked choices per observation.
    //
    // OPTIMIZATION: If model_free_indices is provided, only compute Hessian entries
    // for free model parameters. This significantly speeds up computation when many
    // coefficients are fixed (e.g., via fix_coefficient()).
    // The likelihood is the product of per-rank likelihoods.
    // Gradients are the sum of per-rank log-likelihood gradients.

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    // OPTIMIZATION: Precompute cumulative count of free loadings for each factor index
    // free_loading_cumcount[ifac] = number of free loadings with index < ifac
    // This avoids O(numfac) count_if calls inside O(numchoice) loops
    // PERFORMANCE: Use thread_local to avoid repeated allocation across calls
    static thread_local std::vector<int> free_loading_cumcount;
    if (free_loading_cumcount.size() < static_cast<size_t>(numfac + 1)) {
        free_loading_cumcount.resize(numfac + 1);
    }
    if (facnorm.size() > 0) {
        free_loading_cumcount[0] = 0;
        for (int i = 0; i < numfac; i++) {
            free_loading_cumcount[i + 1] = free_loading_cumcount[i] + (facnorm[i] <= -9998 ? 1 : 0);
        }
    } else {
        // All loadings are free
        for (int i = 0; i <= numfac; i++) {
            free_loading_cumcount[i] = i;
        }
    }

    // Number of parameters per choice
    int nparamchoice = nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

    // Initialize total likelihood to 1.0 (will multiply per-rank likelihoods)
    modEval[0] = 1.0;

    // Number of parameters for gradient/Hessian sizing
    int npar = numfac + (numchoice - 1) * nparamchoice;

    // For exploded logit: track which choices have been made (excluded from later ranks)
    // PERFORMANCE: Use thread_local to avoid repeated allocation across calls
    static thread_local std::vector<bool> chosen;
    if (chosen.size() < static_cast<size_t>(numchoice)) {
        chosen.resize(numchoice);
    }
    std::fill(chosen.begin(), chosen.begin() + numchoice, false);

    // Intermediate gradient storage for Hessian (allocated once, reset per rank)
    // PERFORMANCE: Use thread_local to avoid repeated allocation across calls
    static thread_local std::vector<double> logitgrad;

    // OPTIMIZATION: Build boolean vector marking which Hessian indices are free
    // Factor indices (0..numfac-1) are always needed; model params checked against free list
    // This allows skipping Hessian computation for fixed parameters
    // PERFORMANCE: Use thread_local to avoid repeated allocation across calls
    static thread_local std::vector<bool> hess_idx_free;
    bool use_free_opt = (model_free_indices != nullptr && flag == 3);

    if (flag == 3) {
        size_t logitgrad_size = static_cast<size_t>(numchoice * npar);
        if (logitgrad.size() < logitgrad_size) {
            logitgrad.resize(logitgrad_size);
        }
        // Note: logitgrad is zeroed inside the rank loop at line ~1088, not here

        // Initialize Hessian BEFORE the rank loop so it accumulates across ranks
        // OPTIMIZATION: Use resize like legacy TModel.cc (line 1149) instead of memset
        // resize(n, 0.0) only zeros NEW elements when growing, does nothing if size matches
        // This matches legacy behavior and is much faster for repeated calls with same-size buffer
        // Note: This means hess may contain stale values from previous calls, which is
        // acceptable because the Hessian is accumulated across all observations anyway
        size_t hess_size = static_cast<size_t>(npar * npar);
        hess.resize(hess_size, 0.0);

        // Build free index lookup for optimization
        if (use_free_opt) {
            if (hess_idx_free.size() < static_cast<size_t>(npar)) {
                hess_idx_free.resize(npar);
            }
            std::fill(hess_idx_free.begin(), hess_idx_free.begin() + npar, false);
            // Factor indices are always free
            for (int i = 0; i < numfac; i++) {
                hess_idx_free[i] = true;
            }
            // Model parameters: check against free list
            for (int fi : *model_free_indices) {
                int hess_idx = numfac + fi;
                if (hess_idx < npar) {
                    hess_idx_free[hess_idx] = true;
                }
            }
        }
    }

    // OPTIMIZATION: Pre-allocate vectors outside the rank loop to avoid repeated allocations
    // PERFORMANCE: Use thread_local to avoid repeated allocation across calls to EvalLogit
    static thread_local std::vector<double> rankedChoiceCorr;
    static thread_local std::vector<double> pdf;
    if (rankedChoiceCorr.size() < static_cast<size_t>(numchoice - 1)) {
        rankedChoiceCorr.resize(numchoice - 1);
    }
    if (pdf.size() < static_cast<size_t>(numchoice)) {
        pdf.resize(numchoice);
    }

    // Loop over ranks for exploded logit
    for (int irank = 0; irank < numrank; irank++) {
        // Get outcome for this rank
        double rank_outcome;
        if (numrank == 1) {
            rank_outcome = outcome;  // Standard logit: use passed outcome
        } else {
            // Exploded logit: read from data using stored outcome indices
            rank_outcome = data[iobs_offset + outcome_indices[irank]];
        }

        // Observed choice (0-indexed: 0, 1, ..., numchoice-1)
        int obsCat = int(rank_outcome) - 1;  // Convert from 1-indexed to 0-indexed

        // Skip invalid ranks (missing choices indicated by value not in 1..numchoice)
        if (obsCat < 0 || obsCat >= numchoice) {
            if (numrank == 1) {
                // For standard logit, invalid choice is an error
                Rf_warning("Invalid multinomial choice %d (must be 1 to %d)",
                           int(rank_outcome), numchoice);
                modEval[0] = 1e-100;
                return;
            }
            // For exploded logit, skip this rank (individual didn't use all ranks)
            continue;
        }

        // Load rank-share corrections for this rank (if ranksharevar is provided)
        // Layout: ranksharevar_idx + (numchoice-1)*irank + icat for icat=0..numchoice-2
        // Reset rankedChoiceCorr to 0 for this rank
        std::fill(rankedChoiceCorr.begin(), rankedChoiceCorr.end(), 0.0);
        if (ranksharevar_idx >= 0) {
            for (int icat = 0; icat < numchoice - 1; icat++) {
                int idx = ranksharevar_idx + (numchoice - 1) * irank + icat;
                double corr_val = data[iobs_offset + idx];
                if (corr_val > -9998.0) {
                    rankedChoiceCorr[icat] = corr_val;
                }
            }
        }

        // Compute conditional denominator
        // If exclude_chosen=true (standard exploded logit): exclude already-chosen alternatives
        // If exclude_chosen=false (nested logit): all alternatives remain available
        double logitdenom = 0.0;
        for (int icat = 0; icat < numchoice; icat++) {
            if (!exclude_chosen || !chosen[icat]) {  // Include if not excluding OR not already chosen
                if (icat == 0) {
                    logitdenom += 1.0;  // Reference category (no rankshare correction)
                } else {
                    logitdenom += std::exp(expres[icat - 1] + rankedChoiceCorr[icat - 1]);
                }
            }
        }

        // Compute conditional probabilities for available choices
        // Reset pdf to 0 for this rank (vector pre-allocated outside loop)
        std::fill(pdf.begin(), pdf.end(), 0.0);
        for (int icat = 0; icat < numchoice; icat++) {
            if (!exclude_chosen || !chosen[icat]) {
                if (icat == 0) {
                    pdf[icat] = 1.0 / logitdenom;
                } else {
                    pdf[icat] = std::exp(expres[icat - 1] + rankedChoiceCorr[icat - 1]) / logitdenom;
                }
            }
            // pdf[icat] = 0 for excluded alternatives
        }

        // Per-rank likelihood: P(Y = obsCat | available choices)
        double dens = pdf[obsCat];

        // Multiply into total likelihood
        modEval[0] *= dens;

        // Mark this choice as taken for subsequent ranks (only if excluding)
        if (exclude_chosen) {
            chosen[obsCat] = true;
        }

        if (flag < 2) continue;  // Skip gradient computation for this rank

        // ===== GRADIENT CALCULATION =====
        // Reset logitgrad for this rank (it stores rank-specific demeaned derivatives for Hessian)
        // Use memset for speed and only zero the elements we need (numchoice * npar)
        if (flag == 3) {
            std::memset(logitgrad.data(), 0, static_cast<size_t>(numchoice * npar) * sizeof(double));
        }

        // Gradient for factor variance parameters (theta)
    for (int ifac = 0; ifac < numfac; ifac++) {
        // Get loading value for this factor (either free from param or fixed from facnorm)
        double loading_for_grad;
        bool is_free_loading = (facnorm.size() == 0 || facnorm[ifac] <= -9998);

        if (is_free_loading) {
            // Free loading - will get from parameter vector per choice
        } else {
            // Fixed loading - use the fixed value
            loading_for_grad = facnorm[ifac];
        }

        // Term from observed choice (gradient only - depends on observed category)
        if (obsCat > 0) {
            if (is_free_loading) {
                int param_idx = firstpar + (obsCat - 1) * nparamchoice + nregressors;
                int loading_idx = param_idx + free_loading_cumcount[ifac];
                modEval[1 + ifac] += param[loading_idx];
            } else {
                // Fixed loading - use fixed value for gradient
                modEval[1 + ifac] += loading_for_grad;
            }
        }

        // logitgrad stores dZ_jcat/df_ifac - E[dZ/df_ifac] for Hessian computation
        // The dZ_jcat/df_ifac = lambda_jcat,ifac term is INDEPENDENT of observed category
        // so it must be computed outside the obsCat check
        if (flag == 3) {
            if (is_free_loading) {
                int free_idx_for_ifac = free_loading_cumcount[ifac];
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jparam_idx = firstpar + (jcat - 1) * nparamchoice + nregressors;
                    int jloading_idx = jparam_idx + free_idx_for_ifac;
                    logitgrad[jcat * npar + ifac] += param[jloading_idx];
                }
            } else {
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    logitgrad[jcat * npar + ifac] += loading_for_grad;
                }
            }
        }

        // Sum over all non-reference choices
        // Precompute free loading index offset for this ifac (constant across icat)
        int free_idx_offset = free_loading_cumcount[ifac];
        for (int icat = 1; icat < numchoice; icat++) {
            if (is_free_loading) {
                int param_idx = firstpar + (icat - 1) * nparamchoice + nregressors;
                int loading_idx = param_idx + free_idx_offset;
                modEval[1 + ifac] += -pdf[icat] * param[loading_idx];

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + ifac] += -pdf[icat] * param[loading_idx];
                    }
                }
            } else {
                // Fixed loading - use fixed value for gradient
                modEval[1 + ifac] += -pdf[icat] * loading_for_grad;

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + ifac] += -pdf[icat] * loading_for_grad;
                    }
                }
            }
        }
    }

    // Gradient for factor loadings and regression coefficients
    // Following legacy TModel.cc pattern exactly

    // Factor loadings - obsCat term and logitgrad
    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
    int ifree = 0;
    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
            double fval = fac[ifac];

            // obsCat term for gradient
            if (obsCat > 0) {
                int base_idx = numfac + (obsCat - 1) * nparamchoice;
                modEval[1 + base_idx + nregressors + ifree] += fval;
            }

            // logitgrad update (unconditional, not nested in obsCat check!)
            if (flag == 3) {
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                    logitgrad[jcat * npar + jbase_idx + nregressors + ifree] += fval;
                }
            }

            // All categories term
            for (int icat = 1; icat < numchoice; icat++) {
                int base_idx = numfac + (icat - 1) * nparamchoice;
                modEval[1 + base_idx + nregressors + ifree] += -pdf[icat] * fval;

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + base_idx + nregressors + ifree] += -pdf[icat] * fval;
                    }
                }
            }

            ifree++;
        }
    }

    // Regression coefficients - obsCat term and logitgrad
    for (int ireg = 0; ireg < nregressors; ireg++) {
        double xval = getRegValue(ireg, data, iobs_offset);

        // obsCat term for gradient
        if (obsCat > 0) {
            int base_idx = numfac + (obsCat - 1) * nparamchoice;
            modEval[1 + base_idx + ireg] += xval;
        }

        // logitgrad update (unconditional, not nested in obsCat check!)
        if (flag == 3) {
            for (int jcat = 1; jcat < numchoice; jcat++) {
                int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                logitgrad[jcat * npar + jbase_idx + ireg] += xval;
            }
        }

        // All categories term
        for (int icat = 1; icat < numchoice; icat++) {
            int base_idx = numfac + (icat - 1) * nparamchoice;
            modEval[1 + base_idx + ireg] += -pdf[icat] * xval;

            if (flag == 3) {
                for (int jcat = 0; jcat < numchoice; jcat++) {
                    logitgrad[jcat * npar + base_idx + ireg] += -pdf[icat] * xval;
                }
            }
        }
    }

    // Quadratic factor loadings - obsCat term and logitgrad
    if (n_quadratic_loadings > 0) {
        for (int k = 0; k < numfac; k++) {
            double fk_sq = fac[k] * fac[k];
            int quad_offset = nregressors + ifreefac + k;

            // obsCat term for gradient
            if (obsCat > 0) {
                int base_idx = numfac + (obsCat - 1) * nparamchoice;
                modEval[1 + base_idx + quad_offset] += fk_sq;

                // Add quadratic contribution to factor variance gradient
                double lambda_quad_k = param[firstpar + (obsCat - 1) * nparamchoice + quad_offset];
                modEval[1 + k] += 2.0 * lambda_quad_k * fac[k];
            }

            // logitgrad update
            if (flag == 3) {
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                    logitgrad[jcat * npar + jbase_idx + quad_offset] += fk_sq;

                    // Add quadratic contribution to factor logitgrad (dZ_jcat/df_k includes 2*lambda_quad_k*f_k)
                    double lambda_quad_jcat = param[firstpar + (jcat - 1) * nparamchoice + quad_offset];
                    logitgrad[jcat * npar + k] += 2.0 * lambda_quad_jcat * fac[k];
                }
            }

            // All categories term
            for (int icat = 1; icat < numchoice; icat++) {
                int base_idx = numfac + (icat - 1) * nparamchoice;
                modEval[1 + base_idx + quad_offset] += -pdf[icat] * fk_sq;

                // Add quadratic contribution to factor variance gradient (all categories)
                double lambda_quad_k = param[firstpar + (icat - 1) * nparamchoice + quad_offset];
                modEval[1 + k] += -pdf[icat] * 2.0 * lambda_quad_k * fac[k];

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + base_idx + quad_offset] += -pdf[icat] * fk_sq;
                        // Update factor logitgrad with quadratic contribution (was missing!)
                        logitgrad[jcat * npar + k] += -pdf[icat] * 2.0 * lambda_quad_k * fac[k];
                    }
                }
            }
        }
    }

    // Interaction factor loadings - obsCat term and logitgrad
    if (n_interaction_loadings > 0) {
        int inter_idx = 0;
        for (int j = 0; j < numfac - 1; j++) {
            for (int k = j + 1; k < numfac; k++) {
                double fj_fk = fac[j] * fac[k];
                int inter_offset = nregressors + ifreefac + n_quadratic_loadings + inter_idx;

                // obsCat term for gradient
                if (obsCat > 0) {
                    int base_idx = numfac + (obsCat - 1) * nparamchoice;
                    modEval[1 + base_idx + inter_offset] += fj_fk;

                    // Add interaction contribution to factor variance gradients
                    double lambda_inter = param[firstpar + (obsCat - 1) * nparamchoice + inter_offset];
                    modEval[1 + j] += lambda_inter * fac[k];
                    modEval[1 + k] += lambda_inter * fac[j];
                }

                // logitgrad update
                if (flag == 3) {
                    for (int jcat = 1; jcat < numchoice; jcat++) {
                        int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                        logitgrad[jcat * npar + jbase_idx + inter_offset] += fj_fk;

                        // Add interaction contribution to factor logitgrad (dZ_jcat/df_j includes lambda_inter*f_k, etc.)
                        double lambda_inter_jcat = param[firstpar + (jcat - 1) * nparamchoice + inter_offset];
                        logitgrad[jcat * npar + j] += lambda_inter_jcat * fac[k];
                        logitgrad[jcat * npar + k] += lambda_inter_jcat * fac[j];
                    }
                }

                // All categories term
                for (int icat = 1; icat < numchoice; icat++) {
                    int base_idx = numfac + (icat - 1) * nparamchoice;
                    modEval[1 + base_idx + inter_offset] += -pdf[icat] * fj_fk;

                    // Add interaction contribution to factor variance gradients (all categories)
                    double lambda_inter = param[firstpar + (icat - 1) * nparamchoice + inter_offset];
                    modEval[1 + j] += -pdf[icat] * lambda_inter * fac[k];
                    modEval[1 + k] += -pdf[icat] * lambda_inter * fac[j];

                    if (flag == 3) {
                        for (int jcat_inner = 0; jcat_inner < numchoice; jcat_inner++) {
                            logitgrad[jcat_inner * npar + base_idx + inter_offset] += -pdf[icat] * fj_fk;
                            // Update factor logitgrad with interaction contribution (was missing!)
                            logitgrad[jcat_inner * npar + j] += -pdf[icat] * lambda_inter * fac[k];
                            logitgrad[jcat_inner * npar + k] += -pdf[icat] * lambda_inter * fac[j];
                        }
                    }
                }
                inter_idx++;
            }
        }
    }

    // ===== HESSIAN CALCULATION =====
    // NOTE: Hessian is initialized BEFORE the rank loop (lines ~900-910)
    // and accumulates across ranks. logitgrad is reset per rank.
    if (flag == 3) {
        // OPTIMIZATION: Detect if we can use the fast path (standard logit without complexity)
        // Fast path when: all loadings are free, no quadratic/interaction terms
        // Note: ifreefac == numfac means all factor loadings are free (either facnorm.size()==0
        // or all entries in facnorm are <= -9998)
        // Note: We no longer exclude fixed coefficients (use_free_opt) from fast path because
        // the overhead of checking hess_idx_free on every iteration is worse than just computing
        // the full Hessian. The fixed parameter handling is done at the FactorModel aggregation level.
        bool all_loadings_free = (facnorm.size() == 0) || (ifreefac == numfac);
        bool use_fast_path = all_loadings_free && (n_quadratic_loadings == 0) &&
                             (n_interaction_loadings == 0);

        if (use_fast_path) {
            // ===== FAST PATH: Matches legacy TModel.cc structure exactly =====
            // Access data directly like legacy code - no caching overhead

            // Second-order derivative terms: dZ/dtheta dalpha
            if (obsCat > 0) {
                for (int ifac = 0; ifac < numfac; ifac++) {
                    int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifac;
                    hess[ifac * npar + index] += 1.0;
                }
            }

            // Second term of second-order derivatives
            for (int icat = 1; icat < numchoice; icat++) {
                for (int ifac = 0; ifac < numfac; ifac++) {
                    int index = numfac + (icat - 1) * nparamchoice + nregressors + ifac;
                    hess[ifac * npar + index] += -pdf[icat];
                }
            }

            // First-order derivative Hessian terms: dtheta d...
            for (int ifac = 0; ifac < numfac; ifac++) {
                for (int icat = 1; icat < numchoice; icat++) {
                    double loading_icat = param[firstpar + (icat - 1) * nparamchoice + nregressors + ifac];

                    // dtheta dtheta
                    for (int jfac = ifac; jfac < numfac; jfac++) {
                        hess[ifac * npar + jfac] += -pdf[icat] * logitgrad[icat * npar + jfac] * loading_icat;
                    }

                    // dtheta dalpha and dtheta dbeta for each jcat
                    for (int jcat = 1; jcat < numchoice; jcat++) {
                        int jbase_idx = numfac + (jcat - 1) * nparamchoice;

                        // dtheta dalpha
                        for (int jfac = 0; jfac < numfac; jfac++) {
                            int index = jbase_idx + nregressors + jfac;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_icat;
                        }

                        // dtheta dbeta
                        for (int jreg = 0; jreg < nregressors; jreg++) {
                            int index = jbase_idx + jreg;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_icat;
                        }
                    }
                }
            }

            // Loop over row categories for remaining Hessian terms
            for (int icat = 1; icat < numchoice; icat++) {
                int ibase_idx = numfac + (icat - 1) * nparamchoice;

                for (int jcat = icat; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;

                    // dbeta dbeta
                    for (int ireg = 0; ireg < nregressors; ireg++) {
                        int index1 = ibase_idx + ireg;
                        for (int jreg = 0; jreg < nregressors; jreg++) {
                            if ((jcat > icat) || (jreg >= ireg)) {
                                int index2 = jbase_idx + jreg;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * data[iobs_offset + regressors[ireg]];
                            }
                        }
                    }

                    // dbeta dalpha
                    for (int ireg = 0; ireg < nregressors; ireg++) {
                        int index1 = ibase_idx + ireg;
                        for (int jfac = 0; jfac < numfac; jfac++) {
                            int index2 = jbase_idx + nregressors + jfac;
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * data[iobs_offset + regressors[ireg]];
                        }
                    }

                    // dalpha dbeta (only for jcat > icat)
                    if (jcat > icat) {
                        for (int ifac = 0; ifac < numfac; ifac++) {
                            int index1 = ibase_idx + nregressors + ifac;
                            for (int jreg = 0; jreg < nregressors; jreg++) {
                                int index2 = jbase_idx + jreg;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                        }
                    }

                    // dalpha dalpha
                    for (int ifac = 0; ifac < numfac; ifac++) {
                        int index1 = ibase_idx + nregressors + ifac;
                        for (int jfac = 0; jfac < numfac; jfac++) {
                            if ((jcat > icat) || (jfac >= ifac)) {
                                int index2 = jbase_idx + nregressors + jfac;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                        }
                    }
                }
            }
        } else {
            // ===== SLOW PATH: Full implementation with all features =====
            // Second-order derivative terms: dZ/dtheta dalpha
            if (obsCat > 0) {
                int ifree = 0;
                for (int ifac = 0; ifac < numfac; ifac++) {
                    if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                        int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifree;
                        hess[ifac * npar + index] += 1.0;
                        ifree++;
                    }
                }
            }

            // Second term of second-order derivatives
            for (int icat = 1; icat < numchoice; icat++) {
                int ifree = 0;
                for (int ifac = 0; ifac < numfac; ifac++) {
                    if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                        int index = numfac + (icat - 1) * nparamchoice + nregressors + ifree;
                        hess[ifac * npar + index] += -pdf[icat];
                        ifree++;
                    }
                }
            }

        // Second-order derivative terms for d(theta)/d(lambda_quad): d^2Z/dtheta_k dlambda_quad_k = 2*f_k
        if (n_quadratic_loadings > 0) {
            if (obsCat > 0) {
                for (int k = 0; k < numfac; k++) {
                    int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifreefac + k;
                    hess[k * npar + index] += 2.0 * fac[k];
                }
            }
            for (int icat = 1; icat < numchoice; icat++) {
                for (int k = 0; k < numfac; k++) {
                    int index = numfac + (icat - 1) * nparamchoice + nregressors + ifreefac + k;
                    hess[k * npar + index] += -pdf[icat] * 2.0 * fac[k];
                }
            }

            // Second derivative of linear predictor w.r.t. factor: d^2Z/dtheta_k^2 = 2*lambda_quad_k
            // This contributes to the factor-factor diagonal Hessian
            if (obsCat > 0) {
                for (int k = 0; k < numfac; k++) {
                    double lambda_quad_k = param[firstpar + (obsCat - 1) * nparamchoice + nregressors + ifreefac + k];
                    hess[k * npar + k] += 2.0 * lambda_quad_k;
                }
            }
            for (int icat = 1; icat < numchoice; icat++) {
                for (int k = 0; k < numfac; k++) {
                    double lambda_quad_k = param[firstpar + (icat - 1) * nparamchoice + nregressors + ifreefac + k];
                    hess[k * npar + k] += -pdf[icat] * 2.0 * lambda_quad_k;
                }
            }
        }

        // Second-order derivative terms for d(theta)/d(lambda_inter): d^2Z/dtheta_j dlambda_inter_jk = f_k (and similar for k)
        if (n_interaction_loadings > 0) {
            if (obsCat > 0) {
                int inter_idx = 0;
                for (int j = 0; j < numfac - 1; j++) {
                    for (int k = j + 1; k < numfac; k++) {
                        int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        hess[j * npar + index] += fac[k];  // d(theta_j)/d(lambda_inter_jk)
                        hess[k * npar + index] += fac[j];  // d(theta_k)/d(lambda_inter_jk)
                        inter_idx++;
                    }
                }
            }
            for (int icat = 1; icat < numchoice; icat++) {
                int inter_idx = 0;
                for (int j = 0; j < numfac - 1; j++) {
                    for (int k = j + 1; k < numfac; k++) {
                        int index = numfac + (icat - 1) * nparamchoice + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        hess[j * npar + index] += -pdf[icat] * fac[k];
                        hess[k * npar + index] += -pdf[icat] * fac[j];
                        inter_idx++;
                    }
                }
            }

            // Second derivative of linear predictor w.r.t. factors: d^2Z/df_j df_k = lambda_inter
            // This contributes to the factor-factor cross Hessian (was missing!)
            if (obsCat > 0) {
                int inter_idx = 0;
                for (int jj = 0; jj < numfac - 1; jj++) {
                    for (int kk = jj + 1; kk < numfac; kk++) {
                        int inter_offset = nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        double lambda_inter = param[firstpar + (obsCat - 1) * nparamchoice + inter_offset];
                        hess[jj * npar + kk] += lambda_inter;
                        inter_idx++;
                    }
                }
            }
            for (int icat = 1; icat < numchoice; icat++) {
                int inter_idx = 0;
                for (int jj = 0; jj < numfac - 1; jj++) {
                    for (int kk = jj + 1; kk < numfac; kk++) {
                        int inter_offset = nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        double lambda_inter = param[firstpar + (icat - 1) * nparamchoice + inter_offset];
                        hess[jj * npar + kk] += -pdf[icat] * lambda_inter;
                        inter_idx++;
                    }
                }
            }
        }

        // First-order derivative Hessian terms
        // dtheta d...
        for (int ifac = 0; ifac < numfac; ifac++) {
            // Precompute whether this factor has a free loading and its index offset
            bool ifac_is_free = (facnorm.size() == 0 || facnorm[ifac] <= -9998);
            int ifac_free_offset = free_loading_cumcount[ifac];
            double ifac_fixed_loading = ifac_is_free ? 0.0 : facnorm[ifac];

            for (int icat = 1; icat < numchoice; icat++) {
                int param_idx = firstpar + (icat - 1) * nparamchoice + nregressors;
                // Compute FULL derivative dZ_icat/df_ifac (not just linear loading)
                double loading_val = 0.0;

                if (ifac_is_free) {
                    // Free loading - get from parameter vector
                    int loading_idx = param_idx + ifac_free_offset;
                    loading_val = param[loading_idx];
                } else {
                    // Fixed loading - use the fixed value from facnorm
                    loading_val = ifac_fixed_loading;
                }

                // Add quadratic contribution: dZ/df_ifac includes 2*lambda_quad_ifac*f_ifac
                if (n_quadratic_loadings > 0) {
                    int quad_offset = nregressors + ifreefac + ifac;
                    double lambda_quad_ifac = param[firstpar + (icat - 1) * nparamchoice + quad_offset];
                    loading_val += 2.0 * lambda_quad_ifac * fac[ifac];
                }

                // Add interaction contribution: dZ/df_ifac includes lambda_inter*f_other
                if (n_interaction_loadings > 0) {
                    int inter_idx = 0;
                    for (int jj = 0; jj < numfac - 1; jj++) {
                        for (int kk = jj + 1; kk < numfac; kk++) {
                            int inter_offset = nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                            double lambda_inter = param[firstpar + (icat - 1) * nparamchoice + inter_offset];
                            if (ifac == jj) {
                                loading_val += lambda_inter * fac[kk];
                            } else if (ifac == kk) {
                                loading_val += lambda_inter * fac[jj];
                            }
                            inter_idx++;
                        }
                    }
                }

                // NOTE: dtheta dtheta is computed separately below (outside icat loop)
                // to sum over ALL categories including reference

                // dtheta dalpha, dtheta dbeta, dtheta dquad, dtheta dinter for each category
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;

                    // dtheta dbeta
                    for (int jreg = 0; jreg < nregressors; jreg++) {
                        int index = jbase_idx + jreg;
                        hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                    }

                    // dtheta dalpha
                    int jfree = 0;
                    for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                        if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                            int index = jbase_idx + nregressors + jfree;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                            jfree++;
                        }
                    }

                    // dtheta dquad
                    if (n_quadratic_loadings > 0) {
                        for (int k = 0; k < numfac; k++) {
                            int index = jbase_idx + nregressors + ifreefac + k;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                        }
                    }

                    // dtheta dinter
                    if (n_interaction_loadings > 0) {
                        for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                            int index = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                        }
                    }
                }
            }
        }

        // dtheta dtheta: Factor variance Hessian
        // H_ij = -sum_k P_k * (dZ_k/df_i - μ_i) * (dZ_k/df_j - μ_j)
        // Uses logitgrad which stores the demeaned derivatives
        // MUST sum over ALL categories including reference (k=0)
        for (int ifac = 0; ifac < numfac; ifac++) {
            for (int jfac = ifac; jfac < numfac; jfac++) {
                for (int kcat = 0; kcat < numchoice; kcat++) {
                    hess[ifac * npar + jfac] += -pdf[kcat] * logitgrad[kcat * npar + ifac] * logitgrad[kcat * npar + jfac];
                }
            }
        }

        // Loop over row categories for remaining Hessian terms
        // OPTIMIZATION: use_free_opt allows skipping Hessian entries for fixed parameters
        for (int icat = 1; icat < numchoice; icat++) {
            int ibase_idx = numfac + (icat - 1) * nparamchoice;

            for (int jcat = icat; jcat < numchoice; jcat++) {
                int jbase_idx = numfac + (jcat - 1) * nparamchoice;

                // dbeta dbeta
                for (int ireg = 0; ireg < nregressors; ireg++) {
                    int index1 = ibase_idx + ireg;
                    if (use_free_opt && !hess_idx_free[index1]) continue;  // Skip fixed row
                    double xval_i = getRegValue(ireg, data, iobs_offset);
                    for (int jreg = 0; jreg < nregressors; jreg++) {
                        if ((jcat > icat) || (jreg >= ireg)) {
                            int index2 = jbase_idx + jreg;
                            if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                        }
                    }
                }

                // dbeta dalpha
                for (int ireg = 0; ireg < nregressors; ireg++) {
                    int index1 = ibase_idx + ireg;
                    if (use_free_opt && !hess_idx_free[index1]) continue;  // Skip fixed row
                    double xval_i = getRegValue(ireg, data, iobs_offset);
                    int jfree = 0;
                    for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                        if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                            int index2 = jbase_idx + nregressors + jfree;
                            if (!use_free_opt || hess_idx_free[index2]) {  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                            }
                            jfree++;
                        }
                    }
                }

                // dalpha dbeta (only for jcat > icat)
                if (jcat > icat) {
                    int ifree = 0;
                    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                            int index1 = ibase_idx + nregressors + ifree;
                            if (use_free_opt && !hess_idx_free[index1]) { ifree++; continue; }  // Skip fixed row
                            for (int jreg = 0; jreg < nregressors; jreg++) {
                                int index2 = jbase_idx + jreg;
                                if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                            ifree++;
                        }
                    }
                }

                // dalpha dalpha
                int ifree = 0;
                for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                    if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                        int index1 = ibase_idx + nregressors + ifree;
                        if (use_free_opt && !hess_idx_free[index1]) { ifree++; continue; }  // Skip fixed row
                        int jfree = 0;
                        for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                            if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                                if ((jcat > icat) || (jfac >= ifac)) {
                                    int index2 = jbase_idx + nregressors + jfree;
                                    if (!use_free_opt || hess_idx_free[index2]) {  // Skip fixed col
                                        hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                                    }
                                }
                                jfree++;
                            }
                        }
                        ifree++;
                    }
                }

                // dbeta dquad
                if (n_quadratic_loadings > 0) {
                    for (int ireg = 0; ireg < nregressors; ireg++) {
                        int index1 = ibase_idx + ireg;
                        if (use_free_opt && !hess_idx_free[index1]) continue;  // Skip fixed row
                        double xval_i = getRegValue(ireg, data, iobs_offset);
                        for (int k = 0; k < numfac; k++) {
                            int index2 = jbase_idx + nregressors + ifreefac + k;
                            if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                        }
                    }
                }

                // dbeta dinter
                if (n_interaction_loadings > 0) {
                    for (int ireg = 0; ireg < nregressors; ireg++) {
                        int index1 = ibase_idx + ireg;
                        if (use_free_opt && !hess_idx_free[index1]) continue;  // Skip fixed row
                        double xval_i = getRegValue(ireg, data, iobs_offset);
                        for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                            int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                            if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                        }
                    }
                }

                // dalpha dquad
                if (n_quadratic_loadings > 0) {
                    ifree = 0;
                    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                            int index1 = ibase_idx + nregressors + ifree;
                            if (use_free_opt && !hess_idx_free[index1]) { ifree++; continue; }  // Skip fixed row
                            for (int k = 0; k < numfac; k++) {
                                int index2 = jbase_idx + nregressors + ifreefac + k;
                                if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                            ifree++;
                        }
                    }
                }

                // dalpha dinter
                if (n_interaction_loadings > 0) {
                    ifree = 0;
                    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                            int index1 = ibase_idx + nregressors + ifree;
                            if (use_free_opt && !hess_idx_free[index1]) { ifree++; continue; }  // Skip fixed row
                            for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                                int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                                if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                            ifree++;
                        }
                    }
                }

                // dquad dbeta (for jcat > icat), dquad dalpha, dquad dquad, dquad dinter
                if (n_quadratic_loadings > 0) {
                    for (int ik = 0; ik < numfac; ik++) {
                        int index1 = ibase_idx + nregressors + ifreefac + ik;
                        if (use_free_opt && !hess_idx_free[index1]) continue;  // Skip fixed row
                        double fk_sq_i = fac[ik] * fac[ik];

                        // dquad dbeta (for jcat > icat)
                        if (jcat > icat) {
                            for (int jreg = 0; jreg < nregressors; jreg++) {
                                int index2 = jbase_idx + jreg;
                                if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                            }
                        }

                        // dquad dalpha (for jcat > icat)
                        if (jcat > icat) {
                            int jfree = 0;
                            for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                                if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                                    int index2 = jbase_idx + nregressors + jfree;
                                    if (!use_free_opt || hess_idx_free[index2]) {  // Skip fixed col
                                        hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                                    }
                                    jfree++;
                                }
                            }
                        }

                        // dquad dquad
                        for (int jk = 0; jk < numfac; jk++) {
                            if ((jcat > icat) || (jk >= ik)) {
                                int index2 = jbase_idx + nregressors + ifreefac + jk;
                                if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                            }
                        }

                        // dquad dinter
                        if (n_interaction_loadings > 0) {
                            for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                                int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                                if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                            }
                        }
                    }
                }

                // dinter dbeta, dinter dalpha, dinter dquad, dinter dinter
                if (n_interaction_loadings > 0) {
                    int i_inter_idx = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            int index1 = ibase_idx + nregressors + ifreefac + n_quadratic_loadings + i_inter_idx;
                            if (use_free_opt && !hess_idx_free[index1]) { i_inter_idx++; continue; }  // Skip fixed row
                            double fj_fk_i = fac[ij] * fac[ik];

                            // dinter dbeta (for jcat > icat)
                            if (jcat > icat) {
                                for (int jreg = 0; jreg < nregressors; jreg++) {
                                    int index2 = jbase_idx + jreg;
                                    if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                    hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                }
                            }

                            // dinter dalpha (for jcat > icat)
                            if (jcat > icat) {
                                int jfree = 0;
                                for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                                    if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                                        int index2 = jbase_idx + nregressors + jfree;
                                        if (!use_free_opt || hess_idx_free[index2]) {  // Skip fixed col
                                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                        }
                                        jfree++;
                                    }
                                }
                            }

                            // dinter dquad (for jcat > icat)
                            if (jcat > icat && n_quadratic_loadings > 0) {
                                for (int jk = 0; jk < numfac; jk++) {
                                    int index2 = jbase_idx + nregressors + ifreefac + jk;
                                    if (use_free_opt && !hess_idx_free[index2]) continue;  // Skip fixed col
                                    hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                }
                            }

                            // dinter dinter
                            int j_inter_idx = 0;
                            for (int jj = 0; jj < numfac - 1; jj++) {
                                for (int jk = jj + 1; jk < numfac; jk++) {
                                    if ((jcat > icat) || (j_inter_idx >= i_inter_idx)) {
                                        int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + j_inter_idx;
                                        if (!use_free_opt || hess_idx_free[index2]) {  // Skip fixed col
                                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                        }
                                    }
                                    j_inter_idx++;
                                }
                            }
                            i_inter_idx++;
                        }
                    }
                }
            }
        }
        } // End of slow path else block

        // NOTE: The "HESSIAN CORRECTION FOR FACTOR VARIANCE BLOCK" that was here has been removed.
        // The correction was adding spurious terms because loading_val (computed at lines 1213-1250)
        // already includes the full derivative (linear + quadratic + interaction contributions).
        // The "correction" was double-counting these terms, causing the Hessian to be ~5x too large.
    }

    } // End of rank loop

    // Numerical stability for total likelihood
    if (modEval[0] < 1e-100) modEval[0] = 1e-100;
}

void Model::EvalOprobit(double expres, int outcome_value,
                        const std::vector<double>& fac,
                        const std::vector<double>& param, int firstpar,
                        std::vector<double>& modEval, std::vector<double>& hess,
                        int flag, const std::vector<double>& data, int iobs_offset,
                        const std::vector<int>* model_free_indices)
{
    // Note: model_free_indices is accepted for API consistency but not currently used
    // for optimization in EvalOprobit. The main performance benefit comes from EvalLogit.
    (void)model_free_indices;  // Suppress unused warning

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;

    // Observed category (1, 2, ..., numchoice)
    int obsCat = outcome_value;

    // Threshold parameters come after betas, alphas, and second-order loadings
    int thresh_idx = firstpar + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

    // Build thresholds by accumulating absolute values
    // threshold[0] = lower bound for this category
    // threshold[1] = upper bound for this category
    // Use minimum increment of 0.01 for numerical stability
    const double MIN_THRESH_INCREMENT = 0.01;
    double threshold[2];
    threshold[0] = param[thresh_idx];
    threshold[1] = param[thresh_idx];

    for (int icat = 2; icat <= obsCat; icat++) {
        double incr = std::fabs(param[thresh_idx + icat - 1]);
        if (incr < MIN_THRESH_INCREMENT) incr = MIN_THRESH_INCREMENT;
        if (icat < obsCat) threshold[0] += incr;
        if (icat < numchoice) threshold[1] += incr;
    }

    // Compute CDFs at thresholds
    double CDF[2] = {0.0, 1.0};
    if (obsCat > 1) {
        CDF[0] = normal_cdf(threshold[0] - expres);
    }
    if (obsCat < numchoice) {
        CDF[1] = normal_cdf(threshold[1] - expres);
    }

    // Likelihood: Prob(Y = obsCat) = CDF[upper] - CDF[lower]
    double rawDiffCDF = CDF[1] - CDF[0];

    // Floor the probability to avoid numerical underflow
    // This is critical for multi-factor models where probabilities are multiplied
    const double MIN_PROB = 1.0e-50;
    double diffCDF = rawDiffCDF;
    if (diffCDF < MIN_PROB) diffCDF = MIN_PROB;

    // Return floored probability as the likelihood
    modEval[0] = diffCDF;

    // Fix numerical problems for gradient computation (will divide by CDF later)
    if ((obsCat > 1) && (CDF[0] < MIN_PROB)) CDF[0] = MIN_PROB;
    if (CDF[1] < MIN_PROB) CDF[1] = MIN_PROB;

    if (flag < 2) return;

    // ===== GRADIENT CALCULATION =====

    // Compute PDFs and Z-values at both thresholds
    // Add safeguards for extreme Z-values
    const double MAX_Z = 35.0;  // Beyond this, PDF is essentially 0
    const double MIN_PDF = 1.0e-50;
    double Z[2] = {-9999.0, -9999.0};
    double PDF[2] = {0.0, 0.0};

    if (obsCat > 1) {
        Z[0] = threshold[0] - expres;
        // Bound Z to avoid extreme values in Hessian computation
        if (Z[0] > MAX_Z) Z[0] = MAX_Z;
        if (Z[0] < -MAX_Z) Z[0] = -MAX_Z;
        PDF[0] = normal_pdf(Z[0]);
        if (PDF[0] < MIN_PDF) PDF[0] = MIN_PDF;
    }
    if (obsCat < numchoice) {
        Z[1] = threshold[1] - expres;
        // Bound Z to avoid extreme values in Hessian computation
        if (Z[1] > MAX_Z) Z[1] = MAX_Z;
        if (Z[1] < -MAX_Z) Z[1] = -MAX_Z;
        PDF[1] = normal_pdf(Z[1]);
        if (PDF[1] < MIN_PDF) PDF[1] = MIN_PDF;
    }

    int npar = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings + (numchoice - 1);  // includes thresholds

    // OPTIMIZATION: Only resize if too small, use memset for zeroing
    if (flag == 3) {
        size_t hess_size = static_cast<size_t>(npar * npar);
        if (hess.size() < hess_size) {
            hess.resize(hess_size);
        }
        std::memset(hess.data(), 0, hess_size * sizeof(double));
    }

    // Loop over two terms: lower threshold (iterm=0) and upper threshold (iterm=1)
    for (int iterm = 0; iterm < 2; iterm++) {
        // Skip if they are end categories
        if (((obsCat > 1) && (iterm == 0)) || ((obsCat < numchoice) && (iterm == 1))) {

            std::vector<double> tmpgrad(npar, 0.0);

            // Gradients w.r.t. factor variance parameters and loadings
            int ifree = 0;
            for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                if (facnorm.size() == 0 || facnorm[ifac] <= -9998.0) {
                    // Free loading
                    if (ifac < numfac) {
                        // d/dtheta (factor variance)
                        tmpgrad[ifac] = (-1.0 * param[ifree + firstpar + nregressors]) * PDF[iterm] / diffCDF;
                    }
                    // d/dalpha (factor loading)
                    tmpgrad[numfac + nregressors + ifree] = (-1.0 * fac[ifac]) * PDF[iterm] / diffCDF;
                    ifree++;
                } else {
                    // Fixed loading
                    if (ifac < numfac) {
                        tmpgrad[ifac] = (-1.0 * facnorm[ifac]) * PDF[iterm] / diffCDF;
                    }
                }
            }

            // Gradients w.r.t. regression coefficients
            for (int ireg = 0; ireg < nregressors; ireg++) {
                tmpgrad[ireg + numfac] = (-1.0 * getRegValue(ireg, data, iobs_offset)) * PDF[iterm] / diffCDF;
            }

            // Gradients w.r.t. quadratic factor loadings
            if (n_quadratic_loadings > 0) {
                int quad_grad_start = numfac + nregressors + ifreefac;
                for (int k = 0; k < numfac; k++) {
                    double fk_sq = fac[k] * fac[k];
                    tmpgrad[quad_grad_start + k] = (-1.0 * fk_sq) * PDF[iterm] / diffCDF;

                    // Add quadratic contribution to factor variance gradient
                    double lambda_quad_k = param[firstpar + nregressors + ifreefac + k];
                    tmpgrad[k] += (-1.0 * 2.0 * lambda_quad_k * fac[k]) * PDF[iterm] / diffCDF;
                }
            }

            // Gradients w.r.t. interaction factor loadings
            if (n_interaction_loadings > 0) {
                int inter_grad_start = numfac + nregressors + ifreefac + n_quadratic_loadings;
                int inter_idx = 0;
                for (int j = 0; j < numfac - 1; j++) {
                    for (int k = j + 1; k < numfac; k++) {
                        double fj_fk = fac[j] * fac[k];
                        tmpgrad[inter_grad_start + inter_idx] = (-1.0 * fj_fk) * PDF[iterm] / diffCDF;

                        // Add interaction contribution to factor variance gradients
                        double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx];
                        tmpgrad[j] += (-1.0 * lambda_inter * fac[k]) * PDF[iterm] / diffCDF;
                        tmpgrad[k] += (-1.0 * lambda_inter * fac[j]) * PDF[iterm] / diffCDF;

                        inter_idx++;
                    }
                }
            }

            // Gradients w.r.t. threshold parameters
            int thres_offset = npar - (numchoice - 1);
            int maxthresloop = obsCat;
            if (iterm == 0) maxthresloop--;

            for (int ithres = 0; ithres < maxthresloop; ithres++) {
                tmpgrad[thres_offset + ithres] = PDF[iterm] / diffCDF;
            }

            // Add to total gradient with appropriate sign
            int obsSign = (iterm == 0) ? -1 : 1;
            for (int i = 0; i < npar; i++) {
                modEval[i + 1] += obsSign * tmpgrad[i];
            }
        }
    }

    // ===== HESSIAN CALCULATION =====
    // Using the same factorized approach as the legacy TModel.cc code
    // Initialize tmphess to 1.0 and multiply by factors for each parameter
    if (flag == 3) {
        for (int iterm = 0; iterm < 2; iterm++) {
            // Skip if end categories
            if (((obsCat > 1) && (iterm == 0)) || ((obsCat < numchoice) && (iterm == 1))) {

                std::vector<double> tmphess(npar * npar, 1.0);
                int obsSign = (iterm == 0) ? -1 : 1;

                // Factor-specific Hessian terms (theta and alpha)
                int ifree = 0;
                for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                    // No normalizations (all free) or check if this one is free
                    if (facnorm.size() == 0 || facnorm[ifac] <= -9998.0) {
                        if (ifac < numfac) {
                            // lambda^L(theta) - lambda^Prob(theta) (row)
                            for (int j = ifac; j < npar; j++) {
                                tmphess[ifac*npar + j] *= -Z[iterm] * (-1.0 * param[ifree + firstpar + nregressors]) - modEval[1 + ifac];
                            }
                            // dZ/dtheta_i (col)
                            for (int j = 0; j <= ifac; j++) {
                                tmphess[j*npar + ifac] *= -1.0 * param[ifree + firstpar + nregressors];
                            }
                        }

                        // alpha_i index
                        int index = numfac + nregressors + ifree;
                        // lambda^L(alpha) - lambda^Prob(alpha) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * fac[ifac]) - modEval[1 + numfac + nregressors + ifree];
                        }
                        // dZ/dalpha (col)
                        for (int j = 0; j <= index; j++) {
                            tmphess[j*npar + index] *= -1.0 * fac[ifac];
                        }

                        ifree++;
                    } else {
                        // Fixed loading (facnorm[ifac] > -9998)
                        if (ifac < numfac) {
                            // lambda^L(theta) - lambda^Prob(theta) (row)
                            for (int j = ifac; j < npar; j++) {
                                tmphess[ifac*npar + j] *= -Z[iterm] * (-1.0 * facnorm[ifac]) - modEval[1 + ifac];
                            }
                            // dZ/dtheta_i (col)
                            for (int j = 0; j <= ifac; j++) {
                                tmphess[j*npar + ifac] *= -1.0 * facnorm[ifac];
                            }
                        }
                    }
                }

                // Hessian for regression coefficients (X's)
                for (int ireg = 0; ireg < nregressors; ireg++) {
                    int index = numfac + ireg;
                    // lambda^L(X) - lambda^Prob(X) (row)
                    for (int j = index; j < npar; j++) {
                        tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * getRegValue(ireg, data, iobs_offset)) - modEval[1 + ireg + numfac];
                    }
                    // dZ/dX_i (col)
                    for (int j = 0; j <= index; j++) {
                        tmphess[j*npar + index] *= -1.0 * getRegValue(ireg, data, iobs_offset);
                    }
                }

                // Hessian for quadratic factor loadings
                if (n_quadratic_loadings > 0) {
                    int quad_offset = numfac + nregressors + ifreefac;
                    for (int k = 0; k < numfac; k++) {
                        double fk_sq = fac[k] * fac[k];
                        int index = quad_offset + k;
                        // lambda^L(quad) - lambda^Prob(quad) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * fk_sq) - modEval[1 + quad_offset + k];
                        }
                        // dZ/dquad_k (col)
                        for (int j = 0; j <= index; j++) {
                            tmphess[j*npar + index] *= -1.0 * fk_sq;
                        }
                    }
                }

                // Hessian for interaction factor loadings
                if (n_interaction_loadings > 0) {
                    int inter_offset = numfac + nregressors + ifreefac + n_quadratic_loadings;
                    int inter_idx = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            double fj_fk = fac[ij] * fac[ik];
                            int index = inter_offset + inter_idx;
                            // lambda^L(inter) - lambda^Prob(inter) (row)
                            for (int j = index; j < npar; j++) {
                                tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * fj_fk) - modEval[1 + inter_offset + inter_idx];
                            }
                            // dZ/dinter_jk (col)
                            for (int j = 0; j <= index; j++) {
                                tmphess[j*npar + index] *= -1.0 * fj_fk;
                            }
                            inter_idx++;
                        }
                    }
                }

                // Hessian for thresholds
                int thres_offset = npar - (numchoice - 1);
                int maxthresloop = obsCat;
                if (iterm == 0) maxthresloop--;

                for (int ithres = 0; ithres < numchoice - 1; ithres++) {
                    int index = thres_offset + ithres;
                    if (ithres < maxthresloop) {
                        // lambda^L(chi) - lambda^Prob(chi) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -Z[iterm] - modEval[1 + thres_offset + ithres];
                        }
                    } else {
                        // lambda^L(chi) - lambda^Prob(chi) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -modEval[1 + thres_offset + ithres];
                        }
                        // dZ/dchi_i (col) - set to zero
                        for (int j = 0; j <= index; j++) {
                            tmphess[j*npar + index] = 0.0;
                        }
                    }
                }

                // Add cross-derivative term dZ/dtheta dalpha = -1
                ifree = 0;
                for (int i = 0; i < numfac; i++) {
                    if (facnorm.size() == 0 || facnorm[i] <= -9998.0) {
                        int index = numfac + nregressors + ifree;
                        tmphess[i*npar + index] += -1.0;
                        ifree++;
                    }
                }

                // Add cross-derivative term dZ/dtheta_k dlambda_quad_k = -2*f_k
                if (n_quadratic_loadings > 0) {
                    int quad_offset = numfac + nregressors + ifreefac;
                    for (int k = 0; k < numfac; k++) {
                        int index = quad_offset + k;
                        tmphess[k*npar + index] += -2.0 * fac[k];
                    }
                }

                // Add cross-derivative terms dZ/dtheta dlambda_inter
                if (n_interaction_loadings > 0) {
                    int inter_offset = numfac + nregressors + ifreefac + n_quadratic_loadings;
                    int inter_idx = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            int index = inter_offset + inter_idx;
                            tmphess[ij*npar + index] += -fac[ik];  // d(theta_j)/d(lambda_inter_jk)
                            tmphess[ik*npar + index] += -fac[ij];  // d(theta_k)/d(lambda_inter_jk)
                            inter_idx++;
                        }
                    }
                }

                // Add second derivative term d²Z/df_k² = -2*λ_q_k for factor-factor diagonal
                // Note: Z = threshold - xb, so d²Z/df² = -d²xb/df² = -2*λ_q_k (negative sign!)
                // This comes from: d²L/df² = hess_factor * (dZ/df)² + (PDF/CDF) * obsSign * d²Z/df²
                // The first term is in the HESSIAN CORRECTION block; this is the second term.
                if (n_quadratic_loadings > 0) {
                    for (int k = 0; k < numfac; k++) {
                        double lambda_quad_k = param[firstpar + nregressors + ifreefac + k];
                        tmphess[k*npar + k] += -2.0 * lambda_quad_k;  // Negative because d²Z/df² = -d²xb/df²
                    }
                }

                // Add second derivative term d²Z/df_j df_k = -λ_inter for factor-factor cross
                // Note: Z = threshold - xb, so d²Z/df_j df_k = -d²xb/df_j df_k = -λ_inter
                if (n_interaction_loadings > 0) {
                    int inter_idx2 = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx2];
                            tmphess[ij*npar + ik] += -lambda_inter;  // Negative because d²Z/df = -d²xb/df
                            inter_idx2++;
                        }
                    }
                }

                // ===== HESSIAN CORRECTION FOR FACTOR VARIANCE BLOCK (OPROBIT) =====
                // The factor variance row/column multiplications use λ_k, but should use (λ_k + add_k)
                // Row factor: Z * λ_k - grad_k  ->  Z * (λ_k + add_k) - grad_k
                // Col factor: -λ_k  ->  -(λ_k + add_k)
                // This adds the missing cross-product terms involving additional derivative contributions
                if (n_quadratic_loadings > 0 || n_interaction_loadings > 0) {
                    // Compute additional derivative terms for each factor
                    std::vector<double> additional(numfac, 0.0);
                    std::vector<double> linear_loading(numfac, 0.0);

                    // Get linear loading values
                    int ifree_tmp = 0;
                    for (int k = 0; k < numfac; k++) {
                        if (facnorm.size() == 0 || facnorm[k] <= -9998) {
                            linear_loading[k] = param[firstpar + nregressors + ifree_tmp];
                            ifree_tmp++;
                        } else {
                            linear_loading[k] = facnorm[k];
                        }
                    }

                    // Add quadratic contributions to additional
                    if (n_quadratic_loadings > 0) {
                        for (int k = 0; k < numfac; k++) {
                            double lambda_quad_k = param[firstpar + nregressors + ifreefac + k];
                            additional[k] += 2.0 * lambda_quad_k * fac[k];
                        }
                    }

                    // Add interaction contributions to additional
                    if (n_interaction_loadings > 0) {
                        int inter_idx2 = 0;
                        for (int j = 0; j < numfac - 1; j++) {
                            for (int k = j + 1; k < numfac; k++) {
                                double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx2];
                                additional[j] += lambda_inter * fac[k];
                                additional[k] += lambda_inter * fac[j];
                                inter_idx2++;
                            }
                        }
                    }

                    // Correction for factor variance block
                    // Current: (Z*λ_j - grad_j) * (-λ_k)
                    // Target: (Z*(λ_j+add_j) - grad_j) * (-(λ_k+add_k))
                    // Note: grad already includes add contributions, so row factor error is Z*add_j
                    // Correction: Z*add_j*(-λ_k) + (Z*λ_j - grad_j)*(-add_k) + Z*add_j*(-add_k)
                    //           = -Z*add_j*λ_k - (Z*λ_j - grad_j)*add_k - Z*add_j*add_k
                    for (int j = 0; j < numfac; j++) {
                        for (int k = j; k < numfac; k++) {
                            double row_j = -Z[iterm] * (-linear_loading[j]) - modEval[1 + j];
                            double correction;
                            if (j == k) {
                                // Diagonal: row and col are same factor
                                correction = -Z[iterm] * additional[k] * linear_loading[k]
                                           - row_j * additional[k]
                                           - Z[iterm] * additional[k] * additional[k];
                            } else {
                                correction = -Z[iterm] * additional[j] * linear_loading[k]
                                           - row_j * additional[k]
                                           - Z[iterm] * additional[j] * additional[k];
                            }
                            tmphess[j * npar + k] += correction;
                        }
                    }

                    // Correction for cross-derivatives: θ_k with other parameters
                    // Row k should have factor: Z*λ_k - grad_k -> Z*(λ_k + add_k) - grad_k
                    // So we need to add: Z * add_k * col_j for each param j
                    for (int k = 0; k < numfac; k++) {
                        if (std::abs(additional[k]) < 1e-15) continue;

                        // Cross with regression coefficients
                        for (int ireg = 0; ireg < nregressors; ireg++) {
                            double col_reg = -getRegValue(ireg, data, iobs_offset);
                            int reg_idx = numfac + ireg;
                            // Row correction: Z*add_k * col_reg (col_reg already has -1 factor)
                            tmphess[k * npar + reg_idx] += Z[iterm] * additional[k] * col_reg;
                        }

                        // Cross with free linear loadings
                        int ifree3 = 0;
                        for (int m = 0; m < numfac; m++) {
                            if (facnorm.size() == 0 || facnorm[m] <= -9998) {
                                int load_idx = numfac + nregressors + ifree3;
                                double col_load = -fac[m];
                                tmphess[k * npar + load_idx] += Z[iterm] * additional[k] * col_load;
                                ifree3++;
                            }
                        }

                        // Cross with quadratic loadings
                        if (n_quadratic_loadings > 0) {
                            for (int m = 0; m < numfac; m++) {
                                int quad_idx = numfac + nregressors + ifreefac + m;
                                double col_quad = -fac[m] * fac[m];
                                tmphess[k * npar + quad_idx] += Z[iterm] * additional[k] * col_quad;
                            }
                        }

                        // Cross with interaction loadings
                        if (n_interaction_loadings > 0) {
                            int inter_idx3 = 0;
                            for (int j2 = 0; j2 < numfac - 1; j2++) {
                                for (int k2 = j2 + 1; k2 < numfac; k2++) {
                                    int idx = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx3;
                                    double col_inter = -fac[j2] * fac[k2];
                                    tmphess[k * npar + idx] += Z[iterm] * additional[k] * col_inter;
                                    inter_idx3++;
                                }
                            }
                        }

                        // Cross with thresholds (col factor is 1 for active thresholds)
                        int thres_start = npar - (numchoice - 1);
                        int maxthres = obsCat;
                        if (iterm == 0) maxthres--;
                        for (int ithres = 0; ithres < maxthres; ithres++) {
                            int idx = thres_start + ithres;
                            // Threshold has col factor of 1 (not -1)
                            tmphess[k * npar + idx] += Z[iterm] * additional[k] * 1.0;
                        }
                    }
                }

                // Add tmp hessian to totals, scaled by obsSign * PDF / diffCDF
                for (int i = 0; i < npar; i++) {
                    for (int j = i; j < npar; j++) {
                        hess[i*npar + j] += obsSign * tmphess[i*npar + j] * PDF[iterm] / diffCDF;
                    }
                }
            }
        }
    }
}
