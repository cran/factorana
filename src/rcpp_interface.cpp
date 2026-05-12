#include <Rcpp.h>
#include <RcppEigen.h>
#include <iomanip>
#include "FactorModel.h"
#include "Model.h"
#include "gauss_hermite.h"

using namespace Rcpp;

// Expose FactorModel class to R using Rcpp modules
RCPP_MODULE(factorana_module) {
    // Use function pointers to disambiguate overloaded methods
    void (FactorModel::*setData1)(const std::vector<double>&) = &FactorModel::SetData;
    void (FactorModel::*setData2)(const Eigen::MatrixXd&) = &FactorModel::SetData;
    void (FactorModel::*setConstraints1)(const std::vector<bool>&) = &FactorModel::SetParameterConstraints;
    void (FactorModel::*setConstraints2)(const std::vector<bool>&, const std::vector<double>&) = &FactorModel::SetParameterConstraints;

    class_<FactorModel>("FactorModel")
        .constructor<int, int, int, int, int, bool, int>("Create a FactorModel")
        .method("SetDataVector", setData1)
        .method("SetDataMatrix", setData2)
        .method("SetQuadrature", &FactorModel::SetQuadrature)
        .method("SetParameterConstraints", setConstraints1)
        .method("SetParameterConstraintsWithValues", setConstraints2)
        .method("CalcLogLikelihood", &FactorModel::CalcLogLikelihood)
        .method("GetNObs", &FactorModel::GetNObs)
        .method("GetNParam", &FactorModel::GetNParam)
        .method("GetNParamFree", &FactorModel::GetNParamFree)
        ;
}

//' Compute Gauss-Hermite quadrature nodes and weights
//'
//' @param n Number of quadrature points
//' @return List with nodes and weights
//' @export
// [[Rcpp::export]]
List gauss_hermite_quadrature(int n) {
    std::vector<double> nodes, weights;
    calcgausshermitequadrature(n, nodes, weights);

    return List::create(
        Named("nodes") = nodes,
        Named("weights") = weights
    );
}

// Helper: apply user-fixed factor-distribution parameters from
// factor$fixed_params (set in R via fix_factor_param()). Builds a
// constructor-true name->index map for the factor-level parameter block
// (factor_var, [factor_corr], mix_means, mix_logweight, se_*,
//  typeprob/type_loading, factor_mean, se_cov), looks up each fixed name,
// marks param_fixed_vec[idx] = true, and overrides init_params[idx] with
// the user-supplied value. The same layout is also used by the un-fix
// loop at the top of initialize_factor_model_cpp and by the
// equality_constraints map; keep all three in sync if the layout changes.
static inline void apply_fix_factor_param(
    Rcpp::List factor_model,
    FactorStructure fac_struct, int n_fac, int n_types, int n_mixtures,
    bool any_uses_types, int param_offset,
    std::vector<bool>& param_fixed_vec,
    Rcpp::Nullable<Rcpp::NumericVector> init_params)
{
    if (!factor_model.containsElementNamed("fixed_params") ||
        Rf_isNull(factor_model["fixed_params"])) return;

    Rcpp::NumericVector fp = factor_model["fixed_params"];
    if (fp.size() == 0) return;

    Rcpp::CharacterVector fp_names = fp.names();
    if (fp_names.size() != fp.size()) return;  // unnamed: nothing to do

    bool is_se = (fac_struct == FactorStructure::SE_LINEAR ||
                  fac_struct == FactorStructure::SE_QUADRATIC);
    int n_var_factors = is_se ? (n_fac - 1) : n_fac;

    // Build the same name -> index map used by the un-fix loop.
    std::map<std::string, int> name_to_idx;
    int idx = 0;

    // Block 1: factor variances (per mixture)
    for (int m = 0; m < n_mixtures; m++) {
        for (int k = 0; k < n_var_factors; k++) {
            std::string nm = (n_mixtures == 1)
                ? "factor_var_" + std::to_string(k + 1)
                : "mix" + std::to_string(m + 1) + "_factor_var_" + std::to_string(k + 1);
            name_to_idx[nm] = idx++;
        }
    }

    // Correlation parameter
    if (fac_struct == FactorStructure::CORRELATION && n_fac == 2) {
        name_to_idx["factor_corr_1_2"] = idx++;
    }

    // Block 2: mixture means + log-weights
    for (int m = 0; m < n_mixtures - 1; m++) {
        for (int k = 0; k < n_var_factors; k++) {
            name_to_idx["mix" + std::to_string(m+1) + "_factor_mean_" + std::to_string(k+1)] = idx++;
        }
    }
    for (int m = 0; m < n_mixtures - 1; m++) {
        name_to_idx["mix" + std::to_string(m+1) + "_logweight"] = idx++;
    }

    // Block 3: SE parameters (SE structures only)
    if (is_se) {
        name_to_idx["se_intercept"] = idx++;
        for (int k = 0; k < n_var_factors; k++) {
            name_to_idx["se_linear_" + std::to_string(k+1)] = idx++;
        }
        if (fac_struct == FactorStructure::SE_QUADRATIC) {
            for (int k = 0; k < n_var_factors; k++) {
                name_to_idx["se_quadratic_" + std::to_string(k+1)] = idx++;
            }
        }
        if (n_types > 1) {
            for (int t = 2; t <= n_types; t++) {
                name_to_idx["se_intercept_type_" + std::to_string(t)] = idx++;
            }
        }
        name_to_idx["se_residual_var"] = idx++;
    }

    // Block 4: typeprob_*_intercept + type_*_loading_*
    if (n_types > 1 && any_uses_types) {
        for (int t = 2; t <= n_types; t++) {
            name_to_idx["typeprob_" + std::to_string(t) + "_intercept"] = idx++;
        }
        for (int t = 2; t <= n_types; t++) {
            for (int k = 1; k <= n_fac; k++) {
                name_to_idx["type_" + std::to_string(t) + "_loading_" + std::to_string(k)] = idx++;
            }
        }
    }

    // Block 5: factor_mean_<k>_<cov>
    if (factor_model.containsElementNamed("factor_covariates") &&
        !Rf_isNull(factor_model["factor_covariates"])) {
        Rcpp::CharacterVector fcov = factor_model["factor_covariates"];
        if (fcov.size() > 0) {
            int n_fac_with_mean = is_se ? (n_fac - 1) : n_fac;
            for (int k = 1; k <= n_fac_with_mean; k++) {
                for (int j = 0; j < fcov.size(); j++) {
                    std::string cn = Rcpp::as<std::string>(fcov[j]);
                    name_to_idx["factor_mean_" + std::to_string(k) + "_" + cn] = idx++;
                }
            }
        }
    }

    // Block 6: se_cov_<cov>
    if (factor_model.containsElementNamed("se_covariates") &&
        !Rf_isNull(factor_model["se_covariates"])) {
        Rcpp::CharacterVector secov = factor_model["se_covariates"];
        for (int j = 0; j < secov.size(); j++) {
            std::string cn = Rcpp::as<std::string>(secov[j]);
            name_to_idx["se_cov_" + cn] = idx++;
        }
    }

    // Apply each user fix.
    for (int j = 0; j < fp.size(); j++) {
        std::string nm = Rcpp::as<std::string>(fp_names[j]);
        auto it = name_to_idx.find(nm);
        if (it == name_to_idx.end()) continue;          // not a factor-dist name
        int pos = it->second;
        if (pos < 0 || pos >= (int) param_fixed_vec.size()) continue;
        param_fixed_vec[pos] = true;
        // Override init_params at this position so SetParameterConstraintsWithValues
        // uses the user-supplied value as the FIXED value the C++ side reads.
        if (!init_params.isNull()) {
            Rcpp::NumericVector ip(init_params);
            if (pos < ip.size()) {
                ip[pos] = fp[j];
            }
        }
    }
}


//' Initialize a FactorModel C++ object from R model system
//'
//' @param model_system R model_system object
//' @param data Data frame or matrix with all variables
//' @param n_quad Number of quadrature points
//' @param init_params Optional initial parameter vector (used to set fixed parameter values)
//' @return External pointer to FactorModel object
//' @export
// [[Rcpp::export]]
SEXP initialize_factor_model_cpp(List model_system, SEXP data, int n_quad = 8,
                                  Nullable<NumericVector> init_params = R_NilValue) {

    // Extract factor model information
    List factor_model = model_system["factor"];
    int n_fac = factor_model["n_factors"];
    int n_types = factor_model["n_types"];
    int n_mixtures = factor_model["n_mixtures"];
    bool correlation = factor_model["correlation"];

    // Parse factor_structure
    FactorStructure fac_struct = FactorStructure::INDEPENDENT;
    if (factor_model.containsElementNamed("factor_structure")) {
        std::string fac_struct_str = as<std::string>(factor_model["factor_structure"]);
        if (fac_struct_str == "correlation") {
            fac_struct = FactorStructure::CORRELATION;
        } else if (fac_struct_str == "SE_linear") {
            fac_struct = FactorStructure::SE_LINEAR;
        } else if (fac_struct_str == "SE_quadratic") {
            fac_struct = FactorStructure::SE_QUADRATIC;
        } else if (fac_struct_str == "independent") {
            fac_struct = FactorStructure::INDEPENDENT;
        }
    } else if (correlation) {
        // Backward compatibility: correlation = TRUE maps to CORRELATION
        fac_struct = FactorStructure::CORRELATION;
    }

    // Convert data to matrix (handle both data.frame and matrix)
    NumericMatrix data_mat;
    if (Rf_isMatrix(data)) {
        data_mat = as<NumericMatrix>(data);
    } else {
        // Assume it's a data frame
        DataFrame df = as<DataFrame>(data);
        data_mat = internal::convert_using_rfunction(df, "as.matrix");
    }

    int n_obs = data_mat.nrow();
    int n_var = data_mat.ncol();

    // Get column names for variable indexing
    CharacterVector col_names = colnames(data_mat);
    if (col_names.size() == 0) {
        Rcpp::stop("Data must have column names for variable indexing");
    }

    // Create FactorModel object using new constructor with factor_structure
    Rcpp::XPtr<FactorModel> fm(
        new FactorModel(n_obs, n_var, n_fac, n_types, n_mixtures, fac_struct, n_quad),
        true
    );

    // Set data
    std::vector<double> data_vec(n_obs * n_var);
    for (int i = 0; i < n_obs; i++) {
        for (int j = 0; j < n_var; j++) {
            data_vec[i * n_var + j] = data_mat(i, j);
        }
    }
    fm->SetData(data_vec);

    // Compute and set quadrature
    std::vector<double> nodes, weights;

    // Special case: if n_fac == 0, use single integration point at 0 with weight 1
    if (n_fac == 0) {
        nodes.push_back(0.0);
        weights.push_back(1.0);
    } else {
        calcgausshermitequadrature(n_quad, nodes, weights);

        // Scale for standard normal: x_scaled = sqrt(2) * x, w_scaled = w / sqrt(pi)
        double sqrt2 = std::sqrt(2.0);
        double sqrt_pi = std::sqrt(M_PI);
        for (int i = 0; i < n_quad; i++) {
            nodes[i] *= sqrt2;
            weights[i] /= sqrt_pi;
        }
    }

    fm->SetQuadrature(nodes, weights);

    // Set up factor mean covariates if specified
    // MUST be done BEFORE AddModel() so factor_mean params come before model params
    // This ensures parameter ordering matches R: factor_var, factor_mean, model_params
    if (factor_model.containsElementNamed("factor_covariates") &&
        !Rf_isNull(factor_model["factor_covariates"])) {
        CharacterVector factor_cov_names = factor_model["factor_covariates"];
        int n_cov = factor_cov_names.size();
        if (n_cov > 0) {
            // Find column indices for factor covariates
            std::vector<int> factor_cov_indices(n_cov);
            for (int j = 0; j < n_cov; j++) {
                std::string cov_name = as<std::string>(factor_cov_names[j]);
                int idx = -1;
                for (int k = 0; k < col_names.size(); k++) {
                    if (std::string(col_names[k]) == cov_name) {
                        idx = k;
                        break;
                    }
                }
                if (idx == -1) {
                    Rcpp::stop("Factor mean covariate '" + cov_name + "' not found in data");
                }
                factor_cov_indices[j] = idx;
            }

            // Extract covariate values and compute means
            std::vector<double> cov_means(n_cov, 0.0);
            for (int i = 0; i < n_obs; i++) {
                for (int j = 0; j < n_cov; j++) {
                    cov_means[j] += data_mat(i, factor_cov_indices[j]);
                }
            }
            for (int j = 0; j < n_cov; j++) {
                cov_means[j] /= static_cast<double>(n_obs);
            }

            // Extract, demean, and check variance
            std::vector<std::vector<double>> covariate_data(n_obs, std::vector<double>(n_cov));
            std::vector<double> cov_var(n_cov, 0.0);
            for (int i = 0; i < n_obs; i++) {
                for (int j = 0; j < n_cov; j++) {
                    double demeaned = data_mat(i, factor_cov_indices[j]) - cov_means[j];
                    covariate_data[i][j] = demeaned;
                    cov_var[j] += demeaned * demeaned;
                }
            }

            // Check for near-zero variance (constant covariates)
            for (int j = 0; j < n_cov; j++) {
                cov_var[j] /= static_cast<double>(n_obs);
                if (cov_var[j] < 1e-10) {
                    std::string cov_name = as<std::string>(factor_cov_names[j]);
                    Rcpp::stop("Factor mean covariate '" + cov_name + "' has near-zero variance " +
                              "(variance = " + std::to_string(cov_var[j]) + "). " +
                              "This covariate is constant or near-constant and its coefficient " +
                              "is not identified. Remove it from factor_covariates.");
                }
            }

            fm->SetFactorMeanCovariates(covariate_data);
        }
    }

    // Set up SE covariates if specified (for SE_linear/SE_quadratic)
    // These directly affect the outcome factor: f_k = ... + beta * X + epsilon
    if (factor_model.containsElementNamed("se_covariates") &&
        !Rf_isNull(factor_model["se_covariates"])) {
        CharacterVector se_cov_names = factor_model["se_covariates"];
        int n_cov = se_cov_names.size();
        if (n_cov > 0) {
            // Find column indices for SE covariates
            std::vector<int> se_cov_indices(n_cov);
            for (int j = 0; j < n_cov; j++) {
                std::string cov_name = as<std::string>(se_cov_names[j]);
                int idx = -1;
                for (int k = 0; k < col_names.size(); k++) {
                    if (std::string(col_names[k]) == cov_name) {
                        idx = k;
                        break;
                    }
                }
                if (idx == -1) {
                    Rcpp::stop("SE covariate '" + cov_name + "' not found in data");
                }
                se_cov_indices[j] = idx;
            }

            // Extract covariate values and compute means
            std::vector<double> cov_means(n_cov, 0.0);
            for (int i = 0; i < n_obs; i++) {
                for (int j = 0; j < n_cov; j++) {
                    cov_means[j] += data_mat(i, se_cov_indices[j]);
                }
            }
            for (int j = 0; j < n_cov; j++) {
                cov_means[j] /= static_cast<double>(n_obs);
            }

            // Extract, demean, and check variance
            std::vector<std::vector<double>> se_covariate_data(n_obs, std::vector<double>(n_cov));
            std::vector<double> cov_var(n_cov, 0.0);
            for (int i = 0; i < n_obs; i++) {
                for (int j = 0; j < n_cov; j++) {
                    double demeaned = data_mat(i, se_cov_indices[j]) - cov_means[j];
                    se_covariate_data[i][j] = demeaned;
                    cov_var[j] += demeaned * demeaned;
                }
            }

            // Check for near-zero variance (constant covariates)
            for (int j = 0; j < n_cov; j++) {
                cov_var[j] /= static_cast<double>(n_obs);
                if (cov_var[j] < 1e-10) {
                    std::string cov_name = as<std::string>(se_cov_names[j]);
                    Rcpp::stop("SE covariate '" + cov_name + "' has near-zero variance " +
                              "(variance = " + std::to_string(cov_var[j]) + "). " +
                              "This covariate is constant or near-constant and its coefficient " +
                              "is not identified. Remove it from se_covariates.");
                }
            }

            fm->SetSECovariates(se_covariate_data);
        }
    }

    // Track which factors are identified via fixed non-zero loadings
    // A factor is identified if at least one component has loading_normalization != NA and != 0
    std::vector<bool> factor_identified(n_fac, false);

    // Add model components
    List components = model_system["components"];

    // Pre-scan: check if any component uses types (needed for type loadings)
    bool any_uses_types = false;
    for (int i = 0; i < components.size(); i++) {
        List comp = components[i];
        if (comp.containsElementNamed("use_types") && !Rf_isNull(comp["use_types"])) {
            if (as<bool>(comp["use_types"])) {
                any_uses_types = true;
                break;
            }
        }
    }
    // SE_linear / SE_quadratic with n_types > 1 implies types at the structural level
    // (se_intercept_type_{t}), so the type probability model parameters are needed
    // even when no measurement component has use_types = TRUE.
    if (!any_uses_types && n_types > 1 &&
        (fac_struct == FactorStructure::SE_LINEAR ||
         fac_struct == FactorStructure::SE_QUADRATIC)) {
        any_uses_types = true;
    }
    for (int i = 0; i < components.size(); i++) {
        List comp = components[i];

        // Extract model component information
        std::string model_type_str = as<std::string>(comp["model_type"]);

        // Handle outcome - can be single string or vector (for exploded logit)
        std::vector<std::string> outcome_names;
        SEXP outcome_sexp = comp["outcome"];
        if (TYPEOF(outcome_sexp) == STRSXP) {
            CharacterVector outcome_vec = as<CharacterVector>(outcome_sexp);
            for (int j = 0; j < outcome_vec.size(); j++) {
                outcome_names.push_back(std::string(outcome_vec[j]));
            }
        } else {
            Rcpp::stop("outcome must be a character vector");
        }

        // Handle covariates = NULL case
        std::vector<std::string> covariate_names;
        if (!Rf_isNull(comp["covariates"])) {
            covariate_names = as<std::vector<std::string>>(comp["covariates"]);
        }

        // Check if this is a dynamic model first
        bool is_dynamic = comp.containsElementNamed("is_dynamic") &&
                         as<bool>(comp["is_dynamic"]);

        // Find outcome variable indices in data
        // For exploded logit, we need the index of the first outcome
        // The Model will read consecutive columns: outcome_idx, outcome_idx+1, ..., outcome_idx+nrank-1
        int outcome_idx = -1;

        // Verify all outcome columns exist and are consecutive
        std::vector<int> outcome_indices;
        for (const auto& outcome_name : outcome_names) {
            int idx = -1;
            for (int j = 0; j < col_names.size(); j++) {
                if (std::string(col_names[j]) == outcome_name) {
                    idx = j;
                    break;
                }
            }
            if (idx == -1) {
                if (is_dynamic && outcome_names.size() == 1) {
                    idx = -2;  // Special marker: outcome is always zero
                } else {
                    Rcpp::stop("Outcome variable '" + outcome_name + "' not found in data");
                }
            }
            outcome_indices.push_back(idx);
        }

        // For exploded logit, outcomes don't need to be consecutive in data
        // We'll store first outcome index. The Model will use outcome_indices directly.
        outcome_idx = outcome_indices[0];

        // Find covariate indices
        std::vector<int> regressor_idx;
        for (const auto& cov_name : covariate_names) {
            int idx = -1;
            for (int j = 0; j < col_names.size(); j++) {
                if (std::string(col_names[j]) == cov_name) {
                    idx = j;
                    break;
                }
            }
            // For dynamic models, 'intercept' covariate is always 1
            // Use special marker -3 if not found in data
            if (idx == -1) {
                if (is_dynamic && cov_name == "intercept") {
                    idx = -3;  // Special marker: intercept is always 1
                } else {
                    Rcpp::stop("Covariate '" + cov_name + "' not found in data");
                }
            }
            regressor_idx.push_back(idx);
        }

        // Find missing/evaluation indicator if present
        int missing_idx = -1;
        if (comp.containsElementNamed("evaluation_indicator") &&
            !Rf_isNull(comp["evaluation_indicator"])) {
            std::string eval_name = as<std::string>(comp["evaluation_indicator"]);
            for (int j = 0; j < col_names.size(); j++) {
                if (std::string(col_names[j]) == eval_name) {
                    missing_idx = j;
                    break;
                }
            }
        }

        // Determine model type
        ModelType mtype;
        if (model_type_str == "linear") mtype = ModelType::LINEAR;
        else if (model_type_str == "probit") mtype = ModelType::PROBIT;
        else if (model_type_str == "logit") mtype = ModelType::LOGIT;
        else if (model_type_str == "oprobit") mtype = ModelType::OPROBIT;
        else {
            Rcpp::stop("Unknown model type: " + model_type_str);
        }

        // Extract factor normalizations from COMPONENT (not factor model)
        // This allows each component to have its own loading constraints
        std::vector<double> facnorm;
        if (comp.containsElementNamed("loading_normalization")) {
            SEXP norm_sexp = comp["loading_normalization"];
            if (!Rf_isNull(norm_sexp)) {
                NumericVector norm_vec = as<NumericVector>(norm_sexp);
                for (int j = 0; j < norm_vec.size(); j++) {
                    if (NumericVector::is_na(norm_vec[j])) {
                        facnorm.push_back(-9999.0);  // Free parameter
                    } else {
                        facnorm.push_back(norm_vec[j]);  // Fixed value
                        // Track factor identification: if loading is fixed AND non-zero,
                        // the factor variance is identified (can be estimated)
                        if (std::abs(norm_vec[j]) > 1e-6 && j < n_fac) {
                            factor_identified[j] = true;
                        }
                    }
                }
            }
        }

        int n_choice = comp.containsElementNamed("num_choices") ?
                      int(comp["num_choices"]) : 2;

        // Number of ranks for exploded logit (default 1 = standard logit)
        int n_rank = comp.containsElementNamed("nrank") && !Rf_isNull(comp["nrank"]) ?
                    int(comp["nrank"]) : 1;

        // Check if all parameters are fixed (for multi-stage estimation)
        bool all_params_fixed = false;
        if (comp.containsElementNamed("all_params_fixed")) {
            all_params_fixed = as<bool>(comp["all_params_fixed"]);
        }

        // Extract factor_spec (linear, quadratic, interactions, full)
        FactorSpec fspec = FactorSpec::LINEAR;
        if (comp.containsElementNamed("factor_spec")) {
            std::string fs = as<std::string>(comp["factor_spec"]);
            if (fs == "quadratic") fspec = FactorSpec::QUADRATIC;
            else if (fs == "interactions") fspec = FactorSpec::INTERACTIONS;
            else if (fs == "full") fspec = FactorSpec::FULL;
        }

        // Extract dynamic model outcome_factor index (is_dynamic already set above)
        int outcome_factor_idx = -1;
        if (is_dynamic) {
            // outcome_factor is 1-based in R, convert to 0-based for C++
            outcome_factor_idx = as<int>(comp["outcome_factor"]) - 1;
        }

        // Extract exclude_chosen for exploded nested logit (default true = standard exploded logit)
        bool exclude_chosen = true;
        if (comp.containsElementNamed("exclude_chosen") && !Rf_isNull(comp["exclude_chosen"])) {
            exclude_chosen = as<bool>(comp["exclude_chosen"]);
        }

        // Extract rankshare_var index for exploded nested logit
        int ranksharevar_idx = -1;
        if (comp.containsElementNamed("rankshare_var") && !Rf_isNull(comp["rankshare_var"])) {
            std::string rankshare_name = as<std::string>(comp["rankshare_var"]);
            for (int j = 0; j < col_names.size(); j++) {
                if (std::string(col_names[j]) == rankshare_name) {
                    ranksharevar_idx = j;
                    break;
                }
            }
            if (ranksharevar_idx == -1) {
                Rcpp::warning("rankshare_var '" + rankshare_name + "' not found in data, ignoring.");
            }
        }

        // Extract use_types flag (whether this component uses type-specific intercepts)
        bool comp_use_types = false;
        if (comp.containsElementNamed("use_types") && !Rf_isNull(comp["use_types"])) {
            comp_use_types = as<bool>(comp["use_types"]);
        }

        // Create Model object
        std::shared_ptr<Model> model = std::make_shared<Model>(
            mtype, outcome_idx, missing_idx, regressor_idx,
            n_fac, n_types, facnorm, n_choice, n_rank, all_params_fixed, fspec,
            is_dynamic, outcome_factor_idx, outcome_indices,
            exclude_chosen, ranksharevar_idx, comp_use_types
        );

        // Calculate number of parameters for this model
        int n_free_loadings = 0;
        for (const auto& norm : facnorm) {
            if (norm <= -9998.0) n_free_loadings++;
        }
        if (facnorm.empty()) n_free_loadings = n_fac;

        // Calculate second-order loading counts from model (computed in constructor)
        int n_quad = model->GetNumQuadraticLoadings();
        int n_inter = model->GetNumInteractionLoadings();

        int n_params;
        if (mtype == ModelType::LOGIT && n_choice > 2) {
            // Multinomial logit: each non-reference choice has its own parameters
            n_params = (n_choice - 1) * (regressor_idx.size() + n_free_loadings + n_quad + n_inter);
            // For multinomial logit with use_types, each choice gets (n_types - 1) type-specific intercepts
            if (comp_use_types && n_types > 1) {
                n_params += (n_choice - 1) * (n_types - 1);
            }
        } else if (mtype == ModelType::OPROBIT) {
            // Ordered probit: shared coefficients + thresholds
            n_params = regressor_idx.size() + n_free_loadings + n_quad + n_inter + (n_choice - 1);
            // Add type-specific intercepts only if use_types and n_types > 1
            if (comp_use_types && n_types > 1) {
                n_params += (n_types - 1);
            }
        } else {
            // Binary models (linear, probit, binary logit)
            n_params = regressor_idx.size() + n_free_loadings + n_quad + n_inter;
            if (mtype == ModelType::LINEAR) n_params += 1;  // sigma
            // Add type-specific intercepts only if use_types and n_types > 1
            if (comp_use_types && n_types > 1) {
                n_params += (n_types - 1);
            }
        }

        fm->AddModel(model, n_params);
    }

    // Build parameter constraints based on fixed_coefficients from each component
    int total_params = fm->GetNParam();
    std::vector<bool> param_fixed_vec(total_params, false);

    // Track parameter position as we go through components
    // Start after factor variance parameters (and correlation/SE/type params)
    int param_offset = 0;

    // Handle factor structure-specific parameter counts
    // For mixture models (n_mixtures > 1), parameter layout is:
    // - nmix * n_variance_per_mixture variances
    // - (nmix-1) * n_factors_for_mixture mixture means (for non-reference mixtures)
    // - (nmix-1) log-weights (for non-reference mixtures)
    // - SE params (if SE model)
    int n_factors_for_mixture = n_fac;  // Default: all factors have mixture
    if (fac_struct == FactorStructure::SE_LINEAR || fac_struct == FactorStructure::SE_QUADRATIC) {
        n_factors_for_mixture = n_fac - 1;  // Only input factors for SE models
    }

    if (fac_struct == FactorStructure::SE_LINEAR) {
        // SE_linear: (n_fac - 1) input factor variances + intercept + (n_fac-1) linear + [(n_types-1) type intercepts] + residual var
        int n_input_factors = n_fac - 1;
        param_offset = n_input_factors * n_mixtures;  // Input factor variances (per mixture)
        param_offset += (n_mixtures - 1) * n_input_factors;  // Mixture means
        param_offset += (n_mixtures - 1);  // Mixture log-weights
        param_offset += 1;  // SE intercept
        param_offset += n_input_factors;  // SE linear coefficients
        if (n_types > 1) {
            param_offset += (n_types - 1);  // Type-specific SE intercepts
        }
        param_offset += 1;  // SE residual variance
    } else if (fac_struct == FactorStructure::SE_QUADRATIC) {
        // SE_quadratic: (n_fac - 1) input factor variances + intercept + (n_fac-1) linear + (n_fac-1) quadratic + [(n_types-1) type intercepts] + residual var
        int n_input_factors = n_fac - 1;
        param_offset = n_input_factors * n_mixtures;  // Input factor variances (per mixture)
        param_offset += (n_mixtures - 1) * n_input_factors;  // Mixture means
        param_offset += (n_mixtures - 1);  // Mixture log-weights
        param_offset += 1;  // SE intercept
        param_offset += n_input_factors;  // SE linear coefficients
        param_offset += n_input_factors;  // SE quadratic coefficients
        if (n_types > 1) {
            param_offset += (n_types - 1);  // Type-specific SE intercepts
        }
        param_offset += 1;  // SE residual variance
    } else if (fac_struct == FactorStructure::CORRELATION && n_fac == 2) {
        // Correlated 2-factor: (2 variances + 1 correlation) * nmix + mixture params
        int n_var_per_mix = n_fac * (n_fac + 1) / 2;  // variances + correlation
        param_offset = n_var_per_mix * n_mixtures;
        param_offset += (n_mixtures - 1) * n_fac;  // Mixture means
        param_offset += (n_mixtures - 1);  // Mixture log-weights
    } else {
        // Independent factors
        param_offset = n_fac * n_mixtures;  // Factor variances (per mixture)
        param_offset += (n_mixtures - 1) * n_fac;  // Mixture means
        param_offset += (n_mixtures - 1);  // Mixture log-weights
        // Add correlation parameter offset if present (backward compatibility)
        if (correlation && n_fac == 2 && n_mixtures == 1) {
            param_offset += 1;
        }
    }

    // Parameter-offset layout MUST mirror the actual FactorModel parameter
    // layout produced by the constructor + Set{FactorMean,SE}Covariates +
    // AddModel call sequence in initialize_factor_model_cpp:
    //
    //   factor_var* -> se_* -> typeprob/type_loading -> factor_mean*
    //   -> se_cov* -> model params
    //
    // The constructor (FactorModel.cpp lines ~149-153) appends typeprob/
    // type_loading immediately after the SE block; SetFactorMeanCovariates
    // and SetSECovariates append AFTER typeprob. Computing param_offset in
    // a different order desyncs param_fixed_vec (the constraint vector
    // built here) from the actual parameter positions, silently fixing the
    // wrong slots (e.g., the outcome-factor type loading being marked at
    // an se_cov index).

    // Add type model parameters offset if n_types > 1 and at least one component uses types
    // Type model: log(P(type=t)/P(type=1)) = typeprob_t_intercept + sum_k lambda_t_k * f_k
    // Parameters: (n_types - 1) intercepts + (n_types - 1) * n_fac loadings
    int type_param_start = param_offset;  // capture start BEFORE adding type params
    if (n_types > 1 && any_uses_types) {
        param_offset += (n_types - 1) + (n_types - 1) * n_fac;  // Type intercepts + Type loadings
    }

    // Add factor mean covariate parameters offset if specified
    if (factor_model.containsElementNamed("factor_covariates") &&
        !Rf_isNull(factor_model["factor_covariates"])) {
        CharacterVector factor_cov_names = factor_model["factor_covariates"];
        int n_cov = factor_cov_names.size();
        if (n_cov > 0) {
            // Determine how many factors get mean covariates
            int n_factors_with_mean;
            if (fac_struct == FactorStructure::SE_LINEAR || fac_struct == FactorStructure::SE_QUADRATIC) {
                n_factors_with_mean = n_fac - 1;  // Only input factors
            } else {
                n_factors_with_mean = n_fac;
            }
            param_offset += n_factors_with_mean * n_cov;
        }
    }

    // Add SE covariate parameters offset if specified
    if (factor_model.containsElementNamed("se_covariates") &&
        !Rf_isNull(factor_model["se_covariates"])) {
        CharacterVector se_cov_names = factor_model["se_covariates"];
        int n_se_cov = se_cov_names.size();
        if (n_se_cov > 0) {
            param_offset += n_se_cov;
        }
    }

    // Check if previous_stage_info exists - if so, fix all factor-level parameters
    // (factor variances, correlations, SE params, type params)
    // EXCEPTION: If allow_different_structure is TRUE, the factor-level parameters
    // should remain FREE because Stage 2 uses a different factor structure (SE_linear/SE_quadratic)
    if (model_system.containsElementNamed("previous_stage_info") &&
        !Rf_isNull(model_system["previous_stage_info"])) {
        List prev_stage_info = model_system["previous_stage_info"];
        bool allow_diff_struct = false;
        if (prev_stage_info.containsElementNamed("allow_different_structure") &&
            !Rf_isNull(prev_stage_info["allow_different_structure"])) {
            allow_diff_struct = as<bool>(prev_stage_info["allow_different_structure"]);
        }

        if (!allow_diff_struct) {
            // Standard case: fix all factor-level parameters from previous stage
            for (int j = 0; j < param_offset; j++) {
                param_fixed_vec[j] = true;
            }

            // Selectively un-fix parameters listed in free_param_names.
            // This enables "fix all measurement params but free factor_var"
            // workflows via define_model_system(free_params = c("factor_var_1")).
            if (prev_stage_info.containsElementNamed("free_param_names") &&
                !Rf_isNull(prev_stage_info["free_param_names"])) {
                CharacterVector free_names = prev_stage_info["free_param_names"];
                if (free_names.size() > 0) {
                    // Build a name→index map for factor-level params (indices 0..param_offset-1)
                    // using the SAME layout that param_offset was computed with (and that
                    // the FactorModel constructor + Set{FactorMean,SE}Covariates produce):
                    //   factor_var* -> [mix_means, mix_logweights] -> se_* / corr ->
                    //   typeprob/type_loading -> factor_mean* -> se_cov*
                    //
                    // The map MUST include every factor-distribution parameter type that
                    // a caller could legitimately list in `free_params`. Missing names
                    // (e.g., typeprob_*_intercept, type_*_loading_*, factor_mean_*_*,
                    // se_cov_*) silently fail the un-fix lookup and the C++ side keeps
                    // them frozen while R's setup_parameter_constraints leaves them
                    // free, scrambling the gradient/Hessian/estimates mapping.
                    std::map<std::string, int> fac_name_idx;
                    int idx = 0;
                    bool is_se = (fac_struct == FactorStructure::SE_LINEAR ||
                                  fac_struct == FactorStructure::SE_QUADRATIC);
                    int n_var_factors = is_se ? (n_fac - 1) : n_fac;

                    // Block 1: factor variances (per mixture)
                    for (int m = 0; m < n_mixtures; m++) {
                        for (int k = 0; k < n_var_factors; k++) {
                            std::string nm = (n_mixtures == 1)
                                ? "factor_var_" + std::to_string(k + 1)
                                : "mix" + std::to_string(m + 1) + "_factor_var_" + std::to_string(k + 1);
                            fac_name_idx[nm] = idx++;
                        }
                    }

                    // Correlation parameter (correlation structure only, n_mixtures == 1)
                    if (fac_struct == FactorStructure::CORRELATION && n_fac == 2) {
                        fac_name_idx["factor_corr_1_2"] = idx++;
                    }

                    // Block 2: mixture means + log-weights (for n_mixtures > 1)
                    for (int m = 0; m < n_mixtures - 1; m++) {
                        for (int k = 0; k < n_var_factors; k++)
                            fac_name_idx["mix" + std::to_string(m+1) + "_factor_mean_" + std::to_string(k+1)] = idx++;
                    }
                    for (int m = 0; m < n_mixtures - 1; m++)
                        fac_name_idx["mix" + std::to_string(m+1) + "_logweight"] = idx++;

                    // Block 3: SE parameters (SE structures only)
                    if (is_se) {
                        fac_name_idx["se_intercept"] = idx++;
                        for (int k = 0; k < n_var_factors; k++)
                            fac_name_idx["se_linear_" + std::to_string(k+1)] = idx++;
                        if (fac_struct == FactorStructure::SE_QUADRATIC) {
                            for (int k = 0; k < n_var_factors; k++)
                                fac_name_idx["se_quadratic_" + std::to_string(k+1)] = idx++;
                        }
                        if (n_types > 1) {
                            for (int t = 2; t <= n_types; t++)
                                fac_name_idx["se_intercept_type_" + std::to_string(t)] = idx++;
                        }
                        fac_name_idx["se_residual_var"] = idx++;
                    }

                    // Block 4: type-probability params (typeprob intercepts + type loadings)
                    if (n_types > 1 && any_uses_types) {
                        for (int t = 2; t <= n_types; t++)
                            fac_name_idx["typeprob_" + std::to_string(t) + "_intercept"] = idx++;
                        for (int t = 2; t <= n_types; t++) {
                            for (int k = 1; k <= n_fac; k++) {
                                fac_name_idx["type_" + std::to_string(t) + "_loading_" + std::to_string(k)] = idx++;
                            }
                        }
                    }

                    // Block 5: factor-mean covariate params
                    if (factor_model.containsElementNamed("factor_covariates") &&
                        !Rf_isNull(factor_model["factor_covariates"])) {
                        CharacterVector fcov_names = factor_model["factor_covariates"];
                        int n_fcov = fcov_names.size();
                        if (n_fcov > 0) {
                            int n_fac_with_mean = is_se ? (n_fac - 1) : n_fac;
                            for (int k = 1; k <= n_fac_with_mean; k++) {
                                for (int j = 0; j < n_fcov; j++) {
                                    std::string cname = as<std::string>(fcov_names[j]);
                                    fac_name_idx["factor_mean_" + std::to_string(k) + "_" + cname] = idx++;
                                }
                            }
                        }
                    }

                    // Block 6: SE covariate params
                    if (factor_model.containsElementNamed("se_covariates") &&
                        !Rf_isNull(factor_model["se_covariates"])) {
                        CharacterVector secov_names = factor_model["se_covariates"];
                        int n_secov = secov_names.size();
                        for (int j = 0; j < n_secov; j++) {
                            std::string cname = as<std::string>(secov_names[j]);
                            fac_name_idx["se_cov_" + cname] = idx++;
                        }
                    }

                    // Un-fix each free_param_name. Names not present in the map
                    // (e.g., a user-typo or a measurement-system param) are silently
                    // ignored; only factor-distribution params can be un-fixed here.
                    for (int j = 0; j < free_names.size(); j++) {
                        std::string fn = std::string(free_names[j]);
                        auto it = fac_name_idx.find(fn);
                        if (it != fac_name_idx.end() && it->second < param_offset) {
                            param_fixed_vec[it->second] = false;
                        }
                    }
                }
            }
        }
        // When allow_diff_struct is true, factor-level parameters remain FREE
        // and will use the new factor structure (SE_linear/SE_quadratic)
    }

    // Process each component's fixed coefficients and all_params_fixed
    for (int i = 0; i < components.size(); i++) {
        List comp = components[i];

        std::string model_type_str = as<std::string>(comp["model_type"]);
        std::vector<std::string> covariate_names;
        if (!Rf_isNull(comp["covariates"])) {
            covariate_names = as<std::vector<std::string>>(comp["covariates"]);
        }

        int n_choice = comp.containsElementNamed("num_choices") ?
                      int(comp["num_choices"]) : 2;

        // Count free loadings for this component
        int n_free_loadings = 0;
        if (comp.containsElementNamed("loading_normalization")) {
            SEXP norm_sexp = comp["loading_normalization"];
            if (!Rf_isNull(norm_sexp)) {
                NumericVector norm_vec = as<NumericVector>(norm_sexp);
                for (int j = 0; j < norm_vec.size(); j++) {
                    if (NumericVector::is_na(norm_vec[j])) {
                        n_free_loadings++;
                    }
                }
            }
        } else {
            n_free_loadings = n_fac;
        }

        // Count quadratic and interaction loadings
        int n_quad = 0, n_inter = 0;
        if (comp.containsElementNamed("factor_spec")) {
            std::string fs = as<std::string>(comp["factor_spec"]);
            if (fs == "quadratic" || fs == "full") {
                n_quad = n_free_loadings;  // One quadratic term per free linear loading
            }
            if (fs == "interactions" || fs == "full") {
                // Number of interaction terms = n_fac choose 2
                n_inter = n_fac * (n_fac - 1) / 2;
            }
        }

        // Extract use_types for this component
        bool comp_use_types = false;
        if (comp.containsElementNamed("use_types") && !Rf_isNull(comp["use_types"])) {
            comp_use_types = as<bool>(comp["use_types"]);
        }

        // Calculate n_params for this component (needed for both all_params_fixed and offset)
        int n_params_comp;
        ModelType mtype;
        if (model_type_str == "linear") mtype = ModelType::LINEAR;
        else if (model_type_str == "probit") mtype = ModelType::PROBIT;
        else if (model_type_str == "logit") mtype = ModelType::LOGIT;
        else mtype = ModelType::OPROBIT;

        if (mtype == ModelType::LOGIT && n_choice > 2) {
            n_params_comp = (n_choice - 1) * (covariate_names.size() + n_free_loadings + n_quad + n_inter);
            if (comp_use_types && n_types > 1) n_params_comp += (n_choice - 1) * (n_types - 1);
        } else if (mtype == ModelType::OPROBIT) {
            n_params_comp = covariate_names.size() + n_free_loadings + n_quad + n_inter + (n_choice - 1);
            if (comp_use_types && n_types > 1) n_params_comp += (n_types - 1);
        } else {
            n_params_comp = covariate_names.size() + n_free_loadings + n_quad + n_inter;
            if (mtype == ModelType::LINEAR) n_params_comp += 1;
            if (comp_use_types && n_types > 1) n_params_comp += (n_types - 1);
        }

        // Check if all parameters are fixed for this component (multi-stage estimation)
        bool all_params_fixed = false;
        if (comp.containsElementNamed("all_params_fixed")) {
            all_params_fixed = as<bool>(comp["all_params_fixed"]);
        }

        if (all_params_fixed) {
            // Mark ALL parameters for this component as fixed
            for (int j = 0; j < n_params_comp; j++) {
                int param_idx = param_offset + j;
                if (param_idx < total_params) {
                    param_fixed_vec[param_idx] = true;
                }
            }
        } else {
            // For oprobit models, fix the intercept (absorbed into thresholds)
            if (model_type_str == "oprobit") {
                // Find if 'intercept' is in covariates
                for (int k = 0; k < covariate_names.size(); k++) {
                    if (covariate_names[k] == "intercept") {
                        int param_idx = param_offset + k;
                        if (param_idx < total_params) {
                            param_fixed_vec[param_idx] = true;
                        }
                        break;
                    }
                }
            }

            // Check for individual fixed_coefficients in this component
            if (comp.containsElementNamed("fixed_coefficients") &&
                !Rf_isNull(comp["fixed_coefficients"])) {
                List fixed_coefs = comp["fixed_coefficients"];

                for (int fc_idx = 0; fc_idx < fixed_coefs.size(); fc_idx++) {
                    List fc = fixed_coefs[fc_idx];
                    std::string cov_name = as<std::string>(fc["covariate"]);

                    // Find covariate position
                    int cov_pos = -1;
                    for (int k = 0; k < covariate_names.size(); k++) {
                        if (covariate_names[k] == cov_name) {
                            cov_pos = k;
                            break;
                        }
                    }

                    if (cov_pos >= 0) {
                        // Determine parameter index based on model type and choice
                        int param_idx;

                        if (model_type_str == "logit" && n_choice > 2) {
                            // Multinomial logit: check if choice is specified
                            int choice = 1;  // Default to first non-reference choice
                            if (fc.containsElementNamed("choice") && !Rf_isNull(fc["choice"])) {
                                choice = as<int>(fc["choice"]);
                            }
                            // Each choice has: covariates + loadings + quad + inter
                            int params_per_choice = covariate_names.size() + n_free_loadings + n_quad + n_inter;
                            param_idx = param_offset + (choice - 1) * params_per_choice + cov_pos;
                        } else {
                            // Binary/linear/probit/oprobit: coefficients come first
                            param_idx = param_offset + cov_pos;
                        }

                        // Mark as fixed
                        if (param_idx < total_params) {
                            param_fixed_vec[param_idx] = true;
                        }
                    }
                }
            }

            // Handle fixed_type_intercepts
            if (comp.containsElementNamed("fixed_type_intercepts") &&
                !Rf_isNull(comp["fixed_type_intercepts"])) {
                List fixed_type_ints = comp["fixed_type_intercepts"];

                // Calculate base parameters (everything except type intercepts)
                // Type intercepts are at the end of the component's parameter block
                int n_type_intercepts = (n_types > 1) ? (n_types - 1) : 0;
                int base_params = n_params_comp - n_type_intercepts;

                for (int fti_idx = 0; fti_idx < fixed_type_ints.size(); fti_idx++) {
                    List fti = fixed_type_ints[fti_idx];
                    int type_num = as<int>(fti["type"]);  // 2-indexed (type 2, 3, ...)

                    // Type intercept index within component: base_params + (type_num - 2)
                    // type_num is 2..n_types, so offset is 0..(n_types-2)
                    int type_intercept_offset = type_num - 2;  // 0-indexed
                    int param_idx = param_offset + base_params + type_intercept_offset;

                    // Mark as fixed
                    if (param_idx >= 0 && param_idx < total_params) {
                        param_fixed_vec[param_idx] = true;
                    }
                }
            }
        }

        // Advance param_offset for next component
        param_offset += n_params_comp;
    }

    // Mark unidentified factor variances as fixed
    // Factor variances are the first n_fac parameters (for independent/correlation structures)
    // or first (n_fac - 1) for SE structures
    if (fac_struct == FactorStructure::INDEPENDENT || fac_struct == FactorStructure::CORRELATION) {
        for (int k = 0; k < n_fac; k++) {
            if (!factor_identified[k]) {
                param_fixed_vec[k] = true;  // Fix this factor's variance
            }
        }
    } else if (fac_struct == FactorStructure::SE_LINEAR || fac_struct == FactorStructure::SE_QUADRATIC) {
        // For SE models, only input factors (first n_fac - 1) have variance parameters
        for (int k = 0; k < n_fac - 1; k++) {
            if (!factor_identified[k]) {
                param_fixed_vec[k] = true;
            }
        }

        // Auto-fix type loadings on the outcome factor to 0.
        // Rationale: type probabilities must be a function of the INPUT factors only —
        // the outcome factor is a deterministic function of inputs + type, so a loading
        // on it would create a circular dependency. We silently fix these loadings to 0.
        // The R side errors if the user explicitly set a non-zero initial value.
        //
        // Layout of type params (when n_types > 1 && any_uses_types):
        //   type_param_start + 0 .. (n_types - 2)                       : type_t_intercept (t = 2..n_types)
        //   type_param_start + (n_types - 1) + (t-2)*n_fac + (k-1)      : type_t_loading_k
        if (n_types > 1 && any_uses_types) {
            int type_loading_block_start = type_param_start + (n_types - 1);
            int outcome_k0 = n_fac - 1;  // 0-indexed outcome factor
            for (int t_off = 0; t_off < n_types - 1; t_off++) {
                int idx = type_loading_block_start + t_off * n_fac + outcome_k0;
                param_fixed_vec[idx] = true;
            }
        }
    }

    // Handle equality constraints - mark tied parameters as fixed and build equality mapping
    // Equality constraints tie parameters so that tied params = primary param during optimization
    // C++ handles gradient/Hessian aggregation via equality_mapping in ExtractFreeGradient/Hessian
    if (model_system.containsElementNamed("equality_constraints") &&
        !Rf_isNull(model_system["equality_constraints"])) {
        List eq_constraints = model_system["equality_constraints"];

        // Build parameter name to index mapping
        std::map<std::string, int> param_name_to_idx;

        // Determine n_factors_for_mixture for parameter naming
        int n_factors_for_mixture_eq = n_fac;
        if (fac_struct == FactorStructure::SE_LINEAR || fac_struct == FactorStructure::SE_QUADRATIC) {
            n_factors_for_mixture_eq = n_fac - 1;
        }

        // Factor-level parameters
        // Layout MUST mirror the FactorModel constructor + SetFactorMeanCovariates
        // + SetSECovariates + AddModel call sequence in initialize_factor_model_cpp:
        //
        //   factor_var* -> [factor_corr_*] -> mix_factor_mean_* / mix_logweight
        //   -> se_* (intercept, linear, [quadratic], type intercepts, residual var)
        //   -> typeprob_*_intercept / type_*_loading_*
        //   -> factor_mean_<k>_<cov>
        //   -> se_cov_<cov>
        //   -> component model params (handled by the loop below)
        //
        // Earlier versions of this map were missing typeprob_*_intercept (it was
        // misspelled as type_*_intercept), all type_*_loading_*, and the entire
        // factor_mean_*_* and se_cov_* blocks. Beyond making those equality
        // constraints unrecognised, the missing factor_mean / se_cov slots also
        // shifted every component-level idx that follows backwards, so equality
        // constraints on loadings/sigmas/thresholds got mapped to factor-level
        // slots whenever factor_covariates or se_covariates were used.
        int idx = 0;
        bool is_se_eq = (fac_struct == FactorStructure::SE_LINEAR ||
                         fac_struct == FactorStructure::SE_QUADRATIC);
        int n_var_factors_eq = is_se_eq ? (n_fac - 1) : n_fac;

        // Block 1: factor variances (per mixture) + correlation params
        for (int m = 0; m < n_mixtures; m++) {
            for (int k = 0; k < n_var_factors_eq; k++) {
                std::string nm = (n_mixtures == 1)
                    ? "factor_var_" + std::to_string(k + 1)
                    : "mix" + std::to_string(m + 1) + "_factor_var_" + std::to_string(k + 1);
                param_name_to_idx[nm] = idx++;
            }
            // Correlation params interleave per mixture in the CORRELATION layout
            if (fac_struct == FactorStructure::CORRELATION) {
                for (int j = 0; j < n_fac - 1; j++) {
                    for (int k = j + 1; k < n_fac; k++) {
                        std::string nm = (n_mixtures == 1)
                            ? "factor_corr_" + std::to_string(j + 1) + "_" + std::to_string(k + 1)
                            : "mix" + std::to_string(m + 1) + "_factor_corr_" + std::to_string(j + 1) + "_" + std::to_string(k + 1);
                        param_name_to_idx[nm] = idx++;
                    }
                }
            }
        }

        // Block 2: mixture means + log-weights (n_mixtures > 1)
        for (int m = 0; m < n_mixtures - 1; m++) {
            for (int k = 0; k < n_var_factors_eq; k++) {
                param_name_to_idx["mix" + std::to_string(m + 1) + "_factor_mean_" + std::to_string(k + 1)] = idx++;
            }
        }
        for (int m = 0; m < n_mixtures - 1; m++) {
            param_name_to_idx["mix" + std::to_string(m + 1) + "_logweight"] = idx++;
        }

        // Block 3: SE parameters (SE structures only)
        if (is_se_eq) {
            param_name_to_idx["se_intercept"] = idx++;
            for (int k = 0; k < n_var_factors_eq; k++) {
                param_name_to_idx["se_linear_" + std::to_string(k + 1)] = idx++;
            }
            if (fac_struct == FactorStructure::SE_QUADRATIC) {
                for (int k = 0; k < n_var_factors_eq; k++) {
                    param_name_to_idx["se_quadratic_" + std::to_string(k + 1)] = idx++;
                }
            }
            if (n_types > 1) {
                for (int t = 2; t <= n_types; t++) {
                    param_name_to_idx["se_intercept_type_" + std::to_string(t)] = idx++;
                }
            }
            param_name_to_idx["se_residual_var"] = idx++;
        }

        // Block 4: typeprob_*_intercept + type_*_loading_* (when n_types > 1
        // AND any component uses types — matches the FactorModel constructor's
        // gate, which is `if (ntyp > 1)` only adds these slots).
        if (n_types > 1 && any_uses_types) {
            for (int t = 2; t <= n_types; t++) {
                param_name_to_idx["typeprob_" + std::to_string(t) + "_intercept"] = idx++;
            }
            for (int t = 2; t <= n_types; t++) {
                for (int k = 1; k <= n_fac; k++) {
                    param_name_to_idx["type_" + std::to_string(t) + "_loading_" + std::to_string(k)] = idx++;
                }
            }
        }

        // Block 5: factor_mean_<k>_<cov>
        if (factor_model.containsElementNamed("factor_covariates") &&
            !Rf_isNull(factor_model["factor_covariates"])) {
            CharacterVector fcov_names_eq = factor_model["factor_covariates"];
            int n_fcov_eq = fcov_names_eq.size();
            if (n_fcov_eq > 0) {
                int n_fac_with_mean = is_se_eq ? (n_fac - 1) : n_fac;
                for (int k = 1; k <= n_fac_with_mean; k++) {
                    for (int j = 0; j < n_fcov_eq; j++) {
                        std::string cname = as<std::string>(fcov_names_eq[j]);
                        param_name_to_idx["factor_mean_" + std::to_string(k) + "_" + cname] = idx++;
                    }
                }
            }
        }

        // Block 6: se_cov_<cov>
        if (factor_model.containsElementNamed("se_covariates") &&
            !Rf_isNull(factor_model["se_covariates"])) {
            CharacterVector secov_names_eq = factor_model["se_covariates"];
            int n_secov_eq = secov_names_eq.size();
            for (int j = 0; j < n_secov_eq; j++) {
                std::string cname = as<std::string>(secov_names_eq[j]);
                param_name_to_idx["se_cov_" + cname] = idx++;
            }
        }

        // Component-level parameters (must iterate through components in order)
        for (int i = 0; i < components.size(); i++) {
            List comp = components[i];
            std::string comp_name = as<std::string>(comp["name"]);
            std::string model_type_str = as<std::string>(comp["model_type"]);

            CharacterVector covariate_names;
            if (comp.containsElementNamed("covariates") && !Rf_isNull(comp["covariates"])) {
                covariate_names = comp["covariates"];
            }

            // Get number of free loadings
            int n_free_loadings = 0;
            std::vector<int> free_loading_indices;
            if (comp.containsElementNamed("loading_normalization") && !Rf_isNull(comp["loading_normalization"])) {
                NumericVector norm_vec = comp["loading_normalization"];
                for (int j = 0; j < norm_vec.size(); j++) {
                    if (NumericVector::is_na(norm_vec[j])) {
                        n_free_loadings++;
                        free_loading_indices.push_back(j + 1);
                    }
                }
            }

            int num_choices = 2;  // default for binary
            if (comp.containsElementNamed("num_choices") && !Rf_isNull(comp["num_choices"])) {
                num_choices = as<int>(comp["num_choices"]);
            }

            // Handle different model types
            if (model_type_str == "logit" && num_choices > 2) {
                // Multinomial logit - params per choice
                for (int c = 1; c < num_choices; c++) {
                    std::string choice_suffix = "_c" + std::to_string(c);
                    for (int j = 0; j < covariate_names.size(); j++) {
                        param_name_to_idx[comp_name + choice_suffix + "_" + as<std::string>(covariate_names[j])] = idx++;
                    }
                    for (int j = 0; j < free_loading_indices.size(); j++) {
                        param_name_to_idx[comp_name + choice_suffix + "_loading_" + std::to_string(free_loading_indices[j])] = idx++;
                    }
                    // Type intercepts for this choice
                    if (n_types > 1) {
                        for (int t = 2; t <= n_types; t++) {
                            param_name_to_idx[comp_name + choice_suffix + "_type_" + std::to_string(t) + "_intercept"] = idx++;
                        }
                    }
                }
            } else {
                // Standard component (linear, probit, oprobit, binary logit)
                for (int j = 0; j < covariate_names.size(); j++) {
                    param_name_to_idx[comp_name + "_" + as<std::string>(covariate_names[j])] = idx++;
                }
                for (int j = 0; j < free_loading_indices.size(); j++) {
                    param_name_to_idx[comp_name + "_loading_" + std::to_string(free_loading_indices[j])] = idx++;
                }

                // Linear model has sigma
                if (model_type_str == "linear") {
                    param_name_to_idx[comp_name + "_sigma"] = idx++;
                }

                // Ordered probit has thresholds. The component's field is
                // named "num_choices"; an earlier version of this code
                // checked "n_categories" here, which never matched, so
                // threshold names were never added to param_name_to_idx.
                // That in turn caused equality constraints referencing
                // _thresh_k names to silently be skipped (the name
                // lookup at line ~1132 returned end()).
                if (model_type_str == "oprobit" && comp.containsElementNamed("num_choices")) {
                    int n_cat = as<int>(comp["num_choices"]);
                    for (int j = 1; j < n_cat; j++) {
                        param_name_to_idx[comp_name + "_thresh_" + std::to_string(j)] = idx++;
                    }
                }

                // Type intercepts
                if (n_types > 1) {
                    for (int t = 2; t <= n_types; t++) {
                        param_name_to_idx[comp_name + "_type_" + std::to_string(t) + "_intercept"] = idx++;
                    }
                }
            }
        }

        // Process each equality constraint - mark tied params as fixed and build mapping
        // equality_mapping[i] = j means param i is tied to param j (j is primary)
        // equality_mapping[i] = -1 means param i is not tied
        std::vector<int> equality_map(total_params, -1);

        for (int i = 0; i < eq_constraints.size(); i++) {
            CharacterVector constraint = eq_constraints[i];
            if (constraint.size() >= 2) {
                // First param is primary (free), rest are tied (fixed)
                std::string primary_name = as<std::string>(constraint[0]);
                auto primary_it = param_name_to_idx.find(primary_name);
                int primary_idx = -1;
                if (primary_it != param_name_to_idx.end()) {
                    primary_idx = primary_it->second;
                }

                for (int j = 1; j < constraint.size(); j++) {
                    std::string tied_name = as<std::string>(constraint[j]);
                    auto it = param_name_to_idx.find(tied_name);
                    if (it != param_name_to_idx.end()) {
                        param_fixed_vec[it->second] = true;
                        if (primary_idx >= 0) {
                            equality_map[it->second] = primary_idx;
                        }
                    }
                }
            }
        }

        // Apply user-fixed factor-distribution parameters from
        // factor$fixed_params (set via fix_factor_param() in R). Reuses the
        // same constructor-true layout the un-fix loop and the
        // equality_constraints map use.
        apply_fix_factor_param(factor_model, fac_struct, n_fac, n_types,
                               n_mixtures, any_uses_types, param_offset,
                               param_fixed_vec,
                               init_params /* may be NULL */);

        // Set parameter constraints with optional initial values
        if (init_params.isNotNull()) {
            NumericVector ip(init_params);
            std::vector<double> init_params_vec = as<std::vector<double>>(ip);
            fm->SetParameterConstraints(param_fixed_vec, init_params_vec);
        } else {
            fm->SetParameterConstraints(param_fixed_vec);
        }

        // Set equality constraints mapping
        fm->SetEqualityConstraints(equality_map);
    } else {
        // No equality constraints - just set parameter constraints
        apply_fix_factor_param(factor_model, fac_struct, n_fac, n_types,
                               n_mixtures, any_uses_types, param_offset,
                               param_fixed_vec,
                               init_params /* may be NULL */);

        if (init_params.isNotNull()) {
            NumericVector ip(init_params);
            std::vector<double> init_params_vec = as<std::vector<double>>(ip);
            fm->SetParameterConstraints(param_fixed_vec, init_params_vec);
        } else {
            fm->SetParameterConstraints(param_fixed_vec);
        }
    }

    return fm;
}

//' Evaluate log-likelihood for given parameters
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param params Vector of parameters
//' @param compute_gradient Whether to compute gradient (default FALSE)
//' @param compute_hessian Whether to compute Hessian (default FALSE)
//' @return List with:
//'   - logLikelihood: scalar log-likelihood value
//'   - gradient: vector of length n_param_free (if requested)
//'   - hessian: vector of length n_param_free*(n_param_free+1)/2 stored as
//'              upper-triangular in row-major order (if requested).
//'              To expand to full symmetric matrix in R:
//'              \code{idx <- 1; for(i in 1:n) for(j in i:n) { H[i,j] <- H[j,i] <- hess[idx]; idx <- idx + 1 }}
//' @export
// [[Rcpp::export]]
List evaluate_likelihood_cpp(SEXP fm_ptr, NumericVector params,
                             bool compute_gradient = false,
                             bool compute_hessian = false) {

    // Get FactorModel object
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    // Convert parameters to std::vector
    std::vector<double> params_vec = as<std::vector<double>>(params);

    // Determine flag
    int iflag = 1;  // likelihood only
    if (compute_gradient) iflag = 2;
    if (compute_hessian) iflag = 3;

    // Compute likelihood
    double logLkhd;
    std::vector<double> gradL, hessL;
    fm->CalcLkhd(params_vec, logLkhd, gradL, hessL, iflag);

    // Return results
    List result = List::create(Named("logLikelihood") = logLkhd);

    if (compute_gradient) {
        result["gradient"] = gradL;
    }
    if (compute_hessian) {
        result["hessian"] = hessL;
    }

    return result;
}

//' Evaluate log-likelihood only (for optimization)
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param params Vector of parameters
//' @return Log-likelihood value
//' @export
// [[Rcpp::export]]
double evaluate_loglik_only_cpp(SEXP fm_ptr, NumericVector params) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);
    std::vector<double> params_vec = as<std::vector<double>>(params);
    return fm->CalcLogLikelihood(params_vec);
}

//' Get parameter counts from FactorModel
//'
//' @param fm_ptr External pointer to FactorModel object
//' @return List with parameter count information
//' @export
// [[Rcpp::export]]
List get_parameter_info_cpp(SEXP fm_ptr) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    return List::create(
        Named("n_obs") = fm->GetNObs(),
        Named("n_param") = fm->GetNParam(),
        Named("n_param_free") = fm->GetNParamFree()
    );
}

//' Extract free parameters from full parameter vector
//'
//' Given a full parameter vector (including fixed parameters),
//' extract only the free parameters based on the model's fixed parameter mask.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param full_params Full parameter vector (size n_param)
//' @return Vector of free parameters only (size n_param_free)
//' @export
// [[Rcpp::export]]
NumericVector extract_free_params_cpp(SEXP fm_ptr, NumericVector full_params) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    int n_param = fm->GetNParam();
    int n_param_free = fm->GetNParamFree();

    if (full_params.size() != n_param) {
        Rcpp::stop("Full params size (%d) doesn't match expected (%d)",
                   full_params.size(), n_param);
    }

    // Get the fixed parameter mask
    const std::vector<bool>& param_fixed = fm->GetParamFixed();

    // Extract free parameters
    NumericVector free_params(n_param_free);
    int ifree = 0;
    for (int i = 0; i < n_param; i++) {
        if (!param_fixed[i]) {
            free_params[ifree++] = full_params[i];
        }
    }

    return free_params;
}

//' Set observation weights for weighted likelihood estimation
//'
//' Sets per-observation weights for the likelihood calculation. When weights
//' are set, each observation's contribution to the log-likelihood is multiplied
//' by its weight. This is used for importance sampling in adaptive integration.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param weights Numeric vector of observation weights (length = n_obs)
//' @return No return value. Called for its side effect of setting the
//'   per-observation likelihood weights on the FactorModel pointed to by
//'   \code{fm_ptr}.
//' @export
// [[Rcpp::export]]
void set_observation_weights_cpp(SEXP fm_ptr, NumericVector weights) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);
    std::vector<double> weights_vec = as<std::vector<double>>(weights);
    fm->SetObservationWeights(weights_vec);
}

//' Set up adaptive quadrature based on factor scores and standard errors
//'
//' Enables adaptive integration where the number of quadrature points varies
//' by observation based on the precision of factor score estimates. When factor
//' scores are well-determined (small SE), fewer integration points are used.
//' Importance sampling weights are computed automatically.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param factor_scores Matrix (n_obs x n_factors) of factor score estimates
//' @param factor_ses Matrix (n_obs x n_factors) of standard errors
//' @param factor_vars Vector (n_factors) of factor variances from previous stage
//' @param threshold Threshold for determining quadrature points (default 0.5, matching legacy)
//' @param max_quad Maximum quadrature points per factor (default 16)
//' @param verbose Whether to print summary of adaptive quadrature setup (default TRUE)
//' @return No return value. Called for its side effect of enabling adaptive
//'   quadrature on the FactorModel pointed to by \code{fm_ptr} and (when
//'   \code{verbose = TRUE}) printing a summary of the per-observation
//'   integration-point distribution.
//' @export
// [[Rcpp::export]]
void set_adaptive_quadrature_cpp(SEXP fm_ptr,
                                 NumericMatrix factor_scores,
                                 NumericMatrix factor_ses,
                                 NumericVector factor_vars,
                                 double threshold = 0.5,
                                 int max_quad = 16,
                                 bool verbose = true) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    int n_obs = factor_scores.nrow();
    int n_fac = factor_scores.ncol();

    // Convert to nested vectors
    std::vector<std::vector<double>> scores_vec(n_obs);
    std::vector<std::vector<double>> ses_vec(n_obs);
    for (int i = 0; i < n_obs; i++) {
        scores_vec[i].resize(n_fac);
        ses_vec[i].resize(n_fac);
        for (int j = 0; j < n_fac; j++) {
            scores_vec[i][j] = factor_scores(i, j);
            ses_vec[i][j] = factor_ses(i, j);
        }
    }

    std::vector<double> vars_vec = as<std::vector<double>>(factor_vars);

    fm->SetAdaptiveQuadrature(scores_vec, ses_vec, vars_vec, threshold, max_quad);

    // Print summary if verbose
    if (verbose) {
        Rcpp::Rcout << "\nAdaptive Integration Summary\n";
        Rcpp::Rcout << "----------------------------\n";
        Rcpp::Rcout << "Threshold: " << threshold << ", Max quad points: " << max_quad << "\n\n";

        // Compute per-observation total integration points and collect stats
        std::map<int, int> nquad_counts;  // nquad -> count of observations
        double total_points = 0.0;
        double total_points_standard = std::pow(max_quad, n_fac);

        for (int i = 0; i < n_obs; i++) {
            int obs_total = 1;
            for (int j = 0; j < n_fac; j++) {
                // Recompute nquad using same formula as SetAdaptiveQuadrature
                double f_se = factor_ses(i, j);
                double f_sd = std::sqrt(factor_vars[j]);
                double ratio = f_se / f_sd / threshold;
                int nq = 1 + 2 * static_cast<int>(std::floor(ratio));
                if (nq < 1) nq = 1;
                if (nq > max_quad) nq = max_quad;
                if (f_se > f_sd) nq = max_quad;
                obs_total *= nq;
            }
            nquad_counts[obs_total]++;
            total_points += obs_total;
        }

        double avg_points = total_points / n_obs;
        double reduction = 100.0 * (1.0 - avg_points / total_points_standard);

        // Print distribution table
        Rcpp::Rcout << "Integration points per observation:\n";
        Rcpp::Rcout << "  Points   Observations   Percent\n";
        for (auto& kv : nquad_counts) {
            double pct = 100.0 * kv.second / n_obs;
            Rcpp::Rcout << std::setw(8) << kv.first
                        << std::setw(15) << kv.second
                        << std::setw(10) << std::fixed << std::setprecision(1) << pct << "%\n";
        }

        Rcpp::Rcout << "\nAverage integration points: " << std::fixed << std::setprecision(1)
                    << avg_points << " (vs " << static_cast<int>(total_points_standard)
                    << " standard)\n";
        Rcpp::Rcout << "Computational reduction: " << std::fixed << std::setprecision(1)
                    << reduction << "%\n\n";
    }
}

//' Disable adaptive quadrature
//'
//' Reverts to standard (non-adaptive) quadrature integration.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @return No return value. Called for its side effect of disabling adaptive
//'   quadrature and clearing any observation weights on the FactorModel
//'   pointed to by \code{fm_ptr}.
//' @export
// [[Rcpp::export]]
void disable_adaptive_quadrature_cpp(SEXP fm_ptr) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);
    fm->DisableAdaptiveQuadrature();
}

//' Evaluate log-likelihood for a single observation at given factor values
//'
//' Used for factor score estimation. The model parameters are held fixed,
//' and the factor values are treated as the parameters to optimize.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param iobs Observation index (0-based)
//' @param factor_values Vector of factor values (size n_factors)
//' @param model_params Vector of ALL model parameters (from previous estimation)
//' @param compute_gradient Whether to compute gradient (default FALSE)
//' @param compute_hessian Whether to compute Hessian (default FALSE)
//' @param include_prior Whether to include factor prior in likelihood (default TRUE).
//'        Set to FALSE to match legacy C++ behavior (observation likelihood only).
//' @return List with log-likelihood, gradient (if requested), and Hessian (if requested)
//' @export
// [[Rcpp::export]]
List evaluate_factorscore_likelihood_cpp(SEXP fm_ptr,
                                         int iobs,
                                         NumericVector factor_values,
                                         NumericVector model_params,
                                         bool compute_gradient = false,
                                         bool compute_hessian = false,
                                         bool include_prior = true) {

    // Get FactorModel object
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    // Convert to std::vector
    std::vector<double> fac_vec = as<std::vector<double>>(factor_values);
    std::vector<double> param_vec = as<std::vector<double>>(model_params);

    // Determine flag
    int iflag = 1;  // likelihood only
    if (compute_gradient) iflag = 2;
    if (compute_hessian) iflag = 3;

    // Compute likelihood for single observation
    double logLkhd;
    std::vector<double> gradL, hessL;
    fm->CalcLkhdSingleObs(iobs, fac_vec, param_vec, logLkhd, gradL, hessL, iflag, include_prior);

    // Return results
    List result = List::create(Named("logLikelihood") = logLkhd);

    if (compute_gradient) {
        result["gradient"] = gradL;
    }
    if (compute_hessian) {
        result["hessian"] = hessL;
    }

    return result;
}
