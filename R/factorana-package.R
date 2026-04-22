#' @keywords internal
#' @useDynLib factorana, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats binomial coef glm lm na.omit pnorm sd setNames
#' @importFrom utils head setTxtProgressBar txtProgressBar write.csv
"_PACKAGE"

# Global variables used in foreach loops (to avoid R CMD check notes)
utils::globalVariables(c("obs_idx", ".self_id"))
