#' bigPLScox-package
#'
#' Provides Partial least squares Regression for regular, generalized linear and 
#' Cox models for big data. It allows for missing data in the explanatory variables. 
#' Repeated k-fold cross-validation of such models using various criteria. 
#' Bootstrap confidence intervals constructions are also available.
#'
#' @aliases bigPLScox-package bigPLScox NULL
#' 
#' @references Maumy, M., Bertrand, F. (2023). PLS models and their extension for big data. 
#'   Joint Statistical Meetings (JSM 2023), Toronto, ON, Canada. 
#'   
#'   Maumy, M., Bertrand, F. (2023). bigPLS: Fitting and cross-validating 
#'   PLS-based Cox models to censored big data. BioC2023 — The Bioconductor 
#'   Annual Conference, Dana-Farber Cancer Institute, Boston, MA, USA. 
#'   Poster. https://doi.org/10.7490/f1000research.1119546.1  
#' 
#'   Bastien, P., Bertrand, F., Meyer, N., and Maumy-Bertrand, M.
#'   (2015). Deviance residuals-based sparse PLS and sparse kernel PLS for
#'   binary classification and survival analysis. *BMC Bioinformatics*, 16, 211.
#' 
#' @seealso [big_pls_cox()] and [big_pls_cox_gd()]
#' 
#' @examples
#' set.seed(314)
#' library(bigPLScox)
#' data(sim_data)
#' head(sim_data)
#' 
"_PACKAGE"

#' @import bigSurvSGD
# #' @importFrom bigSurvSGD bigSurvSGD lambdaMaxC oneChunkC oneObsPlugingC
# #' @importClassesFrom bigSurvSGD bigSurvSGD
# #' @importMethodsFrom bigSurvSGD plot.bigSurvSGD print.bigSurvSGD
# #' @importFrom utils head
#' @importFrom graphics abline axis barplot layout legend matplot segments
#' @importFrom grDevices dev.new
#' @importFrom stats as.formula 
#' @importFrom stats is.empty.model 
#' @importFrom stats coef
#' @importFrom stats coefficients
#' @importFrom stats cor
#' @importFrom stats complete.cases
#' @importFrom stats extractAIC
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom stats model.matrix 
#' @importFrom stats model.response 
#' @importFrom stats model.weights 
#' @importFrom stats na.omit
#' @importFrom stats predict
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats quantile
#' @importFrom stats residuals
#' @importFrom stats rexp 
#' @importFrom stats rnorm
#' @importFrom stats runif 
#' @importFrom stats sd
#' @importFrom stats uniroot
#' @importFrom stats var
#' @importFrom utils read.csv
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom survival coxph
#' @importFrom survival coxph.control
NULL

#' @useDynLib bigPLScox, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import bigmemory
#' @import bigalgebra
NULL
