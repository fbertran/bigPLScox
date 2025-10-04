#' Incremental Survival Model Fitting with Pre-Scaled Data 
#' 
#' Loads a previously scaled design matrix and continues the stochastic
#' gradient optimisation for a subset of variables.
#' 
#' 
#' @param name.col Character vector containing the column names that should be 
#' included in the partial fit. 
#' @param datapath File system path or connection where the big-memory backing 
#' file for the scaled design matrix is stored.
#' @param ncores Number of processor cores allocated to the partial fitting
#' procedure. Defaults to \code{1}.
#' @param resBigscale Result object returned by \code{\link{bigscale}}
#' containing scaling statistics to be reused. By default the helper reuses the
#' globally cached \code{resultsBigscale} object created by
#' \code{\link{bigscale}}.
#' @param bigmemory.flag Logical flag determining whether big-memory backed
#' matrices are used when loading and updating the design matrix. Defaults to
#' \code{FALSE}.
#' @param parallel.flag Logical flag toggling the use of parallelised
#' stochastic gradient updates. Defaults to \code{FALSE}.
#' @param inf.mth Inference method requested for the partial fit, such as
#' \code{"none"}, \code{"asymptotic"}, or bootstrap summaries. Defaults to
#' \code{"none"}.
#' 
#' @return Either a numeric vector of log hazard-ratio coefficients or, when
#' inference is requested, a matrix whose columns correspond to the inferred
#' coefficient summaries for each penalisation setting.
#' 
#' @seealso \code{\link{bigscale}} and \code{\link{bigSurvSGD.na.omit}}.
#' 
#' @examples
#' \dontrun{
#' continued <- partialbigSurvSGDv0(
#'   name.col = c("age", "sex"),
#'   datapath = tempfile(),
#'   ncores = 2,
#'   resBigscale = scaled,
#'   bigmemory.flag = TRUE,
#'   parallel.flag = TRUE,
#'   inf.mth = "bootstrap"
#' )
#' }
#' 