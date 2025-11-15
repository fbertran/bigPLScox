#' Transform new data to PLS–Cox scores
#'
#' Project new observations onto previously fitted PLS–Cox components.
#'
#' @param X New data: a numeric matrix or a \code{bigmemory::big.matrix}.
#' @param means Column means used to center the original predictors.
#' @param sds Column standard deviations used to scale the original predictors.
#' @param weights PLS weight matrix (p x ncomp) from a fitted model.
#' @param loadings PLS loading matrix (p x ncomp) from a fitted model.
#' @param comps Integer vector of component indices to return (1-based).
#'
#' @return A numeric matrix of scores with one row per observation in \code{X}
#'   and one column per requested component.
#' @export
big_pls_cox_transform <- function(X,
                                  means,
                                  sds,
                                  weights,
                                  loadings,
                                  comps = seq_len(ncol(weights))) {
  if (is.null(comps) || length(comps) == 0L) {
    stop("'comps' must be a non-empty vector of component indices")
  }
  
  if (!is.numeric(means) || !is.numeric(sds)) {
    stop("'means' and 'sds' must be numeric vectors")
  }
  
  # dispatch on type of X
  if (inherits(X, "big.matrix")) {
    big_pls_cox_transform_cpp(
      xpMat   = X@address,
      means   = means,
      sds     = sds,
      weights = weights,
      loadings = loadings,
      comps   = as.integer(comps)
    )
  } else if (is.matrix(X)) {
    matrix_pls_cox_transform_cpp(
      X       = X,
      means   = means,
      sds     = sds,
      weights = weights,
      loadings = loadings,
      comps   = as.integer(comps)
    )
  } else {
    stop("'X' must be a numeric matrix or a bigmemory::big.matrix")
  }
}