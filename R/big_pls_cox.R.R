#' Partial Least Squares Components for Cox Models with Big Matrices
#'
#' @description Compute Partial Least Squares (PLS) components tailored for 
#' Cox proportional hazards models when predictors are stored as a 
#' \code{big.matrix} from the \pkg{bigmemory} package.
#' 
#' @details   The function standardises each predictor column, iteratively 
#' builds latent scores using martingale residuals from Cox fits, and deflates 
#' the predictors without materialising the full design matrix in memory. 
#' Both in-memory and file-backed \pkg{bigmemory} matrices are supported.
#'
#' @param X A numeric matrix or a [`bigmemory::big.matrix`] object containing the predictors.
#' @param time Numeric vector of survival times.
#' @param status Integer (0/1) vector of event indicators.
#' @param ncomp Number of latent components to compute.
#' @param control Optional list passed to [`survival::coxph.control`].
#' @param keepX Optional integer vector specifying the number of variables to
#'   retain (naive sparsity) in each component. A value of zero keeps all
#'   predictors. If a single integer is supplied it is recycled across
#'   components.
#' @return A list with the computed scores, loadings, weights, scaling information and the
#'   fitted Cox model returned by [`survival::coxph.fit`].
#'   
#' @seealso [big_pls_cox_gd()], [predict.big_pls_cox()], [select_ncomp()],
#'   [computeDR()].
#'
#' @references
#'   Maumy, M., Bertrand, F. (2023). PLS models and their extension for big data. 
#'   Joint Statistical Meetings (JSM 2023), Toronto, ON, Canada. 
#'   
#'   Maumy, M., Bertrand, F. (2023). bigPLS: Fitting and cross-validating 
#'   PLS-based Cox models to censored big data. BioC2023 — The Bioconductor 
#'   Annual Conference, Dana-Farber Cancer Institute, Boston, MA, USA. 
#'   Poster. https://doi.org/10.7490/f1000research.1119546.1  
#'
#'   Bastien, P., Bertrand, F., Meyer, N., & Maumy-Bertrand, M. (2015).
#'   Deviance residuals-based sparse PLS and sparse kernel PLS for censored
#'   data. *Bioinformatics*, 31(3), 397–404. <doi:10.1093/bioinformatics/btu660>
#'   
#'   Bertrand, F., Bastien, P., Meyer, N., & Maumy-Bertrand, M. (2014). PLS
#'   models for censored data. In *Proceedings of UseR! 2014* (p. 152).
#'   
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(100), nrow = 20)
#'   time <- rexp(20)
#'   status <- rbinom(20, 1, 0.5)
#'   fit <- big_pls_cox(X, time, status, ncomp = 2)
#'   str(fit)
#' }
#'
big_pls_cox <- function(X, time, status, ncomp = 2L, control = survival::coxph.control(),
                        keepX = NULL) {
  if (!requireNamespace("bigmemory", quietly = TRUE)) {
    stop("Package 'bigmemory' is required")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required")
  }
  
  if (!is.numeric(time) || !is.numeric(status)) {
    stop("time and status must be numeric vectors")
  }
  if (length(time) != length(status)) {
    stop("time and status must have the same length")
  }
  n <- length(time)
  if (!is.matrix(X) && !inherits(X, "big.matrix")) {
    stop("X must be either a numeric matrix or a big.matrix")
  }
  if (is.matrix(X)) {
    X <- bigmemory::as.big.matrix(X)
  }
  if (nrow(X) != n) {
    stop("Number of rows in X must match the length of the response vectors")
  }
  ncomp <- as.integer(ncomp)
  if (length(ncomp) != 1L || ncomp <= 0L) {
    stop("ncomp must be a positive integer")
  }
  p <- ncol(X)
  if (ncomp > min(n, p)) {
    stop("ncomp must be <= min(nrow(X), ncol(X))")
  }
  
  if (is.null(keepX)) {
    keepX <- rep.int(0L, ncomp)
  } else {
    keepX <- as.integer(keepX)
    if (length(keepX) == 1L) {
      keepX <- rep.int(keepX, ncomp)
    }
    if (length(keepX) != ncomp) {
      stop("keepX must be either NULL, a single integer or a vector of length ncomp")
    }
    if (any(is.na(keepX)) || any(keepX < 0L) || any(keepX > p)) {
      stop("keepX entries must lie between 0 and ncol(X)")
    }
  }
  
  if (!inherits(control, "coxph.control")) {
    control <- do.call(survival::coxph.control, control)
  }
  control$iter.max <- max(50L, control$iter.max %||% 20L)
  
  address <- methods::slot(X, "address")
  stats <- big_pls_cox_col_stats_cpp(address)
  means <- stats$mean
  sds <- stats$sd
  
  scores <- matrix(0, nrow = n, ncol = ncomp)
  loadings <- matrix(0, nrow = p, ncol = ncomp)
  weights <- matrix(0, nrow = p, ncol = ncomp)
  
  y <- cbind(time, status)
  strata <- rep(1, n)
  offset <- rep(0, n)
  wt <- rep(1, n)
  rownms <- as.character(seq_len(n))
  
  current_scores <- matrix(numeric(0), nrow = n, ncol = 0)
  current_loadings <- matrix(numeric(0), nrow = p, ncol = 0)
  
  for (h in seq_len(ncomp)) {
    if (ncol(current_scores) == 0L) {
      cox_fit <- survival::coxph.fit(x = current_scores, y = y, strata = strata,
                                     offset = offset, init = numeric(0),
                                     control = control, weights = wt,
                                     method = "efron", rownames = rownms)
    } else {
      init <- rep(0, ncol(current_scores))
      cox_fit <- survival::coxph.fit(x = current_scores, y = y, strata = strata,
                                     offset = offset, init = init,
                                     control = control, weights = wt,
                                     method = "efron", rownames = rownms)
    }
    residuals <- cox_fit$residuals
    component <- big_pls_cox_component_cpp(address, residuals, current_scores,
                                           current_loadings, means, sds,
                                           keepX[h])
    weights[, h] <- component$weights
    scores[, h] <- component$scores
    loadings[, h] <- component$loadings
    
    current_scores <- scores[, seq_len(h), drop = FALSE]
    current_loadings <- loadings[, seq_len(h), drop = FALSE]
  }
  
  final_fit <- survival::coxph.fit(x = current_scores, y = y, strata = strata,
                                   offset = offset, init = rep(0, ncol(current_scores)),
                                   control = control, weights = wt,
                                   method = "efron", rownames = rownms)
  
  structure(list(
    scores = scores,
    loadings = loadings,
    weights = weights,
    center = means,
    scale = sds,
    cox_fit = final_fit,
    keepX = keepX,
    time = as.numeric(time),
    status = as.numeric(status)
  ), class = "big_pls_cox")
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
