#' Gradient-Descent Solver for Cox Models on Big Matrices
#'
#' Fits a Cox proportional hazards regression model using a gradient-descent
#' optimizer implemented in C++. The function operates directly on a
#' [`bigmemory::big.matrix`][bigmemory::big.matrix-class] object to avoid
#' materialising large design matrices in memory.
#'
#' @param X A [`bigmemory::big.matrix`][bigmemory::big.matrix-class] containing
#'   the design matrix (rows are observations).
#' @param time A numeric vector of follow-up times with length equal to the
#'   number of rows of `X`.
#' @param status A numeric or integer vector of the same length as `time`
#'   containing the event indicators (1 for an event, 0 for censoring).
#' @param ncomp An integer giving the number of components (columns) to use from
#'   `X`. Defaults to `min(5, ncol(X))`.
#' @param max_iter Maximum number of gradient-descent iterations (default 500).
#' @param tol Convergence tolerance on the Euclidean distance between successive
#'   coefficient vectors.
#' @param learning_rate Step size used for the gradient-descent updates.
#' @param keepX Optional integer vector describing the number of predictors to
#'   retain per component (naive sparsity). A value of zero keeps all
#'   predictors.
#'
#' @return A list with components:
#' * `coefficients`: Estimated Cox regression coefficients on the latent scores.
#' * `loglik`: Final partial log-likelihood value.
#' * `iterations`: Number of gradient-descent iterations performed.
#' * `converged`: Logical flag indicating whether convergence was achieved.
#' * `scores`: Matrix of latent score vectors (one column per component).
#' * `loadings`: Matrix of loading vectors associated with each component.
#' * `weights`: Matrix of PLS weight vectors.
#' * `center`: Column means used to centre the predictors.
#' * `scale`: Column scales (standard deviations) used to standardise the predictors.
#'
#' @seealso [big_pls_cox()], [predict.big_pls_cox()], [select_ncomp()],
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
#' @examples
#' \donttest{
#' library(bigmemory)
#' set.seed(1)
#' n <- 50
#' p <- 10
#' X <- bigmemory::as.big.matrix(matrix(rnorm(n * p), n, p))
#' time <- rexp(n, rate = 0.1)
#' status <- rbinom(n, 1, 0.7)
#' fit <- big_pls_cox_gd(X, time, status, ncomp = 3, max_iter = 200)
#' }
#' @export
big_pls_cox_gd <- function(X, time, status, ncomp = NULL, max_iter = 500L,
                           tol = 1e-6, learning_rate = 0.01, keepX = NULL) {
  if (!inherits(X, "big.matrix")) {
    stop("`X` must be a big.matrix object", call. = FALSE)
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required", call. = FALSE)
  }
  n <- nrow(X)
  p <- ncol(X)
  if (length(time) != n) {
    stop("`time` must have length equal to the number of rows of `X`", call. = FALSE)
  }
  if (length(status) != n) {
    stop("`status` must have length equal to the number of rows of `X`", call. = FALSE)
  }
  if (!is.numeric(time)) {
    stop("`time` must be numeric", call. = FALSE)
  }
  if (!is.numeric(status)) {
    stop("`status` must be numeric", call. = FALSE)
  }
  if (is.null(ncomp)) {
    ncomp <- min(5L, p)
  }
  ncomp <- as.integer(ncomp)
  if (length(ncomp) != 1 || is.na(ncomp) || ncomp < 1L || ncomp > p) {
    stop("`ncomp` must be a single integer between 1 and ncol(X)", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (length(max_iter) != 1 || is.na(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be a positive integer", call. = FALSE)
  }
  tol <- as.numeric(tol)
  if (length(tol) != 1 || is.na(tol) || tol <= 0) {
    stop("`tol` must be a strictly positive number", call. = FALSE)
  }
  learning_rate <- as.numeric(learning_rate)
  if (length(learning_rate) != 1 || is.na(learning_rate) || learning_rate <= 0) {
    stop("`learning_rate` must be a strictly positive number", call. = FALSE)
  }
  
  if (is.null(keepX)) {
    keep_vec <- rep.int(0L, ncomp)
  } else {
    keep_vec <- as.integer(keepX)
    if (length(keep_vec) == 1L) {
      keep_vec <- rep.int(keep_vec, ncomp)
    }
    if (length(keep_vec) != ncomp) {
      stop("`keepX` must be NULL, a single integer or have length equal to ncomp", call. = FALSE)
    }
    if (any(is.na(keep_vec)) || any(keep_vec < 0L) || any(keep_vec > p)) {
      stop("`keepX` entries must be between 0 and ncol(X)", call. = FALSE)
    }
  }

  result <- big_pls_cox_gd_cpp(X@address, as.numeric(time), as.numeric(status),
                               ncomp, max_iter, tol, learning_rate, keep_vec)

  result$coefficients <- as.numeric(result$coefficients)
  result$center <- as.numeric(result$center)
  result$scale <- as.numeric(result$scale)
  result$keepX <- as.integer(keep_vec)
  result$time <- as.numeric(time)
  result$status <- as.numeric(status)
  class(result) <- "big_pls_cox_gd"
  
  result
}
