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
#'
#' @return A list with components `coefficients`, `loglik`, `iterations` and
#'   `converged` describing the fitted model.
#'
#' @examples
#' \dontrun{
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
                           tol = 1e-6, learning_rate = 0.01) {
  if (!inherits(X, "big.matrix")) {
    stop("`X` must be a big.matrix object", call. = FALSE)
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
  
  big_pls_cox_gd_cpp(X@address, as.numeric(time), as.numeric(status), ncomp,
                     max_iter, tol, learning_rate)
}
