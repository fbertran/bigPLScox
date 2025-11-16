#' Gradient based PLS Cox for big matrices
#'
#' @description
#' Fit a PLS Cox model where the PLS components are computed once and
#' the Cox regression in the latent space is optimised by gradient based
#' methods.
#'
#' This function is intended as a faster, approximate alternative to
#' [big_pls_cox_fast()] when many fits are required or when the design
#' is stored as a [`bigmemory::big.matrix`].
#'
#' @param X A [`bigmemory::big.matrix`][bigmemory::big.matrix-class] containing
#'   the design matrix (rows are observations).
#' @param time A numeric vector of follow-up times with length equal to the
#'   number of rows of `X`.
#' @param status A numeric or integer vector of the same length as `time`
#'   containing the event indicators (1 for an event, 0 for censoring).
#' @param ncomp An integer giving the number of components (columns) to use from
#'   `X`. Defaults to `min(5, ncol(X))`.
#' @param max_iter Maximum number of gradient iterations.
#' @param tol Convergence tolerance for the optimisation in the Cox space.
#'   Both the change in log-likelihood and the Euclidean change in the
#'   coefficient vector are monitored.
#' @param learning_rate Initial learning rate for first order methods.
#'   This is used by `"gd"` and as a starting scale for `"bb"` and
#'   `"nesterov"`. It is ignored by `"bfgs"`.
#' @param keepX Optional integer vector describing the number of predictors to
#'   retain per component (naive sparsity). A value of zero keeps all
#'   predictors.
#' @param coxfit Optional Boolean to fit a Cox model on the extracted components.
#' @param method Optimisation scheme used in the latent space. One of
#'   \describe{
#'     \item{`"gd"`}{Simple fixed step gradient descent. This is the most
#'       transparent method and is useful for debugging and didactic
#'       purposes.}
#'     \item{`"bb"`}{Barzilai Borwein step size. Uses a diagonal
#'       quasi-second-order update of the step size based on the last
#'       two gradients. Often converges faster than `"gd"` at similar
#'       computational cost.}
#'     \item{`"nesterov"`}{Nesterov type accelerated gradient with a
#'       fixed momentum schedule. Can yield smoother convergence in
#'       early iterations, at the price of slightly more bookkeeping.}
#'     \item{`"bfgs"`}{Quasi Newton update in the latent space, with a
#'       limited memory BFGS type approximation of the Hessian. This
#'       gives the most accurate coefficients in a small number of
#'       iterations but requires more linear algebra per step.}
#'   }
#'
#'   The default is `"bb"`, which provides a good balance between speed
#'   and robustness in most examples.
#'
#' @param diag Logical. If TRUE, store iteration level diagnostics.
#' 
#' @details
#' The function first computes PLS components using the same procedure
#' as [big_pls_cox_fast()], then holds these components fixed and
#' optimises the Cox partial log-likelihood in the reduced space.
#'
#' The coefficients stored in `fit$coefficients` are the result of the
#' chosen optimisation method and are approximate. The field
#' `fit$cox_fit` contains the Cox model refitted by
#' [survival::coxph()] on the final PLS scores. Prediction functions
#' use the coefficients from `cox_fit` so that linear predictors are
#' directly interpretable as Cox risk scores.
#'
#' The optimisation tolerances control the compromise between speed
#' and accuracy. If you need very close agreement with the exact PLS
#' Cox solution, use [big_pls_cox_fast()]. If you only need accurate
#' risk rankings and want to fit many models, the gradient based
#' method with `"bb"` or `"bfgs"` is usually sufficient.
#'
#' @return 
#' 
#' An object of class `"big_pls_cox_gd"` that contains:
#' * `coefficients`: Estimated Cox regression coefficients on the latent scores.
#' * `loglik`: Final partial log-likelihood value.
#' * `iterations`: Number of gradient-descent iterations performed.
#' * `converged`: Logical flag indicating whether convergence was achieved.
#' * `scores`: Matrix of latent score vectors (one column per component).
#' * `loadings`: Matrix of loading vectors associated with each component.
#' * `weights`: Matrix of PLS weight vectors.
#' * `center`: Column means used to centre the predictors.
#' * `scale`: Column scales (standard deviations) used to standardise the predictors.
#' * `keepX`: Vector describing the number of predictors retained per component.
#' * `time`: Numeric vector of follow-up times.
#' * `status`: Numeric or integer vector containing the event indicators.
#' * `loglik_trace`: Trace of the loglik.
#' * `step_trace`: Trace of the steps
#' * `gradnorm_trace`: Trace of the gradnorm.
#' * `cox_fit`: Final Cox model fitted on the components.
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
#' library(bigmemory)
#' set.seed(1)
#' n <- 50
#' p <- 10
#' X <- bigmemory::as.big.matrix(matrix(rnorm(n * p), n, p))
#' time <- rexp(n, rate = 0.1)
#' status <- rbinom(n, 1, 0.7)
#' fit <- big_pls_cox_gd(X, time, status, ncomp = 3, max_iter = 200)
#' str(fit)
#' head(fit$scores)
#' @export
big_pls_cox_gd <- function(X, 
                           time, 
                           status, 
                           ncomp = 2L, 
                           max_iter = 2000L,
                           tol = 1e-8, 
                           learning_rate = 0.05, 
                           keepX = NULL,
                           coxfit = TRUE,
                           method = c("gd", "bb", "nesterov", "bfgs"),
                           diag = TRUE) {
  if (!requireNamespace("bigmemory", quietly = TRUE)) {
    stop("Package 'bigmemory' is required")
  }

  method <- match.arg(method)
  method_code <- match(method, c("gd", "bb", "nesterov", "bfgs")) - 1L
  
  if (is.matrix(X)) {
    X <- bigmemory::as.big.matrix(X)
  }
  if (!inherits(X, "big.matrix")) {
    stop("`X` must be a matrix or a big.matrix object", call. = FALSE)
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

  fit <- big_pls_cox_gd_cpp(
    X_ptr         = X@address,
    time          = as.numeric(time),
    status        = as.numeric(status),
    ncomp         = as.integer(ncomp),
    max_iter      = as.integer(max_iter),
    tol           = tol,
    learning_rate = learning_rate,
    keepX         = keep_vec,
    method_code   = method_code,
    return_diag   = diag
  )
  
  if (isTRUE(coxfit)) {
  # optional exact Cox refit in score space
  scores_df <- as.data.frame(fit$scores)
  colnames(scores_df) <- paste0("comp", seq_len(ncomp))
  scores_df$time   <- as.numeric(time)
  scores_df$status <- as.numeric(status)
  
  cox_fit <- survival::coxph(
    survival::Surv(time, status) ~ .,
    data  = scores_df,
    ties  = "efron",
    x     = FALSE
  )
  }
  
  fit$coefficients <- as.numeric(fit$coefficients)
  fit$center <- as.numeric(fit$center)
  fit$scale <- as.numeric(fit$scale)

  fit$keepX <- as.integer(keep_vec)
  fit$time <- as.numeric(time)
  fit$status <- as.numeric(status)

  fit$cox_fit <- cox_fit
  
  class(fit) <- "big_pls_cox_gd"
  fit
}
