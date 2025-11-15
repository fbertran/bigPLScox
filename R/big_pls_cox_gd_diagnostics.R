#' Extract Diagnostics from a big_pls_cox_gd Model
#'
#' @param object A model returned by \code{big_pls_cox_gd()}.
#' @return A list with log-likelihood, step sizes, gradient norms.
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
#' gd_diagnostics(fit)
#' @export
gd_diagnostics <- function(object) {
  stopifnot(inherits(object, "big_pls_cox_gd"))
  
  list(
    iterations      = seq_along(object$loglik_trace),
    loglik          = object$loglik_trace,
    step_sizes      = object$step_trace,
    gradient_norm   = object$gradnorm_trace,
    coef_trace = object$coef_trace,
    eta_trace = object$eta_trace
  )
}
