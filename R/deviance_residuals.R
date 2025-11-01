#' Cox deviance residuals via C++ backends
#'
#' @description Compute martingale and deviance residuals for Cox models without
#' materialising intermediate survival fits in R. The functions rely on
#' dedicated C++ implementations that operate either on in-memory vectors or on
#' [`bigmemory::big.matrix`] objects to enable streaming computations on large
#' datasets.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric or integer vector of the same length as `time` giving
#'   the event indicators (1 for an event, 0 for censoring).
#' @param weights Optional non-negative case weights. When supplied they must
#'   have the same length as `time`.
#' @param X A [`bigmemory::big.matrix`][bigmemory::big.matrix-class] storing the
#'   survival information column-wise.
#' @param time_col,status_col Integer indices pointing to the columns of `X`
#'   that contain the follow-up time and event indicator respectively.
#' @param coef Numeric vector of regression coefficients used to evaluate the
#'   partial log-likelihood and deviance on a `big.matrix` design.
#' @param iterations Number of iterations used by [`bench::mark`] when
#'   benchmarking the residual computations.
#' @param methods Optional named list of alternative residual implementations to
#'   compare against in [`benchmark_deviance_residuals`].
# #' @param ... Additional arguments passed to lower-level routines.
#'
#' @details
#' * [`cox_deviance_residuals()`] operates on standard R vectors and matches the
#'   output of `residuals(coxph(...), type = "deviance")` for right-censored
#'   data without ties.
#' * [`cox_deviance_residuals_big()`] keeps the computation in C++ while reading
#'   directly from a `big.matrix`, avoiding extra copies.
#' * [`cox_partial_deviance_big()`] evaluates the partial log-likelihood and
#'   deviance for a given coefficient vector and big design matrix. This is
#'   useful when selecting the number of latent components via information
#'   criteria.
#'
#' [`benchmark_deviance_residuals()`] compares the dedicated C++ implementation
#' against reference approaches (for example, the `survival` package) using
#' [`bench::mark`]. The function returns a tibble with iteration statistics.
#'
#' @return
#' * [`cox_deviance_residuals()`] and [`cox_deviance_residuals_big()`] return a
#'   numeric vector of deviance residuals.
#' * [`cox_deviance_details()`] returns a list with cumulative hazard,
#'   martingale, and deviance residuals.
#' * [`cox_partial_deviance_big()`] returns a list containing the partial
#'   log-likelihood, deviance, and the evaluated linear predictor.
#' * [`benchmark_deviance_residuals()`] returns a [`tibble::tibble`].
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(123)
#'   time <- rexp(50)
#'   status <- rbinom(50, 1, 0.6)
#'   dr_cpp <- cox_deviance_residuals(time, status)
#'   dr_surv <- residuals(survival::coxph(survival::Surv(time, status) ~ 1),
#'                        type = "deviance")
#'   all.equal(unname(dr_cpp), unname(dr_surv), tolerance = 1e-6)
#' }
#'
#' @name cox_deviance_residuals
#' @rdname cox_deviance_residuals
#' @export
cox_deviance_residuals <- function(time, status, weights = NULL) {
  if (length(time) != length(status)) {
    stop("`time` and `status` must have the same length", call. = FALSE)
  }
  if (!is.null(weights) && length(weights) != length(time)) {
    stop("`weights` must have the same length as `time`", call. = FALSE)
  }
  cox_deviance_residuals_cpp(as.numeric(time), as.numeric(status),
                             if (is.null(weights)) NULL else as.numeric(weights))
}

#' @rdname cox_deviance_residuals
#' @export
cox_deviance_details <- function(time, status, weights = NULL) {
  if (length(time) != length(status)) {
    stop("`time` and `status` must have the same length", call. = FALSE)
  }
  if (!is.null(weights) && length(weights) != length(time)) {
    stop("`weights` must have the same length as `time`", call. = FALSE)
  }
  cox_deviance_details_cpp(as.numeric(time), as.numeric(status),
                           if (is.null(weights)) NULL else as.numeric(weights))
}

#' @rdname cox_deviance_residuals
#' @export
cox_deviance_residuals_big <- function(X, time_col, status_col, weights = NULL) {
  if (!inherits(X, "big.matrix")) {
    stop("`X` must be a big.matrix", call. = FALSE)
  }
  if (!is.null(weights) && length(weights) != nrow(X)) {
    stop("`weights` must have the same length as the number of rows in `X`", call. = FALSE)
  }
  cox_deviance_residuals_big_cpp(X@address, as.integer(time_col), as.integer(status_col),
                                 if (is.null(weights)) NULL else as.numeric(weights))
}

#' @rdname cox_deviance_residuals
#' @export
cox_partial_deviance_big <- function(X, coef, time, status) {
  if (!inherits(X, "big.matrix")) {
    stop("`X` must be a big.matrix", call. = FALSE)
  }
  if (length(coef) != ncol(X)) {
    stop("`coef` must have length equal to the number of columns in `X`", call. = FALSE)
  }
  cox_partial_deviance_big_cpp(X@address, as.numeric(coef), as.numeric(time),
                               as.numeric(status))
}

#' @rdname cox_deviance_residuals
#' @export
benchmark_deviance_residuals <- function(time, status, iterations = 25,
                                         methods = list()) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("Package 'bench' is required for benchmarking", call. = FALSE)
  }
  base_methods <- list(
    cpp = quote(cox_deviance_residuals(time, status))
  )
  if (requireNamespace("survival", quietly = TRUE)) {
    base_methods$survival <- quote(residuals(survival::coxph(
      survival::Surv(time, status) ~ 1
    ), type = "deviance"))
  }
  eval_env <- list2env(list(
    time = as.numeric(time),
    status = as.numeric(status)
  ), parent = parent.frame())
  all_methods <- c(base_methods, methods)
  do.call(bench::mark, c(all_methods, list(iterations = iterations, check = TRUE)),
          envir = eval_env)
}