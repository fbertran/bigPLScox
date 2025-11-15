#' Summary for big_pls_cox objects
#'
#' @param object A \code{big_pls_cox} object.
#' @param digits Number of digits to print.
#' @param ... Unused, for S3 compatibility.
#'
#' @export
summary.big_pls_cox <- function(object, digits = 3, ...) {
  if (!inherits(object, "big_pls_cox")) {
    stop("summary.big_pls_cox() requires an object of class 'big_pls_cox'.")
  }
  
  n <- nrow(object$scores)
  p <- nrow(object$loadings)
  ncomp <- ncol(object$scores)
  
  cat("Partial Least Squares Cox model (big_pls_cox)\n")
  cat("-------------------------------------------------\n")
  cat("Observations:", n, "\n")
  cat("Predictors  :", p, "\n")
  cat("Components  :", ncomp, "\n\n")
  
  # variance of scores (should be ~1 with your current scaling)
  var_scores <- apply(object$scores, 2, var)
  cat("Score variances (per component):\n")
  print(round(var_scores, digits = digits))
  cat("\n")
  
  if (!is.null(object$keepX)) {
    cat("keepX (sparsity per component):", paste(object$keepX, collapse = ", "), "\n\n")
  }
  
  # Cox fit summary if available
  if (!is.null(object$cox_fit)) {
    cf <- object$cox_fit
    cat("Cox model on PLS scores:\n")
    print(stats::coef(summary(cf)), digits = digits)
    
    # simple concordance estimate
    eta <- stats::predict(cf, type = "link")
    cindex <- try({
      survival::concordance(survival::Surv(object$time, object$status) ~ eta)$concordance
    }, silent = TRUE)
    
    if (!inherits(cindex, "try-error")) {
      cat("\nHarrell's C-index:", round(cindex, digits = digits), "\n")
    }
  }
  
  invisible(object)
}

#' Print method for big_pls_cox objects
#'
#' @param x An object of class \code{big_pls_cox}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @method print big_pls_cox
#' @export
print.big_pls_cox <- function(x, ...) {
  cat("PLS-Cox model with", ncol(x$scores), "components\n")
  cat("Final Cox model coefficients:\n")
  print(x$cox_fit$coefficients)
  invisible(x)
}

#' Summary for big_pls_cox_gd objects
#'
#' @param object A \code{big_pls_cox_gd} object.
#' @param digits Number of digits to print.
#' @param ... Unused, for S3 compatibility.
#'
#' @export
summary.big_pls_cox_gd <- function(object, digits = 3, ...) {
  if (!inherits(object, "big_pls_cox_gd")) {
    stop("summary.big_pls_cox_gd() requires an object of class 'big_pls_cox_gd'.")
  }
  
  n <- nrow(object$scores)
  p <- nrow(object$loadings)
  ncomp <- ncol(object$scores)
  
  cat("Partial Least Squares Cox model (big_pls_cox_gd)\n")
  cat("-------------------------------------------------\n")
  cat("Observations:", n, "\n")
  cat("Predictors  :", p, "\n")
  cat("Components  :", ncomp, "\n\n")
  
  var_scores <- apply(object$scores, 2, var)
  cat("Score variances (per component):\n")
  print(round(var_scores, digits = digits))
  cat("\n")
  
  if (!is.null(object$keepX)) {
    cat("keepX (sparsity per component):", paste(object$keepX, collapse = ", "), "\n\n")
  }
  
  if (!is.null(object$coefficients)) {
    cat("Coefficients (gradient-descent fit):\n")
    print(round(object$coefficients, digits = digits))
    cat("\n")
  }
  
  if (!is.null(object$loglik) && is.finite(object$loglik)) {
    cat("Partial log-likelihood:", round(object$loglik, digits = digits), "\n")
  }
  if (!is.null(object$iterations)) {
    cat("Iterations           :", object$iterations, "\n")
  }
  if (!is.null(object$converged)) {
    cat("Converged            :", object$converged, "\n")
  }
  cat("\n")
  
  # if center/scale present, briefly report
  if (!is.null(object$center) && !is.null(object$scale)) {
    cat("Predictors centered and scaled.\n")
  }
  
  # quick C-index using GD coefficients (if finite)
  if (!any(is.na(object$coefficients)) && is.finite(object$loglik)) {
    eta <- drop(object$scores %*% object$coefficients)
    cindex <- try({
      survival::concordance(survival::Surv(object$time, object$status) ~ eta)$concordance
    }, silent = TRUE)
    if (!inherits(cindex, "try-error")) {
      cat("Harrell's C-index    :", round(cindex, digits = digits), "\n")
    }
  }
  
  invisible(object)
}

#' Print method for big_pls_cox_gd objects
#'
#' @param x An object of class \code{big_pls_cox_gd}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @method print big_pls_cox_gd
#' @export
print.big_pls_cox_gd <- function(x, ...) {
  cat("PLS-Cox model with", ncol(x$scores), "components\n")
  cat("Final Cox model coefficients:\n")
  print(x$cox_fit$coefficients)
  invisible(x)
}

#' @noRd
print.summary.big_pls_cox_gd <- function(x, ...) {
  cat("bigPLScox GD model\n")
  cat("  Observations (n):", x$n, "\n")
  cat("  Predictors   (p):", x$p, "\n")
  cat("  Components     :", x$ncomp, "\n")
  if (!is.null(x$keepX)) {
    cat("  keepX per component:", paste(x$keepX, collapse = ", "), "\n")
  }
  cat("  Iterations      :", x$iterations, "\n")
  cat("  Converged       :", x$converged, "\n")
  cat("  Log-likelihood  :", x$loglik, "\n\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  invisible(x)
}
