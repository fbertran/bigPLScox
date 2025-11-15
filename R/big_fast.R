#' Partial Least Squares Components for Cox Models (fast backend)
#'
#' @description
#' Compute PLS components for Cox models, using a fast C++ backend for both
#' in-memory matrices and \code{bigmemory::big.matrix} objects.
#'
#' @inherit big_pls_cox params return references
#' @export
big_pls_cox_fast <- function(X, time, status, ncomp = 2L,
                             control = survival::coxph.control(),
                             keepX = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required", call. = FALSE)
  }
  
  if (!is.numeric(time) || !is.numeric(status)) {
    stop("time and status must be numeric vectors")
  }
  if (length(time) != length(status)) {
    stop("time and status must have the same length")
  }
  n <- length(time)
  
  is_big <- inherits(X, "big.matrix")
  
  if (!is.matrix(X) && !is_big) {
    stop("X must be either a numeric matrix or a big.matrix")
  }
  if (is.matrix(X)) {
    if (nrow(X) != n) {
      stop("Number of rows in X must match length of time/status")
    }
  } else {
    if (!requireNamespace("bigmemory", quietly = TRUE)) {
      stop("Package 'bigmemory' is required for big.matrix input",
           call. = FALSE)
    }
    if (nrow(X) != n) {
      stop("Number of rows in X must match length of time/status")
    }
  }
  
  ncomp <- as.integer(ncomp)
  if (length(ncomp) != 1L || ncomp <= 0L) {
    stop("ncomp must be a positive integer")
  }
  
  p <- if (is_big) ncol(X) else ncol(X)
  if (ncomp > min(n, p)) {
    stop("ncomp must be <= min(nrow(X), ncol(X))")
  }
  
  ## keepX handling
  if (is.null(keepX)) {
    keepX <- rep.int(0L, ncomp)
  } else {
    keepX <- as.integer(keepX)
    if (length(keepX) == 1L) {
      keepX <- rep.int(keepX, ncomp)
    }
    if (length(keepX) != ncomp) {
      stop("keepX must be NULL, a single integer, or a vector of length ncomp")
    }
    if (any(is.na(keepX)) || any(keepX < 0L) || any(keepX > p)) {
      stop("keepX entries must lie between 0 and ncol(X)")
    }
  }
  
  ## Control is only used for final Cox fit
  if (!inherits(control, "coxph.control")) {
    control <- do.call(survival::coxph.control, control)
  }
  
  ## Column stats
  if (is_big) {
    address <- methods::slot(X, "address")
    stats <- big_pls_cox_col_stats_cpp(address)
    means <- stats$mean
    sds   <- stats$sd
    res <- big_pls_cox_fast_big_cpp(
      xpMat  = address,
      time   = as.numeric(time),
      status = as.numeric(status),
      ncomp  = ncomp,
      means  = means,
      sds    = sds,
      keepX  = as.integer(keepX)
    )
  } else {
    X <- as.matrix(X)
    means <- colMeans(X)
    sds <- apply(X, 2L, stats::sd)
    bad <- !is.finite(sds) | sds == 0
    if (any(bad)) sds[bad] <- 1
    
    res <- big_pls_cox_fast_dense_cpp(
      X      = X,
      time   = as.numeric(time),
      status = as.numeric(status),
      ncomp  = ncomp,
      means  = means,
      sds    = sds,
      keepX  = as.integer(keepX)
    )
  }
  
  scores   <- res$scores
  loadings <- res$loadings
  weights  <- res$weights
  
  ## Final Cox fit on components (scores)
  scores_df <- as.data.frame(scores)
  colnames(scores_df) <- paste0("comp", seq_len(ncol(scores)))
  cox_fit <- survival::coxph(
    survival::Surv(time, status) ~ .,
    data = scores_df,
    ties = "efron",
    x = FALSE
  )
  
  structure(
    list(
      scores  = scores,
      loadings = loadings,
      weights = weights,
      center  = as.numeric(means),
      scale   = as.numeric(sds),
      cox_fit = cox_fit,
      keepX   = as.integer(keepX),
      time    = as.numeric(time),
      status  = as.numeric(status)
    ),
    class = "big_pls_cox_fast"
  )
}

#' Summary for big_pls_cox objects
#'
#' @param object A \code{big_pls_cox} object.
#' @param ... Not used.
#'
#' @export
summary.big_pls_cox_fast <- function(object, ...) {
  n     <- nrow(object$scores)
  p     <- nrow(object$loadings)
  ncomp <- ncol(object$scores)
  
  s_cox <- summary(object$cox_fit)
  
  res <- list(
    n      = n,
    p      = p,
    ncomp  = ncomp,
    keepX  = object$keepX,
    center = object$center,
    scale  = object$scale,
    cox    = s_cox
  )
  class(res) <- "summary.big_pls_cox"
  res
}

#' Print method for summary.big_pls_cox_fast objects
#'
#' @param x An object of class \code{summary.big_pls_cox_fast}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @method print summary.big_pls_cox_fast
#' @export
print.summary.big_pls_cox_fast <- function(x, ...) {
  cat("bigPLScox model\n")
  cat("  Observations (n):", x$n, "\n")
  cat("  Predictors   (p):", x$p, "\n")
  cat("  Components     :", x$ncomp, "\n")
  if (!is.null(x$keepX)) {
    cat("  keepX per component:", paste(x$keepX, collapse = ", "), "\n")
  }
  cat("\nCox model summary:\n\n")
  print(x$cox)
  invisible(x)
}


#' Predictions for fast big PLS–Cox fits
#'
#' @param object An object of class \code{"big_pls_cox_fast"}.
#' @param newdata Optional matrix or \code{big.matrix} with predictors.
#' @param type Type of prediction: \code{"link"}, \code{"risk"},
#'   \code{"response"}, or \code{"components"}.
#' @param comps Optional components to use. Defaults to all.
#' @param coef Optional coefficient vector. Defaults to \code{object$coefficients}.
#' @param ... Ignored.
#'
#' @export
predict.big_pls_cox_fast <- function(object,
                                     newdata = NULL,
                                     type = c("link", "risk", "response", "components"),
                                     comps = NULL,
                                     coef = NULL,
                                     ...) {
  type <- match.arg(type)
  
  if (is.null(comps)) {
    comps <- seq_len(ncol(object$scores))
  }
  
  if (is.null(coef)) {
    coef <- object$cox_fit$coefficients
  }
  coef <- as.numeric(coef)[comps]
  
  scores_new <- if (is.null(newdata)) {
      scores_new <- object$scores[, comps, drop = FALSE]
  } else { 
    if (inherits(newdata, "big.matrix")) {
      address <- methods::slot(newdata, "address")
      big_pls_cox_transform_big_cpp(
        address,          # slot(X_new, "address")
        means  = object$center,
        sds    = object$scale,
        weights = object$weights,
        loadings = object$loadings,
        comps  = comps
      )
    } else { if (inherits(newdata, "matrix")) {
      big_pls_cox_transform_dense_cpp(
        newdata,
        means   = object$center,
        sds     = object$scale,
        weights = object$weights,
        loadings = object$loadings,
        comps   = comps
      )
    } else {
      stop("`newdata` must be NULL, a numeric matrix, or a big.matrix.")
    }
    }
  }  
  scores_new <- as.matrix(scores_new)
  
  if (type == "components") {
    return(scores_new)
  }
  
  eta <- drop(scores_new %*% coef)
  
  switch(
    type,
    link     = eta,
    risk     = exp(eta),
    response = exp(eta)  # for Cox, 'risk' and 'response' are usually identical
  )
}
