#' Information criteria for component selection
#'
#' Computes log-likelihood, AIC and BIC values for nested models using the
#' latent components estimated by [big_pls_cox()] or [big_pls_cox_gd()].
#'
#' @param object A fitted object of class `big_pls_cox` or `big_pls_cox_gd`.
#' @param max_comp Maximum number of components to consider. Defaults to all
#'   components stored in the model.
#' 
#' @aliases component_information, component_information.big_pls_cox, 
#' component_information.big_pls_cox_gd, select_ncomp
#' 
#' @return A data frame with columns `ncomp`, `loglik`, `AIC`, and `BIC`.
#' @export
component_information <- function(object, max_comp = ncol(object$scores)) {
  UseMethod("component_information")
}

#' Component_information for [big_pls_cox()]
#' @rdname component_information
#' @export
component_information.big_pls_cox <- function(object, max_comp = ncol(object$scores)) {
  compute_component_information(object$scores, object$time, object$status, max_comp)
}

#' Component_information for [big_pls_cox_gd()]
#' @rdname component_information
#' @export
component_information.big_pls_cox_gd <- function(object, max_comp = ncol(object$scores)) {
  compute_component_information(object$scores, object$time, object$status, max_comp)
}

#' Select number of components
#' @rdname component_information
#' @param criterion Criterion to optimise: `"AIC"`, `"BIC"` or `"loglik"`.
#' @param ... Passed to [component_information()].
#' @return A list with the table of information criteria and the recommended
#'   number of components.
#' @export
select_ncomp <- function(object, criterion = c("AIC", "BIC", "loglik"), ...) {
  criterion <- match.arg(criterion)
  info <- component_information(object, ...)
  if (!nrow(info)) {
    stop("No component information available")
  }
  if (criterion == "loglik") {
    idx <- which.max(info$loglik)
  } else if (criterion == "AIC") {
    idx <- which.min(info$AIC)
  } else {
    idx <- which.min(info$BIC)
  }
  list(
    summary = info,
    criterion = criterion,
    opt_ncomp = info$ncomp[idx]
  )
}

compute_component_information <- function(scores, time, status, max_comp) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required to compute information criteria", call. = FALSE)
  }
  if (length(time) != nrow(scores) || length(status) != nrow(scores)) {
    stop("Length of `time` and `status` must match number of rows in `scores`")
  }
  max_comp <- min(max_comp, ncol(scores))
  if (max_comp < 1L) {
    stop("`max_comp` must be at least 1")
  }
  n <- length(time)
  y <- cbind(time, status)
  strata <- rep(1, n)
  offset <- rep(0, n)
  wt <- rep(1, n)
  rownms <- as.character(seq_len(n))
  control <- survival::coxph.control(iter.max = 50)
  
  res <- data.frame(
    ncomp = seq_len(max_comp),
    loglik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_
  )
  
  for (k in seq_len(max_comp)) {
    xk <- scores[, seq_len(k), drop = FALSE]
    fit <- tryCatch(
      survival::coxph.fit(
        x = xk,
        y = y,
        strata = strata,
        offset = offset,
        init = rep(0, k),
        control = control,
        weights = wt,
        method = "efron",
        rownames = rownms
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    loglik <- fit$loglik[2]
    res$loglik[k] <- loglik
    res$AIC[k] <- -2 * loglik + 2 * k
    res$BIC[k] <- -2 * loglik + log(n) * k
  }
  res
}
