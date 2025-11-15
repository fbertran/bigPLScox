#' Predict method for big_pls_cox_gd
#'
#' @param object A big_pls_cox_gd model.
#' @param newdata A numeric matrix or big.matrix. If NULL, uses training scores.
#' @param type One of "link", "risk", "response", "components".
#' @param comps Integer vector of components to use. Default: all.
#' @param coef Optional coefficient vector. If NULL, uses coef(object$cox_fit).
#' @param ... Unused.
#'
#' @export
predict.big_pls_cox_gd <- function(object,
                                   newdata = NULL,
                                   type = c("link", "risk", "response", "components"),
                                   comps = NULL,
                                   coef = NULL,
                                   ...) {
  
  type <- match.arg(type)
  
  total_comp <- ncol(object$scores)
  
  # ---- components to use ----
  if (is.null(comps)) {
    comps <- seq_len(total_comp)
  }
  comps <- as.integer(comps)
  if (any(comps < 1L) || any(comps > total_comp)) {
    stop("Component indices out of bounds")
  }
  
  # ---- compute component scores ----
  if (is.null(newdata)) {
    scores <- object$scores[, comps, drop = FALSE]
  } else {
    scores <- big_pls_cox_transform(
      X        = newdata,
      means    = object$center,
      sds      = object$scale,
      weights  = object$weights,
      loadings = object$loadings,
      comps    = comps
    )
  }
  
  # return only scores if requested
  if (type == "components") {
    return(scores)
  }
  
  # ---- coefficients in Cox space ----
  if (is.null(coef)) {
    if (!is.null(object$cox_fit)) {
      coef <- stats::coef(object$cox_fit)
    } else if (!is.null(object$coefficients)) {
      coef <- object$coefficients
      warning("Using internal GD coefficients; consider providing `coef`.")
    } else {
      stop("No coefficients available in object; supply `coef` manually.")
    }
  }
  
  coef <- as.numeric(coef)
  if (length(coef) < max(comps)) {
    stop("Supplied coef has fewer entries than required by `comps`.")
  }
  
  # ---- linear predictor ----
  eta <- as.numeric(scores %*% coef[comps])
  
  switch(type,
         link     = eta,
         risk     = exp(eta),
         response = exp(eta))
}