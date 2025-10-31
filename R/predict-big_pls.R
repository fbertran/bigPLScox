#' Predict method for big-memory PLS-Cox models
#'
#' @param object A model fitted with [big_pls_cox()].
#' @param newdata Optional matrix, data frame or [`bigmemory::big.matrix`]
#'   containing predictors to project on the latent space. When `NULL` the
#'   training scores are used.
#' @param type Type of prediction: `"link"` for the linear predictor, `"risk"`
#'   or `"response"` for the exponential of the linear predictor, or
#'   `"components"` to obtain latent scores.
#' @param comps Integer vector indicating which components to use. Defaults to
#'   all available components.
#' @param coef Optional coefficient vector overriding the fitted Cox model
#'   coefficients.
#' @param ... Unused.
#'
#' @return Depending on `type`, either a numeric vector of predictions or a
#'   matrix of component scores.
#' @export
predict.big_pls_cox <- function(object, newdata = NULL,
                                type = c("link", "risk", "response", "components"),
                                comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  total_comp <- ncol(object$scores)
  if (is.null(comps)) {
    comps <- seq_len(total_comp)
  }
  comps <- as.integer(comps)
  if (length(comps) == 0L) {
    stop("`comps` must contain at least one component")
  }
  if (any(comps < 1L) || any(comps > total_comp)) {
    stop("`comps` indices are out of bounds")
  }
  
  scores <- compute_big_pls_scores(object, newdata, comps)
  
  if (type == "components") {
    return(scores)
  }
  
  if (is.null(coef)) {
    coef <- object$cox_fit$coefficients
    if (is.null(coef)) {
      stop("Cox model coefficients are not available; provide them explicitly via `coef`")
    }
  }
  coef <- as.numeric(coef)
  if (length(coef) < max(comps)) {
    stop("`coef` must have at least max(comps) entries")
  }
  eta <- as.numeric(scores %*% coef[comps])
  switch(type,
         link = eta,
         risk = exp(eta),
         response = exp(eta))
}

compute_big_pls_scores <- function(object, newdata, comps) {
  comps <- as.integer(comps)
  if (is.null(newdata)) {
    return(object$scores[, comps, drop = FALSE])
  }
  if (inherits(newdata, "big.matrix")) {
    return(big_pls_cox_transform_cpp(newdata@address, object$center, object$scale,
                                     object$weights, object$loadings, comps))
  }
  if (is.data.frame(newdata)) {
    newdata <- as.matrix(newdata)
  }
  if (!is.matrix(newdata)) {
    stop("`newdata` must be a matrix, data frame, or big.matrix")
  }
  storage.mode(newdata) <- "double"
  matrix_pls_cox_transform_cpp(newdata, object$center, object$scale,
                               object$weights, object$loadings, comps)
}

#' @rdname predict.big_pls_cox
#' @export
predict.big_pls_cox_gd <- function(object, newdata = NULL,
                                   type = c("link", "risk", "response", "components"),
                                   comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  total_comp <- ncol(object$scores)
  if (is.null(comps)) {
    comps <- seq_len(total_comp)
  }
  comps <- as.integer(comps)
  if (length(comps) == 0L) {
    stop("`comps` must contain at least one component")
  }
  if (any(comps < 1L) || any(comps > total_comp)) {
    stop("`comps` indices are out of bounds")
  }
  
  scores <- compute_big_pls_scores(object, newdata, comps)
  if (type == "components") {
    return(scores)
  }
  
  if (is.null(coef)) {
    coef <- object$coefficients
    if (is.null(coef)) {
      stop("Coefficients are not stored in the gradient-descent fit; supply them via `coef`")
    }
  }
  coef <- as.numeric(coef)
  if (length(coef) < max(comps)) {
    stop("`coef` must have at least max(comps) entries")
  }
  eta <- as.numeric(scores %*% coef[comps])
  switch(type,
         link = eta,
         risk = exp(eta),
         response = exp(eta))
}
