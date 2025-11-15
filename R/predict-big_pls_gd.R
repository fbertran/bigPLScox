#' @export
predict.big_pls_cox_gd <- function(object, newdata, type = c("link", "risk", "components"), ...) {
  type <- match.arg(type)
  
  if (missing(newdata))
    stop("newdata must be provided")
  
  if (!inherits(newdata, "matrix") && !inherits(newdata, "big.matrix"))
    stop("newdata must be a matrix or big.matrix")
  
  center <- object$center
  scale  <- object$scale
  weights <- object$weights
  loadings <- object$loadings
  
  if (inherits(newdata, "big.matrix")) {
    Xmat <- as.matrix(newdata)
  } else {
    Xmat <- newdata
  }
  
  # Standardize
  Xs <- sweep(Xmat, 2, center, "-")
  Xs <- sweep(Xs, 2, scale, "/")
  
  # Compute components EXACTLY as in training
  Tnew <- matrix_pls_cox_transform_cpp(
    Xs, center, scale,
    weights, loadings,
    comps = seq_len(ncol(weights))
  )
  
  if (type == "components")
    return(Tnew)
  
  lp <- drop(Tnew %*% object$coefficients)
  
  switch(type,
         link = lp,
         risk = exp(lp))
}
