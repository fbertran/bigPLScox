#' Predict responses and latent scores from PLS fits
#'
#' These prediction helpers reconstruct the response matrix and latent
#' component scores for partial least squares (PLS) models fitted inside the
#' Cox-PLS toolbox. They support group PLS, sparse PLS, sparse-group PLS, and
#' classical PLS models created by [sgPLS::gPLS()], [sgPLS::sPLS()],
#' [sgPLS::sgPLS()], or [pls.cox()].
#'
#' @param object A fitted PLS model returned by [sgPLS::gPLS()],
#'   [sgPLS::sPLS()], [sgPLS::sgPLS()], or [pls.cox()].
#' @param newdata Numeric matrix or data frame with the same number of columns
#'   as the training design matrix used when fitting `object`.
#' @param scale.X,scale.Y Logical flags indicating whether the predictors and
#'   responses supplied in `newdata` should be centred and scaled according to
#'   the training statistics stored in `object`.
#' @param ... Unused arguments included for compatibility with the generic
#'   [stats::predict()] signature.
#'
#' @return A list containing reconstructed responses, latent component scores,
#'   and regression coefficients. The exact elements depend on the specific PLS
#'   algorithm but always include components named `predict`, `variates`, and
#'   `B.hat`.
#'
#' @seealso [coxgpls()], [coxsgpls()], [coxspls_sgpls()], and
#'   [coxDKgplsDR()] for Cox model wrappers that return PLS fits using these
#'   prediction methods.
#'
#' @references
#'   Bastien, P., Bertrand, F., Meyer, N., & Maumy-Bertrand, M. (2015).
#'   Deviance residuals-based sparse PLS and sparse kernel PLS for censored
#'   data. *Bioinformatics*, 31(3), 397–404. <doi:10.1093/bioinformatics/btu660>
#'
#' @examples
#'
#' n <- 100
#' sigma.gamma <- 1
#' sigma.e <- 1.5
#' p <- 400
#' q <- 500
#' theta.x1 <- c(rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5), rep(1.5,15), 
#'               rep(0, 5), rep(-1.5, 15), rep(0, 325))
#' theta.x2 <- c(rep(0, 320), rep(1, 15), rep(0, 5), rep(-1, 15), rep(0, 5),
#'               rep(1.5, 15), rep(0, 5), rep(-1.5, 15), rep(0, 5))
#' theta.y1 <- 1
#' theta.y2 <- 1
#'
#' Sigmax <- matrix(0, nrow = p, ncol = p)
#' diag(Sigmax) <- sigma.e ^ 2
#' Sigmay <- matrix(0,nrow = 1, ncol = 1)
#' diag(Sigmay) <- sigma.e ^ 2
#'
#' set.seed(125)
#'
#' gam1 <- rnorm(n)
#' gam2 <- rnorm(n)
#'
#' X <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.x1, theta.x2),
#'  nrow = 2, byrow = TRUE) + mvtnorm::rmvnorm(n, mean = rep(0, p), sigma =
#'  Sigmax, method = "svd")
#' Y <- matrix(c(gam1, gam2), ncol = 2, byrow = FALSE) %*% matrix(c(theta.y1, theta.y2), 
#' nrow = 2, byrow = TRUE) + rnorm(n,0,sd=sigma.e)
#'
#' ind.block.x <- seq(20, 380, 20)
#'
#' model.gPLS <- sgPLS::gPLS(X, Y, ncomp = 2, mode = "regression", keepX = c(4, 4), 
#'                    keepY = c(4, 4), ind.block.x = ind.block.x)
#' head(predict(model.gPLS, newdata = X)$variates)
#'
#' @name predict_pls_latent
#' @rdname predict_pls_latent
#' @exportS3Method stats::predict
predict.gPLS=function (object, newdata, scale.X=TRUE, scale.Y=TRUE, ...)
{
  if (missing(newdata))
    stop("No new data available.")
  
  X <- object$X
  if (!is.null(X)) {
    X <- as.matrix(X)
  }
  Y <- object$Y
  if (!is.null(Y)) {
    Y <- as.matrix(Y)
  }
  
  a <- object$loadings$X
  b <- object$loadings$Y
  if (is.null(a)) {
    stop("The gPLS fit does not contain X loadings; cannot proceed with prediction.")
  }
  if (is.null(b)) {
    stop("The gPLS fit does not contain Y loadings; cannot proceed with prediction.")
  }
  
  if (is.null(Y)) {
    stop("The gPLS fit does not contain the response used during training; cannot compute predictions.")
  }
  
  p <- if (!is.null(X)) {
    ncol(X)
  } else {
    nrow(as.matrix(a))
  }
  q <- if (!is.null(Y)) {
    ncol(Y)
  } else {
    nrow(as.matrix(b))
  }
  
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ",
           p, " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ",
           p, " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p)
  }
  
  means.X <- if (!is.null(X)) attr(X, "scaled:center") else NULL
  if (is.null(means.X)) {
    means.X <- rep(0, p)
  }
  sigma.X <- if (!is.null(X)) attr(X, "scaled:scale") else NULL
  if (is.null(sigma.X)) {
    sigma.X <- rep(1, p)
  }
  means.Y <- if (!is.null(Y)) attr(Y, "scaled:center") else NULL
  if (is.null(means.Y)) {
    means.Y <- rep(0, q)
  }
  sigma.Y <- if (!is.null(Y)) attr(Y, "scaled:scale") else NULL
  if (is.null(sigma.Y)) {
    sigma.Y <- rep(1, q)
  }
  
  ncomp = object$ncomp
  c = object$mat.c
  newdata = as.matrix(newdata)
  ones = matrix(rep(1, nrow(newdata)), ncol = 1)
  B.hat = array(0, dim = c(p, q, ncomp))
  Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  Y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  variates.X = object$variates$X
  betay = list()
  for (h in 1:ncomp) {
    dd = coefficients(lm(Y ~ variates.X[, 1:h, drop = FALSE]))
    if (q == 1) {
      betay[[h]] = (dd[-1])
    }
    if (q >= 2) {
      betay[[h]] = (dd[-1, ])
    }
    W = a[, 1:h, drop = FALSE] %*% solve(t(c[, 1:h, drop = FALSE]) %*%
                                           a[, 1:h, drop = FALSE])
    B = W %*% drop(betay[[h]])
    Y.temp = scale(newdata, center = ifelse(scale.X,means.X,FALSE), scale = ifelse(scale.X,sigma.X,FALSE)) %*%
      as.matrix(B)
    Y.temp2 = scale(Y.temp, center = FALSE, scale = if(scale.Y){1/sigma.Y} else FALSE)
    Y.temp3 = scale(Y.temp2, center = if(scale.Y){-means.Y} else FALSE, scale = FALSE)
    Y.hat[, , h] = Y.temp3
    t.pred[, h] = scale(newdata, center = if(scale.X){means.X} else FALSE, scale = if(scale.X){sigma.X} else FALSE) %*%
      W[, h]
    B.hat[, , h] = B
  }
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(Y.hat) = rownames(newdata)
  colnames(Y.hat) = if (!is.null(Y)) colnames(Y) else rownames(as.matrix(b))
  return(invisible(list(predict = Y.hat, variates = t.pred,
                        B.hat = B.hat, betay = betay)))
}
