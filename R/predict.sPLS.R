#' @exportS3Method stats::predict
predict.sPLS <- function (object, newdata, scale.X=TRUE, scale.Y=TRUE, ...) 
{
  if (missing(newdata)) 
    stop("No new data available.")
  X = object$X
  Y = object$Y
  q = ncol(Y)
  p = ncol(X)
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
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  means.X = attr(X, "scaled:center")
  means.Y = attr(Y, "scaled:center")
  sigma.X = attr(X, "scaled:scale")
  sigma.Y = attr(Y, "scaled:scale")
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
    Y.temp2 = scale(Y.temp, center = FALSE, scale = ifelse(scale.Y,1/sigma.Y,FALSE))
    Y.temp3 = scale(Y.temp2, center = ifelse(scale.Y,-means.Y,FALSE), scale = FALSE)
    Y.hat[, , h] = Y.temp3
    t.pred[, h] = scale(newdata, center = ifelse(scale.X,means.X,FALSE), scale = ifelse(scale.X,sigma.X,FALSE)) %*% 
      W[, h]
    B.hat[, , h] = B
  }
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(Y.hat) = rownames(newdata)
  colnames(Y.hat) = colnames(Y)
  return(invisible(list(predict = Y.hat, variates = t.pred, 
                        B.hat = B.hat, betay = betay)))
}
