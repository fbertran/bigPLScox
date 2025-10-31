#' @title Internal bigPLScox functions
#' 
#' @name internal-bigPLScox
#' 
#' @description These are not to be called by the user.
#' 
#' @aliases ust spls.dv correctp.cox
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @references plsRcox, Cox-Models in a high dimensional setting in R, Frederic
#' Bertrand, Philippe Bastien, Nicolas Meyer and Myriam Maumy-Bertrand (2014).
#' Proceedings of User2014!, Los Angeles, page 152.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords internal
NULL

logplik = function (x, time, status, b, method = c("breslow", "efron"), return.all = FALSE) {
  method <- match.arg(method)
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(b)) {
    b <- as.matrix(b)
  }
  res <- cox_partial_loglik_cpp(x, as.numeric(time), as.numeric(status), b, method, return.all)
  res
}

getmin2 = function (lambda, cvm, cvsd) {
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idminl = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idminl]
  idmin2 = cvm >= semin
  #    if(lambda.min==-Inf){
  #    lambda.1se = -Inf
  #    } else {
  idmin2[idminl:length(idmin2)] = FALSE 
  lambda.1se = min(c(lambda[idmin2],min(lambda)), na.rm = TRUE)
  #    }
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}


correctp.cox=function (x, y, eta, K, kappa, select, fit, verbose=FALSE)  {
  force(K)
  if (min(eta) < 0 | max(eta) >= 1) {
    if (max(eta) == 1) {
      stop("eta should be strictly less than 1!")
    }
    if (length(eta) == 1) {
      stop("eta should be between 0 and 1!")
    }
    else {
      stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
    }
  }
  if (max(K) > ncol(x)) {
    stop("K cannot exceed the number of predictors! Pick up smaller K!")
  }
  if (max(K) >= nrow(x)) {
    stop("K cannot exceed the sample size! Pick up smaller K!")
  }
  if (min(K) <= 0 | !all(K%%1 == 0)) {
    if (length(K) == 1) {
      stop("K should be a positive integer!")
    }
    else {
      stop("K should be a positive integer! \n  Choose appropriate range of K!")
    }
  }
  if (kappa > 0.5 | kappa < 0) {
    if(verbose){cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")}
    kappa <- 0.5
  }
  if (select != "pls2" & select != "simpls") {
    if(verbose){cat("Invalid PLS algorithm for variable selection.\n")}
    if(verbose){cat("pls2 algorithm is used. \n\n")}
    select <- "pls2"
  }
  fits <- c("regression", "canonical", "invariant", "classic")
  if (!any(fit == fits)) {
    if(verbose){cat("Invalid PLS algorithm for model fitting\n")}
    if(verbose){cat("regression algorithm is used. \n\n")}
    fit <- "regression"
  }
  list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
}

spls.dv <- function (Z, eta, kappa, eps, maxstep) 
{
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median(abs(Z))
  Z <- Z/Znorm1
  if (q == 1) {
    c <- ust(Z, eta)
  }
  if (q > 1) {
    M <- Z %*% t(Z)
    dis <- 10
    i <- 1
    if (kappa == 0.5) {
      c <- matrix(10, p, 1)
      c.old <- c
      while (dis > eps & i <= maxstep) {
        mcsvd <- svd(M %*% c)
        a <- mcsvd$u %*% t(mcsvd$v)
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
    if (kappa > 0 & kappa < 0.5) {
      kappa2 <- (1 - kappa)/(1 - 2 * kappa)
      c <- matrix(10, p, 1)
      c.old <- c
      h <- function(lambda) {
        alpha <- solve(M + lambda * diag(p)) %*% M %*% 
          c
        obj <- t(alpha) %*% alpha - 1/kappa2^2
        return(obj)
      }
      if (h(eps) * h(1e+30) > 0) {
        while (h(eps) <= 1e+05) {
          M <- 2 * M
          c <- 2 * c
        }
      }
      while (dis > eps & i <= maxstep) {
        if (h(eps) * h(1e+30) > 0) {
          while (h(eps) <= 1e+05) {
            M <- 2 * M
            c <- 2 * c
          }
        }
        lambdas <- uniroot(h, c(eps, 1e+30))$root
        a <- kappa2 * solve(M + lambdas * diag(p)) %*% 
          M %*% c
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
  }
  return(c)
}


# Helper to coerce prediction inputs while supporting big.matrix objects
.bigPLScox_coerce_newdata <- function(newdata, p) {
  if (inherits(newdata, "big.matrix")) {
    if (!requireNamespace("bigmemory", quietly = TRUE)) {
      stop("Package 'bigmemory' is required to handle big.matrix inputs")
    }
    newdata <- bigmemory::as.matrix(newdata[, , drop = FALSE])
  }
  if (is.data.frame(newdata)) {
    newdata <- data.matrix(newdata)
  }
  if (is.null(dim(newdata))) {
    if (length(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ",
           p, " or a vector of length = ", p, ".")
    dim(newdata) <- c(1, p)
  } else {
    newdata <- as.matrix(newdata)
    if (ncol(newdata) != p)
      stop("'newdata' must be a numeric matrix with ncol = ",
           p, " or a vector of length = ", p, ".")
  }
  storage.mode(newdata) <- "double"
  newdata
}
