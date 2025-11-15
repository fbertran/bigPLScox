#' Align a GD fit to a PLS fit (optional refit)
#' @param fit_gd  object from big_pls_cox_gd()
#' @param fit_pls object from big_pls_cox()
#' @param time,status Surv parts used for refit
#' @param rotate  logical; Procrustes-rotate GD scores to the PLS basis
#' @param refit   logical; refit a Cox model on the rotated GD scores
#' @return fit_gd with $scores/$coefficients/$cox_fit possibly updated
align_big_plscox <- function(fit_gd, fit_pls, time, status, rotate = TRUE, refit = TRUE) {
  Sg <- fit_gd$scores
  Sp <- fit_pls$scores
  stopifnot(is.matrix(Sg), is.matrix(Sp), nrow(Sg) == nrow(Sp), ncol(Sg) == ncol(Sp))
  
  if (rotate) {
    C  <- crossprod(Sg, Sp)            # p x p
    sv <- svd(C)
    R  <- sv$u %*% t(sv$v)             # orthogonal rotation
    Sg <- Sg %*% R
  }
  
  # Optional refit to put coefficients on the same scale
  if (refit) {
    cf <- survival::coxph(survival::Surv(time, status) ~ Sg, ties = "efron", x = FALSE)
    fit_gd$coefficients <- unname(coef(cf))
    fit_gd$cox_fit      <- cf
  }
  
  # Replace scores (and keep other slots untouched)
  fit_gd$scores <- Sg
  
  # Harmonize overall sign (doesn't change Cox partial likelihood)
  eta_gd  <- drop(Sg %*% fit_gd$coefficients)
  eta_pls <- drop(Sp %*% fit_pls$cox_fit$coefficients)
  if (cor(eta_gd, eta_pls) < 0) {
    fit_gd$scores       <- -fit_gd$scores
    fit_gd$coefficients <- -fit_gd$coefficients
    if (!is.null(fit_gd$cox_fit)) {
      fit_gd$cox_fit$linear.predictors <- -fit_gd$cox_fit$linear.predictors
      fit_gd$cox_fit$coefficients      <- -fit_gd$cox_fit$coefficients
    }
  }
  
  fit_gd
}