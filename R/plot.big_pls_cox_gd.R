#' Basic diagnostic plots for big_pls_cox_gd
#' @param what Character string specifying what to plot.
#'   Typically `"scores"`, `"loadings"` or `"coefficients"`.
#'   
#' @export
plot.big_pls_cox_gd <- function(x, what = c("loglik", "gradient", "steps", "coef", "risk", "all"), ...) {
  what <- match.arg(what)
  d <- gd_diagnostics(x)
  
  if (what == "loglik" || what == "all") {
    plot(d$iterations, d$loglik, type="l", col="blue", lwd=2,
         xlab="Iteration", ylab="Log-likelihood",
         main="Log-likelihood progression")
    if (what != "all") return(invisible(NULL))
  }
  
  if (what == "gradient" || what == "all") {
    plot(d$iterations, d$gradnorm, type="l", col="red", lwd=2,
         xlab="Iteration", ylab="||Gradient||",
         main="Gradient norm progression")
    if (what != "all") return(invisible(NULL))
  }
  
  if (what == "steps" || what == "all") {
    plot(d$iterations, d$steps, type="l", col="darkgreen", lwd=2,
         xlab="Iteration", ylab="Step size",
         main="Line-search step size")
    if (what != "all") return(invisible(NULL))
  }
  
  if (what == "coef" || what == "all") {
    matplot(do.call(rbind, d$coef_trace), type="l", lwd=2,
            xlab="Iteration", ylab="Coefficient value",
            main="Coefficient trace (per component)")
    if (what != "all") return(invisible(NULL))
  }
  
  if (what == "risk" || what == "all") {
    eta_list <- d$eta_trace
    corrs <- sapply(2:length(eta_list), function(i)
      cor(eta_list[[i]], eta_list[[i-1]])
    )
    plot(2:length(eta_list), corrs, type="l", col="purple", lwd=2,
         xlab="Iteration", ylab="Correlation",
         main="Correlation of linear predictors across iterations")
    if (what != "all") return(invisible(NULL))
  }
  
  invisible(NULL)
}
