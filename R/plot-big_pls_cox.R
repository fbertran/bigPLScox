#' Plot method for big_pls_cox objects
#'
#' @param x A \code{big_pls_cox} object.
#' @param which Type of plot: "scores", "loadings", or "risk_groups".
#' @param comps Components to use (for "scores" and "loadings").
#' @param groups Optional grouping factor for scores plot; if NULL and
#'   \code{status} available, groups are derived from event status.
#' @param breaks Number of risk groups for "risk_groups" (default 3).
#' @param ... Further graphical arguments (passed to base plotting functions).
#'
#' @export
plot.big_pls_cox <- function(x,
                             which = c("scores", "loadings", "risk_groups"),
                             comps = c(1, 2),
                             groups = NULL,
                             breaks = 3,
                             ...) {
  if (!inherits(x, "big_pls_cox")) {
    stop("plot.big_pls_cox() requires an object of class 'big_pls_cox'.")
  }
  which <- match.arg(which)
  
  scores  <- x$scores
  loadings <- x$loadings
  
  if (which == "scores") {
    if (length(comps) != 2L) {
      stop("'comps' must be a length-2 vector for 'scores' plot.")
    }
    c1 <- comps[1]; c2 <- comps[2]
    if (c1 < 1 || c1 > ncol(scores) || c2 < 1 || c2 > ncol(scores)) {
      stop("Component indices in 'comps' out of bounds.")
    }
    
    if (is.null(groups)) {
      if (!is.null(x$status)) {
        groups <- factor(x$status, labels = c("censored", "event"))
      } else {
        groups <- factor(rep(1, nrow(scores)))
      }
    } else {
      groups <- factor(groups)
    }
    
    cols <- grDevices::rainbow(length(levels(groups)))
    plot(scores[, c1], scores[, c2],
         col = cols[groups],
         pch = 19,
         xlab = paste0("Component ", c1),
         ylab = paste0("Component ", c2),
         ...)
    legend("topright", legend = levels(groups), col = cols, pch = 19, bty = "n")
    
  } else if (which == "loadings") {
    c1 <- comps[1]
    if (c1 < 1 || c1 > ncol(loadings)) {
      stop("Component index in 'comps' out of bounds.")
    }
    barplot(
      height = loadings[, c1],
      border = NA,
      main = paste("Loadings - Component", c1),
      xlab = "Predictor index",
      ylab = "Loading",
      ...
    )
    
  } else if (which == "risk_groups") {
    if (is.null(x$cox_fit)) {
      stop("No 'cox_fit' element found; cannot plot risk groups.")
    }
    eta <- stats::predict(x$cox_fit, type = "link")
    if (is.null(x$time) || is.null(x$status)) {
      stop("Object must contain 'time' and 'status' for risk_groups plot.")
    }
    
    # group by linear predictor quantiles
    g <- cut(eta,
             breaks = stats::quantile(eta, probs = seq(0, 1, length.out = breaks + 1)),
             include.lowest = TRUE)
    sf <- survival::survfit(survival::Surv(x$time, x$status) ~ g)
    plot(sf, col = 1:breaks, lwd = 2, xlab = "Time", ylab = "Survival probability", ...)
    legend("topright", legend = levels(g), col = 1:breaks, lwd = 2, bty = "n",
           title = "Risk group")
  }
  
  invisible(x)
}

#' Plot method for big_pls_cox_gd objects
#'
#' @inheritParams plot.big_pls_cox
#' @export
plot.big_pls_cox_gd <- function(x,
                                which = c("scores", "loadings", "risk_groups"),
                                comps = c(1, 2),
                                groups = NULL,
                                breaks = 3,
                                ...) {
  if (!inherits(x, "big_pls_cox_gd")) {
    stop("plot.big_pls_cox_gd() requires an object of class 'big_pls_cox_gd'.")
  }
  which <- match.arg(which)
  
  scores   <- x$scores
  loadings <- x$loadings
  
  if (which == "scores") {
    if (length(comps) != 2L) {
      stop("'comps' must be a length-2 vector for 'scores' plot.")
    }
    c1 <- comps[1]; c2 <- comps[2]
    if (c1 < 1 || c1 > ncol(scores) || c2 < 1 || c2 > ncol(scores)) {
      stop("Component indices in 'comps' out of bounds.")
    }
    
    if (is.null(groups)) {
      if (!is.null(x$status)) {
        groups <- factor(x$status, labels = c("censored", "event"))
      } else {
        groups <- factor(rep(1, nrow(scores)))
      }
    } else {
      groups <- factor(groups)
    }
    
    cols <- grDevices::rainbow(length(levels(groups)))
    plot(scores[, c1], scores[, c2],
         col = cols[groups],
         pch = 19,
         xlab = paste0("Component ", c1),
         ylab = paste0("Component ", c2),
         ...)
    legend("topright", legend = levels(groups), col = cols, pch = 19, bty = "n")
    
  } else if (which == "loadings") {
    c1 <- comps[1]
    if (c1 < 1 || c1 > ncol(loadings)) {
      stop("Component index in 'comps' out of bounds.")
    }
    barplot(
      height = loadings[, c1],
      border = NA,
      main = paste("Loadings - Component", c1),
      xlab = "Predictor index",
      ylab = "Loading",
      ...
    )
    
  } else if (which == "risk_groups") {
    # For GD we only have coefficients; build a linear predictor directly
    if (is.null(x$coefficients)) {
      stop("No 'coefficients' element found; cannot plot risk groups.")
    }
    if (is.null(x$time) || is.null(x$status)) {
      stop("Object must contain 'time' and 'status' for risk_groups plot.")
    }
    eta <- drop(x$scores %*% x$coefficients)
    
    g <- cut(eta,
             breaks = stats::quantile(eta, probs = seq(0, 1, length.out = breaks + 1)),
             include.lowest = TRUE)
    sf <- survival::survfit(survival::Surv(x$time, x$status) ~ g)
    plot(sf, col = 1:breaks, lwd = 2, xlab = "Time", ylab = "Survival probability", ...)
    legend("topright", legend = levels(g), col = 1:breaks, lwd = 2, bty = "n",
           title = "Risk group")
  }
  
  invisible(x)
}
