#' Compute deviance residuals
#'
#' This function computes deviance residuals from a null Cox model. By default
#' it delegates to [`survival::coxph()`], but a high-performance C++ engine is
#' also available for large in-memory or [`bigmemory::big.matrix`] design
#' matrices.
#' 
#' @param time for right censored data, this is the follow up time. For
#' interval data, the first argument is the starting time for the interval.
#' @param time2 The status indicator, normally 0=alive, 1=dead. Other choices
#' are \code{TRUE/FALSE} (\code{TRUE} = death) or 1/2 (2=death). For interval
#' censored data, the status indicator is 0=right censored, 1=event at
#' \code{time}, 2=left censored, 3=interval censored. Although unusual, the
#' event indicator can be omitted, in which case all subjects are assumed to
#' have an event.
#' @param event ending time of the interval for interval censored or counting
#' process data only. Intervals are assumed to be open on the left and closed
#' on the right, \code{(start, end]}. For counting process data, event
#' indicates whether an event occurred at the end of the interval.
#' @param type character string specifying the type of censoring. Possible
#' values are \code{"right"}, \code{"left"}, \code{"counting"},
#' \code{"interval"}, or \code{"interval2"}. The default is \code{"right"} or
#' \code{"counting"} depending on whether the \code{time2} argument is absent
#' or present, respectively.
#' @param origin for counting process data, the hazard function origin. This
#' option was intended to be used in conjunction with a model containing time
#' dependent strata in order to align the subjects properly when they cross
#' over from one strata to another, but it has rarely proven useful.
#' @param typeres character string indicating the type of residual desired.
#' Possible values are \code{"martingale"}, \code{"deviance"}, \code{"score"},
#' \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and
#' \code{"scaledsch"}. Only enough of the string to determine a unique match is
#' required.
#' @param collapse vector indicating which rows to collapse (sum) over. In
#' time-dependent models more than one row data can pertain to a single
#' individual. If there were 4 individuals represented by 3, 1, 2 and 4 rows of
#' data respectively, then \code{collapse=c(1,1,1,2,3,3,4,4,4,4)} could be used
#' to obtain per subject rather than per observation residuals.
#' @param weighted if \code{TRUE} and the model was fit with case weights, then
#' the weighted residuals are returned.
#' @param scaleY Should the \code{time} values be standardized ?
#' @param plot Should the survival function be plotted ?
#' @param engine Either `"survival"` (default) to call
#'   [`survival::coxph()`] or `"cpp"` to use the C++ implementation.
#' @param method Tie handling to use with `engine = "cpp"`: either
#'   `"efron"` (default) or `"breslow"`.
#' @param X Optional design matrix used to compute the linear predictor when
#'   `engine = "cpp"`. Supports base matrices, data frames, and
#'   [`bigmemory::big.matrix`] objects.
#' @param coef Optional coefficient vector associated with `X` when
#'   `engine = "cpp"`.
#' @param eta Optional precomputed linear predictor passed directly to the C++
#'   engine.
#' @param center,scale Optional centring and scaling vectors applied to `X`
#'   before computing the linear predictor with the C++ engine.
#'   
#' @return Residuals from a null model fit. When `engine = "cpp"`, the returned
#'   vector has attributes `"martingale"`, `"cumhaz"`, and
#'   `"linear_predictor"`.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[survival]{coxph}}
#' @references
#'   Bastien, P., Bertrand, F., Meyer, N., and Maumy-Bertrand, M.
#'   (2015). Deviance residuals-based sparse PLS and sparse kernel PLS for
#'   binary classification and survival analysis. *BMC Bioinformatics*, 16, 211.
#'   
#' Therneau, T.M., Grambsch, P.M. (2000). *Modeling Survival Data: Extending the
#'   Cox Model*. Springer.
#'
#' @keywords models regression
#' @examples
#' 
#' data(micro.censure, package = "bigPLScox")
#' 
#' Y_train_micro <- micro.censure$survyear[1:80]
#' C_train_micro <- micro.censure$DC[1:80]
#' 
#' Y_DR <- computeDR(Y_train_micro,C_train_micro)
#' Y_DR <- computeDR(Y_train_micro,C_train_micro,plot=TRUE)
#' 
#' Y_cpp <- computeDR(
#'   Y_train_micro,
#'   C_train_micro,
#'   engine = "cpp",
#'   eta = rep(0, length(Y_train_micro))
#' )
#' 
#' Y_qcpp <- computeDR(
#'   Y_train_micro,
#'   C_train_micro,
#'   engine = "qcpp"
#' )
#' 
#' @export computeDR
computeDR <- function (time, time2, event, type, origin, typeres = "deviance",
                       collapse, weighted, scaleY = TRUE, plot = FALSE,
                       engine = c("survival", "cpp", "qcpp"),
                       method = c("efron", "breslow"),
                       X = NULL, coef = NULL, eta = NULL,
                       center = NULL, scale = NULL)
{
  try(attachNamespace("survival"), silent = TRUE)
  engine_missing <- missing(engine)
  engine <- match.arg(engine)
  method <- match.arg(method)
  
  simple_status <- if (missing(time2)) rep(1, length(time)) else time2
  simple_case <- missing(event) && missing(origin) && (missing(type) || type == "right") &&
    missing(collapse) && missing(weighted) && identical(typeres, "deviance") && !plot
  if (simple_case && (engine_missing || engine == "qcpp")) {
    time_use <- if (scaleY) {
      as.numeric(scale(time))
    } else {
      as.numeric(time)
    }
    return(cox_deviance_residuals(time_use, as.numeric(simple_status)))
  }
  
  if (simple_case && !(engine == "cpp" && !is.null(eta))) {
    warning("'engine' is set to '", engine, "' so the fast C++ backend is not used.",
            call. = FALSE)
  }

  if ((scaleY & missing(time2))) {
    time <- scale(time)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("time", "time2", "event", "type", "origin"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("Surv")
  YCsurv <- eval(mf, parent.frame())
  if (plot) {
    plot(survival::survfit(YCsurv ~ 1))
  }
  
  if (engine == "cpp") {
    surv_mat <- as.matrix(YCsurv)
    time_vec <- surv_mat[, 1]
    status_vec <- surv_mat[, ncol(surv_mat)]
    if (!is.null(eta)) {
      eta_vec <- as.numeric(eta)
      if (length(eta_vec) != length(time_vec)) {
        stop("`eta` must have the same length as `time`", call. = FALSE)
      }
      # if (simple_case) {
      #   details <- cox_deviance_details(time_vec, status_vec)
      #   dev <- details$deviance
      #   attr(dev, "martingale") <- details$martingale
      #   attr(dev, "cumhaz") <- details$cumulative_hazard
      #   attr(dev, "linear_predictor") <- eta_vec
      #   attr(dev, "names") <- 1:length(YCsurv)
      #   return(dev)
      # }
      res <- deviance_residuals_cpp(time_vec, status_vec, eta_vec, method)
    } else {
      if (is.null(X) || is.null(coef)) {
        stop("`X` and `coef` must be supplied when `eta` is not provided and engine = 'cpp'")
      }
      if (inherits(X, "big.matrix")) {
        res <- big_deviance_residuals_cpp(X@address, as.numeric(coef),
                                          time_vec, status_vec,
                                          if (!is.null(center)) as.numeric(center) else NULL,
                                          if (!is.null(scale)) as.numeric(scale) else NULL,
                                          method)
      } else {
        if (is.data.frame(X)) {
          X <- as.matrix(X)
        }
        if (!is.matrix(X)) {
          stop("`X` must be a matrix, data frame, or big.matrix")
        }
        storage.mode(X) <- "double"
        res <- matrix_deviance_residuals_cpp(X, as.numeric(coef),
                                             time_vec, status_vec,
                                             if (!is.null(center)) as.numeric(center) else NULL,
                                             if (!is.null(scale)) as.numeric(scale) else NULL,
                                             method)
      }
    }
    dev <- res$deviance
    attr(dev, "martingale") <- res$martingale
    attr(dev, "cumhaz") <- res$cumhaz
    attr(dev, "linear_predictor") <- res$linear_predictor
    attr(dev, "names") <- 1:length(YCsurv)
    return(dev)
  }
  mf1 <- match.call(expand.dots = TRUE)
  m1 <- match(c(head(names(as.list(args(survival::coxph))), -2), head(names(as.list(args(survival::coxph.control))),
                                                            -1)), names(mf1), 0L)
  mf1 <- mf1[c(1L, m1)]
  mf1$formula <- as.formula(YCsurv ~ 1)
  mf1[[1L]] <- as.name("coxph")
  coxDR <- eval(mf1, parent.frame())
  mf2 <- match.call(expand.dots = FALSE)
  m2 <- match(c("weighted", "collapse", "origin"), names(mf2), 
              0L)
  mf2 <- mf2[c(1L, m2)]
  mf2$type <- typeres
  mf2$object <- coxDR
  mf2[[1L]] <- as.name("residuals")
  DR_coxph <- eval(mf2, parent.frame())
  return(DR_coxph)  
}
