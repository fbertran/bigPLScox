#' Fitting a Cox-Model on group sparse PLSR components
#' 
#' This function computes the Cox Model based on PLSR components computed model
#' with \itemize{\item as the response: the Survival time \item as explanatory
#' variables: Xplan.  } It uses the package \code{sgPLS} to perform group PLSR
#' fit.
#' 
#' If \code{allres=FALSE} returns only the final Cox-model. If
#' \code{allres=TRUE} returns a list with the PLS components, the final
#' Cox-model and the group PLSR model. \code{allres=TRUE} is useful for evluating
#' model prediction accuracy on a test sample.
#' 
#' @aliases coxsgpls coxsgpls.default coxsgpls.formula
#' @param Xplan a formula or a matrix with the eXplanatory variables (training)
#' dataset
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
#' @param scaleX Should the \code{Xplan} columns be standardized ?
#' @param scaleY Should the \code{time} values be standardized ?
#' @param ncomp The number of components to include in the model. It this is
#' not supplied, min(7,maximal number) components is used.
#' @param ind.block.x a vector of integers describing the grouping of the 
#' X-variables. \code{ind.block.x <- c(3,10,15)} means that \code{X} is 
#' structured into 4 groups: \code{X1} to \code{X3}; \code{X4} to \code{X10}, 
#' \code{X11} to \code{X15} and \code{X16} to \code{Xp} where \code{p} is the 
#' number of variables in the X matrix.
#' @param keepX	numeric vector of length ncomp, the number of variables to keep 
#' in X-loadings. By default all variables are kept in the model.
#' @param alpha.x The mixing parameter (value between 0 and 1) related to the sparsity 
#' within group for the X dataset.
#' @param upper.lambda By default \code{upper.lambda=10^5}. A large value 
#' specifying the upper bound of the intervall of lambda values for searching 
#' the value of the tuning parameter (lambda) corresponding to a non-zero group 
#' of variables.
#' @param modepls character string. What type of algorithm to use, (partially)
#' matching one of "regression", "canonical". See
#' \code{\link[sgPLS]{gPLS}} for details
#' @param plot Should the survival function be plotted ?)
#' @param allres FALSE to return only the Cox model and TRUE for additionnal
#' results. See details. Defaults to FALSE.
#' @param dataXplan an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in \code{dataXplan}, the variables are
#' taken from \code{environment(Xplan)}, typically the environment from which
#' \code{coxpls} is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights an optional vector of 'prior weights' to be used in the
#' fitting process. Should be \code{NULL} or a numeric vector.
#' @param model_frame If \code{TRUE}, the model frame is returned.
#' @param model_matrix If \code{TRUE}, the model matrix is returned.
#' @param contrasts.arg a list, whose entries are values (numeric matrices, 
#' functions or character strings naming functions) to be used as replacement 
#' values for the contrasts replacement function and whose names are the names 
#' of columns of data containing factors.
#' @param \dots Arguments to be passed on to \code{survival::coxph}.
#' @return If \code{allres=FALSE} : \item{cox_sgpls}{Final Cox-model.} If
#' \code{allres=TRUE} : \item{tt_sgpls}{PLSR components.} \item{cox_sgpls}{Final
#' Cox-model.} \item{sgpls_mod}{The PLSR model.}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[survival]{coxph}}, \code{\link[sgPLS]{gPLS}}
#' @references A group and Sparse Group Partial Least Square approach applied 
#' in Genomics context, Liquet Benoit, Lafaye de Micheaux, Boris Hejblum, 
#' Rodolphe Thiebaut (2016). Bioinformatics.\cr
#' 
#' Deviance residuals-based sparse PLS and sparse kernel PLS regression for
#' censored data, Philippe Bastien, Frederic Bertrand, Nicolas Meyer and Myriam
#' Maumy-Bertrand (2015), Bioinformatics, 31(3):397-404,
#' doi:10.1093/bioinformatics/btu660.
#' @keywords models regression
#' @examples
#' 
#' data(micro.censure)
#' data(Xmicro.censure_compl_imp)
#' 
#' X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
#' X_train_micro_df <- data.frame(X_train_micro)
#' Y_train_micro <- micro.censure$survyear[1:80]
#' C_train_micro <- micro.censure$DC[1:80]
#' 
#' (coxsgpls_fit=coxsgpls(X_train_micro,Y_train_micro,C_train_micro,
#' ncomp=6,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 6)))
#' (coxsgpls_fit=coxsgpls(~X_train_micro,Y_train_micro,C_train_micro,
#' ncomp=6,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 6)))
#' (coxsgpls_fit=coxsgpls(~.,Y_train_micro,C_train_micro,ncomp=6,
#' dataXplan=X_train_micro_df,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 6)))
#' 
#' rm(X_train_micro,Y_train_micro,C_train_micro,cox_sgpls_sgfit)
#' 
#' @export coxsgpls
coxsgpls <- function (Xplan, ...) UseMethod("coxsgpls")

#' @rdname coxsgpls
#' @export
coxsgpls.formula <-
  function (Xplan, time, time2, event, type, origin, typeres = "deviance", 
            collapse, weighted, scaleX = TRUE, scaleY = TRUE, ncomp = 
              min(7, ncol(Xplan)), modepls = "regression", ind.block.x, keepX, 
            alpha.x, upper.lambda = 10 ^ 5, plot = FALSE, allres = FALSE, 
            dataXplan = NULL, subset, weights, model_frame = FALSE, 
            model_matrix = FALSE, contrasts.arg=NULL, ...) 
  {
    if (missing(dataXplan)) 
      dataXplan <- environment(Xplan)
    mf0 <- match.call(expand.dots = FALSE)
    m0 <- match(c("subset", "weights"), names(mf0), 0L)
    mf0 <- mf0[c(1L, m0)]
    mf0$data <- dataXplan
    mf0$formula <- as.formula(paste(c(as.character(Xplan),"+0"),collapse=""))
    mf0$drop.unused.levels <- TRUE
    mf0[[1L]] <- as.name("model.frame")
    mf0 <- eval(mf0, parent.frame())
    if (model_frame) 
      return(mf0)
    mt0 <- attr(mf0, "terms")
    Y <- model.response(mf0, "any")
    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
        names(Y) <- nm
    }
    Xplan <- if (!is.empty.model(mt0)) 
      model.matrix(mt0, mf0, 
                   contrasts.arg=
                     contrasts.arg)
    else matrix(, NROW(Y), 0L)
    if (model_matrix) 
      return(model.matrix(mt0, mf0, 
                          contrasts.arg=
                            contrasts.arg))
#    ind.block.x <- sapply(ind.block.x, function(x) {sum(attr(bbb,"assign") <= x)})
    weights <- as.vector(model.weights(mf0))
    if (!is.null(weights) && !is.numeric(weights)) 
      stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
      stop("negative weights not allowed")
    NextMethod("coxsgpls")
  }


#' @rdname coxsgpls
#' @export
coxsgpls.default <- 
  function (Xplan, time, time2, event, type, origin, typeres = "deviance", 
            collapse, weighted, scaleX = TRUE, scaleY = TRUE, 
            ncomp = min(7,ncol(Xplan)), modepls = "regression", ind.block.x, 
            keepX, alpha.x, upper.lambda = 10 ^ 5, plot = FALSE, 
            allres = FALSE, ...) 
  {
    if (scaleX) {
      Xplan <- scale(Xplan)
      XplanScal <- attr(Xplan, "scaled:scale")
      XplanCent <- attr(Xplan, "scaled:center")
      Xplan <- as.data.frame(Xplan)
    }
    else {
      Xplan <- as.data.frame(Xplan)
      XplanScal <- rep(1, ncol(Xplan))
      XplanCent <- rep(0, ncol(Xplan))
    }
    if ((scaleY & missing(time2))) {
      time <- scale(time)
    }
    try(attachNamespace("survival"), silent = TRUE)
    suppressMessages(try(attachNamespace("sgPLS"), silent = TRUE))
    on.exit(try(unloadNamespace("sgPLS"), silent = TRUE), 
            add = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("time", "time2", "event", "type", "origin"), 
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("Surv")
    YCsurv <- eval(mf, parent.frame())
    mf2 <- match.call(expand.dots = FALSE)
    m2 <- match(c("ncomp"), names(mf2), 0L)
    mf2 <- mf2[c(1L, m2)]
    mf2$ncomp <- eval.parent(mf2$ncomp)
    mf2$X <- eval.parent(Xplan)
    mf2$Y <- eval.parent(time)
    mf2$mode <- eval.parent(modepls)
    mf2$ind.block.x <- eval.parent(ind.block.x)
    mf2$keepX <- if (missing(keepX)) {rep(length(mf2$ind.block.x)+1,mf2$ncomp)} else {eval.parent(keepX)}
    mf2$alpha.x <- eval.parent(alpha.x)
    mf2$upper.lambda <- eval.parent(upper.lambda)
    mf2$scale = FALSE
    mf2[[1L]] <- as.name("sgPLS")
    if (mf2$ncomp == 0) {
      sgpls_mod <- NULL
    }
    else {
      sgpls_mod <- eval(mf2)
    }
    tt_sgpls <- data.frame(sgpls_mod$variates$X)
    if (mf2$ncomp > 0) {
      colnames(tt_sgpls) <- paste("dim", 1:ncol(tt_sgpls), sep = ".")
    }
    if (mf2$ncomp == 0) {
      mf2b <- match.call(expand.dots = TRUE)
      m2b <- match(c(head(names(as.list(args(survival::coxph))), -2), 
                     head(names(as.list(args(survival::coxph.control))), -1)), 
                   names(mf2b), 0L)
      mf2b <- mf2b[c(1L, m2b)]
      mf2b$formula <- as.formula(YCsurv ~ 1)
      mf2b$data <- tt_sgpls
      mf2b[[1L]] <- as.name("coxph")
      cox_sgpls <- eval(mf2b, parent.frame())
      cox_sgpls$call$data <- as.name("tt_sgpls")
    }
    else {
      mf2b <- match.call(expand.dots = TRUE)
      m2b <- match(c(head(names(as.list(args(survival::coxph))), -2), 
                     head(names(as.list(args(survival::coxph.control))), -1)), 
                   names(mf2b), 0L)
      mf2b <- mf2b[c(1L, m2b)]
      mf2b$formula <- as.formula(YCsurv ~ .)
      mf2b$data <- tt_sgpls
      mf2b[[1L]] <- as.name("coxph")
      cox_sgpls <- eval(mf2b, parent.frame())
      cox_sgpls$call$data <- as.name("tt_sgpls")
    }
    if (!allres) {
      return(cox_sgpls)
    }
    else {
      CoeffCFull = matrix(NA, nrow = ncomp, ncol = ncomp)
      if (mf2$ncomp > 0) {
        for (iii in 1:ncomp) {
          mf2b <- match.call(expand.dots = TRUE)
          m2b <- match(c(head(names(as.list(args(survival::coxph))),
                              -2), head(names(as.list(args(survival::coxph.control))),
                                        -1)), names(mf2b), 0L)
          mf2b <- mf2b[c(1L, m2b)]
          mf2b$formula <- as.formula(YCsurv ~ .)
          mf2b$data <- tt_sgpls[, 1:iii, drop = FALSE]
          mf2b[[1L]] <- as.name("coxph")
          cox_sgpls <- eval(mf2b, parent.frame())
          cox_sgpls$call$data <- as.name("tt_sgpls")
          CoeffCFull[, iii] <- c(cox_sgpls$coefficients,
                                 rep(NA, ncomp - iii))
        }
      }
      res <- list(tt_sgpls = tt_sgpls, cox_sgpls = cox_sgpls, sgpls_mod = sgpls_mod,
                  XplanScal = XplanScal, XplanCent = XplanCent,
                  CoeffCFull = CoeffCFull)
      res$XplanTrain <- as.matrix(Xplan)
      class(res) <- c("coxsgpls", "cox_pls_legacy")
      return(res)
    }
  }
