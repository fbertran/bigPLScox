# Internal helpers and S3 predict methods for legacy Cox-PLS fits

legacy_cox_scale_newdata <- function(object, newdata) {
  if (is.null(object$XplanCent) || is.null(object$XplanScal)) {
    stop("The fitted object does not contain centering/scaling information required for prediction.")
  }
  p <- length(object$XplanCent)
  new_mat <- .bigPLScox_coerce_newdata(newdata, p)
  scale(new_mat, center = object$XplanCent, scale = object$XplanScal)
}

legacy_cox_training_scores <- function(object, score_name) {
  scores <- object[[score_name]]
  if (is.null(scores)) {
    n_obs <- if (!is.null(object$XplanTrain)) nrow(as.matrix(object$XplanTrain)) else 0L
    matrix(numeric(0), nrow = n_obs, ncol = 0,
           dimnames = list(NULL, NULL))
  } else {
    as.matrix(scores)
  }
}

legacy_cox_compute_scores <- function(object, newdata, comps, score_name, mod_name,
                                      predict_fun, kernel_slot = NULL) {
  training_scores <- legacy_cox_training_scores(object, score_name)
  total_comp <- if (is.null(dim(training_scores))) {
    if (length(training_scores)) 1L else 0L
  } else {
    ncol(training_scores)
  }
  if (total_comp == 1L && is.null(dim(training_scores))) {
    training_scores <- matrix(training_scores, ncol = 1L)
  }
  if (is.null(comps)) {
    comps <- if (total_comp > 0L) seq_len(total_comp) else integer(0L)
  }
  comps <- as.integer(comps)
  if (length(comps) > 0L) {
    if (any(comps < 1L) || any(comps > total_comp)) {
      stop("`comps` indices are out of bounds")
    }
  } else if (total_comp == 0L && length(comps) > 0L) {
    stop("No components were fitted, cannot select from `comps`.")
  }
  
  if (is.null(newdata)) {
    if (length(comps) == 0L) {
      scores <- matrix(numeric(0), nrow = nrow(training_scores), ncol = 0,
                       dimnames = list(rownames(training_scores), NULL))
    } else {
      scores <- training_scores[, comps, drop = FALSE]
    }
    return(list(scores = scores, comps = comps, total = total_comp))
  }
  
  if (is.null(object[[mod_name]]) || total_comp == 0L) {
    scaled <- legacy_cox_scale_newdata(object, newdata)
    scores <- matrix(numeric(0), nrow = nrow(scaled), ncol = 0,
                     dimnames = list(rownames(scaled), NULL))
    return(list(scores = scores, comps = comps, total = total_comp))
  }
  
  scaled <- legacy_cox_scale_newdata(object, newdata)
  if (!is.null(kernel_slot)) {
    if (!requireNamespace("kernlab", quietly = TRUE)) {
      stop("Package 'kernlab' is required to score kernel-based models")
    }
    if (is.null(object$XplanTrain)) {
      stop("The fitted object does not store the scaled training design matrix required for kernel predictions")
    }
    kernel_obj <- object[[kernel_slot]]
    if (is.null(kernel_obj)) {
      stop("Kernel definition is missing from the fitted object")
    }
    train_scaled <- as.matrix(object$XplanTrain)
    kernel_mat <- kernlab::kernelMatrix(kernel_obj, as.matrix(scaled), train_scaled)
    preds <- predict_fun(object[[mod_name]], newdata = kernel_mat,
                         scale.X = FALSE, scale.Y = FALSE)
  } else {
    preds <- predict_fun(object[[mod_name]], newdata = scaled,
                         scale.X = FALSE, scale.Y = FALSE)
  }
  variates <- if (is.list(preds) && !is.null(preds$variates)) preds$variates else preds
  variates <- as.matrix(variates)
  if (length(comps) == 0L) {
    scores <- matrix(numeric(0), nrow = nrow(variates), ncol = 0,
                     dimnames = list(rownames(variates), NULL))
  } else {
    scores <- variates[, comps, drop = FALSE]
  }
  list(scores = scores, comps = comps, total = total_comp)
}

legacy_predict_cox <- function(object, newdata, type, comps, coef,
                               score_name, cox_name, mod_name, predict_fun,
                               kernel_slot = NULL) {
  comp_info <- legacy_cox_compute_scores(object, newdata, comps, score_name,
                                         mod_name, predict_fun, kernel_slot)
  scores <- comp_info$scores
  comps <- comp_info$comps
  if (type == "components") {
    return(scores)
  }
  n_obs <- if (length(scores)) nrow(scores) else {
    if (is.null(newdata)) {
      nrow(legacy_cox_training_scores(object, score_name))
    } else {
      nrow(legacy_cox_scale_newdata(object, newdata))
    }
  }
  if (length(comps) == 0L) {
    eta <- rep(0, n_obs)
  } else {
    if (is.null(coef)) {
      cox_fit <- object[[cox_name]]
      if (is.null(cox_fit) || is.null(cox_fit$coefficients)) {
        stop("Cox model coefficients are not available; supply them via `coef`")
      }
      coef <- cox_fit$coefficients
    }
    coef <- as.numeric(coef)
    if (length(coef) < max(comps)) {
      stop("`coef` must have at least max(comps) entries")
    }
    eta <- as.numeric(scores %*% coef[comps])
  }
  if (!is.null(rownames(scores)) && length(eta) == nrow(scores)) {
    names(eta) <- rownames(scores)
  }
  switch(type,
         link = eta,
         risk = exp(eta),
         response = exp(eta))
}

#' Predict survival summaries from legacy Cox-PLS fits
#'
#' These methods extend [stats::predict()] for Cox models fitted with the
#' original PLS engines exposed by [coxgpls()], [coxsgpls()], and their
#' deviance-residual or kernel variants. They provide access to latent component
#' scores alongside linear predictors and risk estimates, ensuring consistent
#' behaviour with the newer big-memory solvers.
#'
#' @param object A fitted model returned by [coxgpls()], [coxsgpls()],
#'   [coxspls_sgpls()], or any of their deviance-residual/kernel counterparts
#'   with `allres = TRUE`.
#' @param newdata Optional matrix or data frame of predictors. When `NULL`, the
#'   training components stored in `object` are reused.
#' @param type Type of prediction requested: `"link"` for linear predictors,
#'   `"risk"`/`"response"` for exponentiated scores, or `"components"` to return
#'   latent PLS scores.
#' @param comps Optional integer vector specifying which latent components to
#'   retain. Defaults to all available components.
#' @param coef Optional coefficient vector overriding the Cox model
#'   coefficients stored in `object`.
#' @param ... Unused arguments for future extensions.
#'
#' @return When `type` is `"components"`, a matrix of latent scores; otherwise a
#'   numeric vector containing the requested prediction with names inherited from
#'   the supplied data.
#'
#' @seealso [coxgpls()], [coxsgpls()], [coxspls_sgpls()],
#'   [coxDKgplsDR()], [predict.big_pls_cox()], [computeDR()].
#'
#' @references
#'   Bastien, P., Bertrand, F., Meyer, N., & Maumy-Bertrand, M. (2015).
#'   Deviance residuals-based sparse PLS and sparse kernel PLS for censored
#'   data. *Bioinformatics*, 31(3), 397–404. <doi:10.1093/bioinformatics/btu660>
#'
#'   Bertrand, F., Bastien, P., & Maumy-Bertrand, M. (2018).
#'   Cross validating extensions of kernel, sparse or regular partial least
#'   squares regression models to censored data. <https://arxiv.org/abs/1810.01005>.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   data(micro.censure, package = "bigPLScox")
#'   data(Xmicro.censure_compl_imp, package = "bigPLScox")
#'
#'   X <- as.matrix(Xmicro.censure_compl_imp[1:60, 1:10])
#'   time <- micro.censure$survyear[1:60]
#'   status <- micro.censure$DC[1:60]
#'
#'   set.seed(321)
#'   fit <- coxgpls(
#'     Xplan = X,
#'     time = time,
#'     status = status,
#'     ncomp = 2,
#'     allres = TRUE
#'   )
#'
#'   predict(fit, newdata = X[1:5, ], type = "risk")
#'   head(predict(fit, type = "components"))
#' }
#'
#' @name predict_cox_pls
#' @rdname predict_cox_pls
#' @export
predict.coxgpls <- function(object, newdata = NULL,
                            type = c("link", "risk", "response", "components"),
                            comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_gpls", cox_name = "cox_gpls",
                     mod_name = "gpls_mod", predict_fun = predict.gPLS)
}

#' @rdname predict_cox_pls
#' @export
predict.coxgplsDR <- function(object, newdata = NULL,
                              type = c("link", "risk", "response", "components"),
                              comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_gplsDR", cox_name = "cox_gplsDR",
                     mod_name = "gplsDR_mod", predict_fun = predict.gPLS)
}

#' @rdname predict_cox_pls
#' @export
predict.coxsgpls <- function(object, newdata = NULL,
                             type = c("link", "risk", "response", "components"),
                             comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_sgpls", cox_name = "cox_sgpls",
                     mod_name = "sgpls_mod", predict_fun = predict.gPLS)
}

#' @rdname predict_cox_pls
#' @export
predict.coxsgplsDR <- function(object, newdata = NULL,
                               type = c("link", "risk", "response", "components"),
                               comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_sgplsDR", cox_name = "cox_sgplsDR",
                     mod_name = "sgplsDR_mod", predict_fun = predict.gPLS)
}

#' @rdname predict_cox_pls
#' @export
predict.coxspls_sgpls <- function(object, newdata = NULL,
                                  type = c("link", "risk", "response", "components"),
                                  comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_spls_sgpls", cox_name = "cox_spls_sgpls",
                     mod_name = "spls_sgpls_mod", predict_fun = predict.gPLS)
}

#' @rdname predict_cox_pls
#' @export
predict.coxDKgplsDR <- function(object, newdata = NULL,
                                type = c("link", "risk", "response", "components"),
                                comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_DKgplsDR", cox_name = "cox_DKgplsDR",
                     mod_name = "DKgplsDR_mod", predict_fun = predict.gPLS,
                     kernel_slot = "kernDKgplsDR_mod")
}

#' @rdname predict_cox_pls
#' @export
predict.coxDKsgplsDR <- function(object, newdata = NULL,
                                 type = c("link", "risk", "response", "components"),
                                 comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_DKsgplsDR", cox_name = "cox_DKsgplsDR",
                     mod_name = "DKsgplsDR_mod", predict_fun = predict.gPLS,
                     kernel_slot = "kernDKsgplsDR_mod")
}

#' @rdname predict_cox_pls
#' @export
predict.coxDKspls_sgplsDR <- function(object, newdata = NULL,
                                      type = c("link", "risk", "response", "components"),
                                      comps = NULL, coef = NULL, ...) {
  type <- match.arg(type)
  legacy_predict_cox(object, newdata, type, comps, coef,
                     score_name = "tt_DKspls_sgplsDR",
                     cox_name = "cox_DKspls_sgplsDR",
                     mod_name = "DKspls_sgplsDR_mod", predict_fun = predict.gPLS,
                     kernel_slot = "kernDKspls_sgplsDR_mod")
}
