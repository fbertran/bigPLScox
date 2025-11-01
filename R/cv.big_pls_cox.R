#' Cross-validation for big-memory PLS-Cox models
#'
#' @description Performs K-fold cross-validation for models fitted with
#'   [big_pls_cox()] or [big_pls_cox_gd()]. The routine mirrors the behaviour of
#'   the cross-validation helpers available in the original \pkg{plsRcox}
#'   package while operating on `big.matrix` inputs.
#'
#' @param data A list with entries `x`, `time` and `status` matching the
#'   arguments of [big_pls_cox()] or [big_pls_cox_gd()]. `x` can be either a
#'   numeric matrix/data frame or a [`bigmemory::big.matrix`].
#' @param nfold Integer giving the number of folds to use.
#' @param nt Number of latent components to evaluate.
#' @param keepX Optional integer vector passed to the modelling function to
#'   enforce naive sparsity (see [big_pls_cox()]).
#' @param givefold Optional list of fold indices. When supplied, it must contain
#'   `nfold` integer vectors whose union is `seq_len(nrow(data$x))`.
#' @param allCVcrit Logical; when `FALSE` (default) only the recommended
#'   integrated AUC computed with \pkg{survivalROC} is returned. When `TRUE`, the
#'   13 additional criteria from \pkg{plsRcox} are also evaluated.
#' @param times.auc Optional time grid used for time-dependent AUC computations.
#'   Defaults to an equally spaced grid between zero and the maximum observed
#'   time.
#' @param times.prederr Optional time grid used for prediction error curves.
#'   Defaults to the same grid as `times.auc` without the last ten evaluation
#'   points to avoid instabilities.
#' @param method Ties handling method passed to [`survival::coxph`].
#' @param verbose Logical; print progress information.
#' @param ... Additional arguments forwarded to the underlying modelling
#'   function.
#'
#' @details The function returns cross-validated estimates for each component
#'   (including the null model) using either [big_pls_cox()] or
#'   [big_pls_cox_gd()], depending on the `engine` argument. The implementation
#'   reuses the internal indicators (`getIndicCV`, `getIndicCViAUCSurvROCTest`)
#'   to provide consistent metrics with the legacy \pkg{plsRcox} helpers.
#'
#' @return A list containing cross-validation summaries. When `allCVcrit =
#'   FALSE`, the list holds
#'     \item{`nt`}{Number of components assessed.}
#'     \item{`cv.error10`}{Mean iAUC of \pkg{survivalROC} across folds for 0 to
#'       `nt` components.}
#'     \item{`cv.se10`}{Estimated standard errors for `cv.error10`.}
#'     \item{`folds`}{Fold assignments.}
#'     \item{`lambda.min10`}{Component minimising the cross-validated error.}
#'     \item{`lambda.1se10`}{Largest component within one standard error of the
#'       optimum.}
#'   When `allCVcrit = TRUE`, the full set of 14 criteria (log partial
#'   likelihood, iAUC variants and Brier scores) is returned together with their
#'   associated standard errors and one-standard-error selections.
#'
#' @export
cv.big_pls_cox <- function(data, nfold = 5L, nt = 5L, keepX = NULL,
                           givefold, allCVcrit = FALSE,
                           times.auc = NULL,
                           times.prederr = NULL,
                           method = c("efron", "breslow"),
                           verbose = TRUE, ...) {
  method <- match.arg(method)
  engine <- "pls"
  cv_big_pls_core(data = data, nfold = nfold, nt = nt, keepX = keepX,
                  givefold = givefold, allCVcrit = allCVcrit,
                  times.auc = times.auc, times.prederr = times.prederr,
                  method = method, verbose = verbose, engine = engine, ...)
}

#' @rdname cv.big_pls_cox
#' @export
cv.big_pls_cox_gd <- function(data, nfold = 5L, nt = NULL, keepX = NULL,
                              givefold, allCVcrit = FALSE,
                              times.auc = NULL,
                              times.prederr = NULL,
                              method = c("efron", "breslow"),
                              verbose = TRUE, ...) {
  method <- match.arg(method)
  engine <- "gd"
  cv_big_pls_core(data = data, nfold = nfold, nt = nt, keepX = keepX,
                  givefold = givefold, allCVcrit = allCVcrit,
                  times.auc = times.auc, times.prederr = times.prederr,
                  method = method, verbose = verbose, engine = engine, ...)
}

cv_big_pls_core <- function(data, nfold, nt, keepX, givefold, allCVcrit,
                            times.auc, times.prederr, method, verbose,
                            engine, ...) {
  if (missing(data) || !is.list(data)) {
    stop("`data` must be a list with elements x, time and status", call. = FALSE)
  }
  if (!all(c("x", "time", "status") %in% names(data))) {
    stop("`data` must contain x, time and status entries", call. = FALSE)
  }
  
  x <- data$x
  time <- as.numeric(data$time)
  status <- as.numeric(data$status)
  
  if (length(time) != length(status)) {
    stop("`time` and `status` must have the same length", call. = FALSE)
  }
  n <- length(time)
  if (n == 0L) {
    stop("no observations supplied", call. = FALSE)
  }
  
  nfold <- as.integer(nfold)
  if (length(nfold) != 1L || is.na(nfold) || nfold < 2L) {
    stop("`nfold` must be an integer >= 2", call. = FALSE)
  }
  
  if (is.null(nt)) {
    nt <- min(5L, ncol(x))
  }
  nt <- as.integer(nt)
  if (length(nt) != 1L || is.na(nt) || nt < 1L) {
    stop("`nt` must be a positive integer", call. = FALSE)
  }
  if (nt > ncol(x)) {
    stop("`nt` cannot exceed the number of predictors", call. = FALSE)
  }
  
  if (inherits(x, "big.matrix")) {
    if (!requireNamespace("bigmemory", quietly = TRUE)) {
      stop("Package 'bigmemory' is required to handle big.matrix inputs", call. = FALSE)
    }
    x_mat <- bigmemory::as.matrix(x[, , drop = FALSE])
  } else if (is.data.frame(x)) {
    x_mat <- as.matrix(x)
  } else if (is.matrix(x)) {
    x_mat <- x
  } else {
    x_mat <- as.matrix(x)
  }
  storage.mode(x_mat) <- "double"
  
  if (nrow(x_mat) != n) {
    stop("Number of rows in `x` must match the length of `time`", call. = FALSE)
  }
  
  if (is.null(times.auc)) {
    times.auc <- seq(0, max(time), length.out = 1000L)
  }
  if (is.null(times.prederr)) {
    if (length(times.auc) > 10L) {
      times.prederr <- times.auc[-seq(max(1L, length(times.auc) - 9L), length(times.auc))]
    } else {
      times.prederr <- times.auc
    }
  }
  
  if (missing(givefold)) {
    folds <- split(sample.int(n), rep(seq_len(nfold), length.out = n))
  } else {
    folds <- givefold
    if (!is.list(folds) || length(folds) != nfold) {
      stop("`givefold` must be a list with `nfold` elements", call. = FALSE)
    }
  }
  
  number_ind <- if (allCVcrit) 14L else 1L
  errormats <- replicate(number_ind, matrix(NA_real_, nrow = nt + 1L, ncol = nfold), simplify = FALSE)
  
  for (fold_id in seq_len(nfold)) {
    omit <- folds[[fold_id]]
    if (!is.integer(omit)) omit <- as.integer(omit)
    keep <- setdiff(seq_len(n), omit)
    n_train <- length(keep)
    if (nt > min(n_train, ncol(x_mat))) {
      stop("`nt` is too large for at least one training fold", call. = FALSE)
    }
    if (!requireNamespace("bigmemory", quietly = TRUE)) {
      stop("Package 'bigmemory' is required", call. = FALSE)
    }
    train_big <- bigmemory::as.big.matrix(x_mat[keep, , drop = FALSE])
    fit_args <- list(X = train_big, time = time[keep], status = status[keep],
                     ncomp = nt, keepX = keepX)
    if (engine == "pls") {
      fit <- do.call(big_pls_cox, c(fit_args, list(...)))
    } else if (engine == "gd") {
      fit <- do.call(big_pls_cox_gd, c(fit_args, list(...)))
    } else {
      stop("Unknown engine", call. = FALSE)
    }
    
    train_scores <- fit$scores
    if (!is.matrix(train_scores)) {
      train_scores <- as.matrix(train_scores)
    }
    test_scores <- predict(fit, newdata = x_mat[omit, , drop = FALSE], type = "components")
    if (!is.matrix(test_scores)) {
      test_scores <- as.matrix(test_scores)
    }
    if (allCVcrit) {
      full_scores <- predict(fit, newdata = x_mat, type = "components")
      if (!is.matrix(full_scores)) {
        full_scores <- as.matrix(full_scores)
      }
    } else {
      full_scores <- NULL
    }
    
    fold_metrics <- compute_fold_metrics(train_scores = train_scores,
                                         test_scores = test_scores,
                                         full_scores = full_scores,
                                         time = time, status = status,
                                         keep = keep, omit = omit,
                                         nt = nt,
                                         method = method,
                                         allCVcrit = allCVcrit,
                                         times.auc = times.auc,
                                         times.prederr = times.prederr)
    
    if (allCVcrit) {
      for (idx in seq_len(number_ind)) {
        errormats[[idx]][, fold_id] <- fold_metrics[[idx]]
      }
    } else {
      errormats[[1L]][, fold_id] <- fold_metrics[[1L]]
    }
    
    if (verbose) {
      message("Completed fold ", fold_id)
    }
  }
  
  if (allCVcrit) {
    titles_sign <- c(1, 1, rep(-1, 8), rep(1, 4))
    result <- list(nt = nt)
    for (idx in seq_len(number_ind)) {
      cv_mean <- apply(errormats[[idx]], 1L, mean, na.rm = TRUE)
      cv_se <- sqrt(apply(errormats[[idx]], 1L, var, na.rm = TRUE)) / nfold
      lambda <- getmin2(0:nt, titles_sign[idx] * cv_mean, cv_se)
      result[[paste0("cv.error", idx)]] <- cv_mean
      result[[paste0("cv.se", idx)]] <- cv_se
      result[[paste0("lambda.min", idx)]] <- lambda$lambda.min
      result[[paste0("lambda.1se", idx)]] <- lambda$lambda.1se
    }
    result$folds <- folds
    return(result)
  } else {
    cv_mean <- apply(errormats[[1L]], 1L, mean, na.rm = TRUE)
    cv_se <- sqrt(apply(errormats[[1L]], 1L, var, na.rm = TRUE)) / nfold
    lambda <- getmin2(0:nt, -cv_mean, cv_se)
    list(nt = nt, cv.error10 = cv_mean, cv.se10 = cv_se,
         folds = folds, lambda.min10 = lambda$lambda.min,
         lambda.1se10 = lambda$lambda.1se)
  }
}

compute_fold_metrics <- function(train_scores, test_scores, full_scores,
                                 time, status, keep, omit, nt,
                                 method, allCVcrit, times.auc, times.prederr) {
  test_size <- length(omit)
  n_train <- length(keep)
  pred_list <- vector("list", if (allCVcrit) 14L else 1L)
  for (idx in seq_along(pred_list)) {
    pred_list[[idx]] <- rep(NA_real_, nt + 1L)
  }
  
  time_train <- time[keep]
  status_train <- status[keep]
  time_test <- time[omit]
  status_test <- status[omit]
  
  has_event_train <- any(status_train == 1)
  has_event_test <- any(status_test == 1)
  
  zero_lp <- rep(0, n_train)
  zero_lp_test <- rep(0, test_size)
  if (allCVcrit) {
    pred_list[[1L]][1L] <- if (has_event_test) {
      logplik(matrix(0, nrow = test_size, ncol = 1L), time_test, status_test,
              matrix(0, nrow = 1L, ncol = 1L), method = method, return.all = FALSE)
    } else NA_real_
    pred_list[[2L]][1L] <- NA_real_
  }
  
  metric_values <- evaluate_metrics(lp_train = zero_lp, lp_test = zero_lp_test,
                                    time_train = time_train,
                                    status_train = status_train,
                                    time_test = time_test,
                                    status_test = status_test,
                                    method = method,
                                    allCVcrit = allCVcrit,
                                    times.auc = times.auc,
                                    times.prederr = times.prederr)
  if (allCVcrit) {
    pred_list[[3L]][1L] <- metric_values$AUC_CD
    pred_list[[4L]][1L] <- metric_values$AUC_hc
    pred_list[[5L]][1L] <- metric_values$AUC_sh
    pred_list[[6L]][1L] <- metric_values$AUC_Uno
    pred_list[[7L]][1L] <- metric_values$AUC_hz_train
    pred_list[[8L]][1L] <- metric_values$AUC_hz_test
    pred_list[[9L]][1L] <- metric_values$AUC_survROC_train
    pred_list[[10L]][1L] <- metric_values$AUC_survROC_test
    pred_list[[11L]][1L] <- metric_values$Brier_unw
    pred_list[[12L]][1L] <- metric_values$Schmid_unw
    pred_list[[13L]][1L] <- metric_values$Brier_w
    pred_list[[14L]][1L] <- metric_values$Schmid_w
  } else {
    pred_list[[1L]][1L] <- metric_values$AUC_survROC_test
  }
  
  for (comp in seq_len(nt)) {
    train_mat <- train_scores[, seq_len(comp), drop = FALSE]
    test_mat <- test_scores[, seq_len(comp), drop = FALSE]
    if (allCVcrit) {
      full_mat <- full_scores[, seq_len(comp), drop = FALSE]
    } else {
      full_mat <- NULL
    }
    
    control_fit <- survival::coxph.control()
    control_fit$iter.max <- max(50L, control_fit$iter.max)
    cox_fit <- tryCatch(
      survival::coxph.fit(x = train_mat,
                          y = cbind(time_train, status_train),
                          strata = rep(1, n_train),
                          offset = rep(0, n_train),
                          init = rep(0, comp),
                          control = control_fit,
                          weights = rep(1, n_train),
                          method = method,
                          rownames = as.character(seq_len(n_train))),
      error = function(e) NULL
    )
    if (is.null(cox_fit) || any(!is.finite(cox_fit$coefficients))) {
      next
    }
    coef_vec <- as.numeric(cox_fit$coefficients)
    
    lp_train <- drop(train_mat %*% coef_vec)
    lp_test <- drop(test_mat %*% coef_vec)
    
    if (allCVcrit) {
      pred_list[[1L]][comp + 1L] <- logplik(test_mat, time_test, status_test,
                                            matrix(coef_vec, ncol = 1L),
                                            method = method, return.all = FALSE)
      pred_list[[2L]][comp + 1L] <- logplik(full_mat, time, status,
                                            matrix(coef_vec, ncol = 1L),
                                            method = method, return.all = FALSE) -
        logplik(train_mat, time_train, status_train,
                matrix(coef_vec, ncol = 1L), method = method, return.all = FALSE)
    }
    
    metric_values <- evaluate_metrics(lp_train = lp_train, lp_test = lp_test,
                                      time_train = time_train,
                                      status_train = status_train,
                                      time_test = time_test,
                                      status_test = status_test,
                                      method = method,
                                      allCVcrit = allCVcrit,
                                      times.auc = times.auc,
                                      times.prederr = times.prederr)
    if (allCVcrit) {
      pred_list[[3L]][comp + 1L] <- metric_values$AUC_CD
      pred_list[[4L]][comp + 1L] <- metric_values$AUC_hc
      pred_list[[5L]][comp + 1L] <- metric_values$AUC_sh
      pred_list[[6L]][comp + 1L] <- metric_values$AUC_Uno
      pred_list[[7L]][comp + 1L] <- metric_values$AUC_hz_train
      pred_list[[8L]][comp + 1L] <- metric_values$AUC_hz_test
      pred_list[[9L]][comp + 1L] <- metric_values$AUC_survROC_train
      pred_list[[10L]][comp + 1L] <- metric_values$AUC_survROC_test
      pred_list[[11L]][comp + 1L] <- metric_values$Brier_unw
      pred_list[[12L]][comp + 1L] <- metric_values$Schmid_unw
      pred_list[[13L]][comp + 1L] <- metric_values$Brier_w
      pred_list[[14L]][comp + 1L] <- metric_values$Schmid_w
    } else {
      pred_list[[1L]][comp + 1L] <- metric_values$AUC_survROC_test
    }
  }
  
  if (allCVcrit) {
    pred_list[[1L]] <- ifelse(is.finite(pred_list[[1L]]), -pred_list[[1L]] / test_size, NA_real_)
    pred_list[[2L]] <- ifelse(is.finite(pred_list[[2L]]), -pred_list[[2L]] / test_size, NA_real_)
    for (idx in 3:10) {
      pred_list[[idx]] <- ifelse(is.finite(pred_list[[idx]]), pred_list[[idx]], NA_real_)
    }
    for (idx in 11:14) {
      pred_list[[idx]] <- ifelse(is.finite(pred_list[[idx]]), pred_list[[idx]], NA_real_)
    }
  } else {
    pred_list[[1L]] <- ifelse(is.finite(pred_list[[1L]]), pred_list[[1L]], NA_real_)
  }
  
  pred_list
}

evaluate_metrics <- function(lp_train, lp_test, time_train, status_train,
                             time_test, status_test, method, allCVcrit,
                             times.auc, times.prederr) {
  if (!any(status_train == 1) || !any(status_test == 1)) {
    return(list(AUC_CD = NA_real_, AUC_hc = NA_real_, AUC_sh = NA_real_,
                AUC_Uno = NA_real_, AUC_hz_train = NA_real_,
                AUC_hz_test = NA_real_, AUC_survROC_train = NA_real_,
                AUC_survROC_test = NA_real_, Brier_unw = NA_real_,
                Schmid_unw = NA_real_, Brier_w = NA_real_,
                Schmid_w = NA_real_))
  }
  train_df <- data.frame(time = time_train, status = status_train, Xlp = lp_train)
  test_df <- data.frame(time = time_test, status = status_test, Xlp = lp_test)
  Surv.rsp <- survival::Surv(time_train, status_train)
  Surv.rsp.new <- survival::Surv(time_test, status_test)
  
  train_fit <- survival::coxph(Surv(time, status) ~ Xlp, data = train_df,
                               method = method, iter.max = 0, init = 1,
                               x = TRUE, y = TRUE)
  lp <- as.numeric(predict(train_fit))
  lpnew <- as.numeric(predict(train_fit, newdata = test_df))
  
  if (allCVcrit) {
    AUCs <- getIndicCV(lp, lpnew, Surv.rsp, Surv.rsp.new,
                       times.auc = times.auc, times.prederr = times.prederr,
                       train.fit = train_fit, plot.it = FALSE)
    list(AUC_CD = AUCs$AUC_CD$iauc,
         AUC_hc = AUCs$AUC_hc$iauc,
         AUC_sh = AUCs$AUC_sh$iauc,
         AUC_Uno = AUCs$AUC_Uno$iauc,
         AUC_hz_train = AUCs$AUC_hz.train$iauc,
         AUC_hz_test = AUCs$AUC_hz.test$iauc,
         AUC_survROC_train = AUCs$AUC_survivalROC.train$iauc,
         AUC_survROC_test = AUCs$AUC_survivalROC.test$iauc,
         Brier_unw = AUCs$prederr$brier.unw$ierror,
         Schmid_unw = AUCs$prederr$robust.unw$ierror,
         Brier_w = AUCs$prederr$brier.w$ierror,
         Schmid_w = AUCs$prederr$robust.w$ierror)
  } else {
    AUCs <- getIndicCViAUCSurvROCTest(lp, lpnew, Surv.rsp, Surv.rsp.new,
                                      times.auc = times.auc,
                                      times.prederr = times.prederr,
                                      train.fit = train_fit, plot.it = FALSE)
    list(AUC_CD = NA_real_,
         AUC_hc = NA_real_,
         AUC_sh = NA_real_,
         AUC_Uno = NA_real_,
         AUC_hz_train = NA_real_,
         AUC_hz_test = NA_real_,
         AUC_survROC_train = NA_real_,
         AUC_survROC_test = AUCs$AUC_survivalROC.test$iauc,
         Brier_unw = NA_real_,
         Schmid_unw = NA_real_,
         Brier_w = NA_real_,
         Schmid_w = NA_real_)
  }
}
