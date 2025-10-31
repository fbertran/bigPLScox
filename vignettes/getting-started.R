## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4,
  dpi = 150
)

## ----load-data----------------------------------------------------------------
library(bigPLScox)

data(micro.censure)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]

data(Xmicro.censure_compl_imp)
X_train_micro <- apply(
  as.matrix(Xmicro.censure_compl_imp),
  MARGIN = 2,
  FUN = as.numeric
)[1:80, ]
X_train_micro_df <- data.frame(X_train_micro)

## ----original-design----------------------------------------------------------
X_train_micro_orig <- Xmicro.censure_compl_imp[1:80, ]
X_train_micro_orig_df <- data.frame(X_train_micro_orig)

## ----deviance-residuals-------------------------------------------------------
dev_res <- computeDR(Y_train_micro, C_train_micro, plot = FALSE)
head(dev_res)

## ----model-matrix-------------------------------------------------------------
coxgpls(~ ., Y_train_micro, C_train_micro,
        ncomp = 6,
        trace = TRUE,
        model_matrix = TRUE,
        dataXplan = X_train_micro_orig_df,
        ind.block.x = c(3, 10, 20))[1:10, 1:6]

## ----coxgpls-fits-------------------------------------------------------------
cox_gpls_fit <- coxgpls(
  X_train_micro,
  Y_train_micro,
  C_train_micro,
  ncomp = 6,
  ind.block.x = c(3, 10, 15)
)
cox_gpls_fit

## ----coxgpls-formula----------------------------------------------------------
cox_gpls_fit_formula <- coxgpls(
  ~ X_train_micro,
  Y_train_micro,
  C_train_micro,
  ncomp = 6,
  ind.block.x = c(3, 10, 15)
)
cox_gpls_fit_formula

## ----cv-coxgpls---------------------------------------------------------------
set.seed(123456)
cv_coxgpls_res <- cv.coxgpls(
  list(x = X_train_micro, time = Y_train_micro, status = C_train_micro),
  nt = 10,
  ind.block.x = c(3, 10, 15)
)
cv_coxgpls_res

## ----coxgplsdr----------------------------------------------------------------
cox_gplsDR_fit <- coxgplsDR(
  X_train_micro,
  Y_train_micro,
  C_train_micro,
  ncomp = 6,
  ind.block.x = c(3, 10, 15)
)
cox_gplsDR_fit

## ----cv-coxgplsdr-------------------------------------------------------------
set.seed(123456)
cv_coxgplsDR_res <- cv.coxgplsDR(
  list(x = X_train_micro, time = Y_train_micro, status = C_train_micro),
  nt = 10,
  ind.block.x = c(3, 10, 15)
)
cv_coxgplsDR_res

## ----coxdkgplsdr--------------------------------------------------------------
cox_DKsplsDR_fit <- coxDKgplsDR(
  X_train_micro,
  Y_train_micro,
  C_train_micro,
  ncomp = 6,
  validation = "CV",
  ind.block.x = c(3, 10, 15),
  verbose = FALSE
)
cox_DKsplsDR_fit

## ----cv-coxdkgplsdr-----------------------------------------------------------
set.seed(123456)
cv_coxDKgplsDR_res <- cv.coxDKgplsDR(
  list(x = X_train_micro, time = Y_train_micro, status = C_train_micro),
  nt = 10,
  ind.block.x = c(3, 10, 15)
)
cv_coxDKgplsDR_res

