## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/getting-started-",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

## ----load-data----------------------------------------------------------------
library(bigPLScox)

data(micro.censure)
data(Xmicro.censure_compl_imp)

Y_train <- micro.censure$survyear[1:80]
status_train <- micro.censure$DC[1:80]
X_train <- apply(
  as.matrix(Xmicro.censure_compl_imp),
  MARGIN = 2,
  FUN = as.numeric
)[1:80, ]

## ----original-design----------------------------------------------------------
X_train_raw <- Xmicro.censure_compl_imp[1:80, ]

## ----deviance-residuals-------------------------------------------------------
residuals_overview <- computeDR(Y_train, status_train, plot = TRUE)
head(residuals_overview)

## ----fit-coxgpls--------------------------------------------------------------
set.seed(123)
cox_pls_fit <- coxgpls(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 6,
    ind.block.x = c(3, 10, 20)
)
cox_pls_fit

## ----fit-formula--------------------------------------------------------------
cox_pls_fit_formula <- coxgpls(
  ~ ., Y_train, status_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 20),
  dataXplan = data.frame(X_train_raw)
)
cox_pls_fit_formula

## ----cv-coxgpls---------------------------------------------------------------
set.seed(123456)
cv_results <- suppressWarnings(cv.coxgpls(
  list(x = X_train, time = Y_train, status = status_train),
  nt = 6,
  ind.block.x = c(3, 10, 20)
))
cv_results

## ----fit-coxgplsdr------------------------------------------------------------
set.seed(123456)
cox_pls_dr <- coxgplsDR(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = cv_results$nt,
  ind.block.x = c(3, 10, 20)
)
cox_pls_dr

## ----dk-splines---------------------------------------------------------------
cox_DKsplsDR_fit <- coxDKgplsDR(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 6,
  validation = "CV",
  ind.block.x = c(3, 10, 20),
  verbose = FALSE
)
cox_DKsplsDR_fit

## ----cv-dk-splines------------------------------------------------------------
set.seed(123456)
cv_coxDKgplsDR_res <- suppressWarnings(cv.coxDKgplsDR(
  list(x = X_train, time = Y_train, status = status_train),
  nt = 6,
  ind.block.x = c(3, 10, 20)
))
cv_coxDKgplsDR_res

