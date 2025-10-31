## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/overview-",
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

train_idx <- seq_len(80)
Y_train <- micro.censure$survyear[train_idx]
C_train <- micro.censure$DC[train_idx]
X_train <- Xmicro.censure_compl_imp[train_idx, -40]

## ----fit-coxgpls--------------------------------------------------------------
fit <- coxgpls(
  X_train,
  Y_train,
  C_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 15)
)
fit

## ----cv-coxgpls---------------------------------------------------------------
set.seed(123)
cv_res <- cv.coxgpls(
  list(x = X_train, time = Y_train, status = C_train),
  nt = 10,
  ind.block.x = c(3, 10, 15)
)
cv_res

## ----alternative-estimators---------------------------------------------------
dr_fit <- coxgplsDR(
  X_train,
  Y_train,
  C_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 15)
)
dr_fit

## ----bigmemory-example--------------------------------------------------------
X_big <- bigmemory::as.big.matrix(X_train)
big_fit <- big_pls_cox(
  X_big,
  time = Y_train,
  status = C_train,
  ncomp = 6
)
big_fit

