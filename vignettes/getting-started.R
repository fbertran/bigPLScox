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

Y_all <- micro.censure$survyear[1:80]
status_all <- micro.censure$DC[1:80]
X_all <- apply(
  as.matrix(Xmicro.censure_compl_imp),
  MARGIN = 2,
  FUN = as.numeric
)[1:80, ]

set.seed(2024)
train_id <- 1:60
test_id <- 61:80

X_train <- X_all[train_id, ]
X_test <- X_all[test_id, ]
Y_train <- Y_all[train_id]
Y_test <- Y_all[test_id]
status_train <- status_all[train_id]
status_test <- status_all[test_id]

## ----original-design----------------------------------------------------------
X_train_raw <- Xmicro.censure_compl_imp[train_id, ]
X_test_raw <- Xmicro.censure_compl_imp[test_id, ]

## ----deviance-residuals-------------------------------------------------------
residuals_overview <- computeDR(Y_train, status_train, plot = TRUE)
eta_null <- rep(0, length(Y_train))
head(residuals_overview)

if (requireNamespace("bench", quietly = TRUE)) {
  benchmark_dr <- bench::mark(
    survival = computeDR(Y_train, status_train, engine = "survival"),
    cpp = computeDR(Y_train, status_train, engine = "cpp", eta = eta_null),
    iterations = 10,
    check = FALSE
  )
  benchmark_dr
}

all.equal(
  as.numeric(computeDR(Y_train, status_train, engine = "survival")),
  as.numeric(computeDR(Y_train, status_train, engine = "cpp", eta = eta_null)),
  tolerance = 1e-7
)

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

