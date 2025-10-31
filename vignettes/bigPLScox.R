## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/bigmemory-",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

## ----simulate-bigmatrix-------------------------------------------------------
library(bigPLScox)
library(bigmemory)

set.seed(2024)
n_obs <- 5000
n_pred <- 100

X_dense <- matrix(rnorm(n_obs * n_pred), nrow = n_obs)
time <- rexp(n_obs, rate = 0.2)
status <- rbinom(n_obs, 1, 0.7)

big_dir <- tempfile("bigPLScox-")
dir.create(big_dir)
X_big <- filebacked.big.matrix(
  nrow = n_obs,
  ncol = n_pred,
  backingpath = big_dir,
  backingfile = "X.bin",
  descriptorfile = "X.desc",
  init = X_dense
)

## ----big-pls-cox--------------------------------------------------------------
fit_big <- big_pls_cox(
  X = X_big,
  time = time,
  status = status,
  ncomp = 5
)

head(fit_big$scores)
str(fit_big)

## ----big-pls-cox-gd-----------------------------------------------------------
fit_big_gd <- big_pls_cox_gd(
  X = X_big,
  time = time,
  status = status,
  ncomp = 5,
  max_iter = 100,
  tol = 1e-4
  )

head(fit_big$scores)
str(fit_big)

## ----big-cv, eval = FALSE-----------------------------------------------------
# set.seed(2024)
# data_big <- list(x = X_big, time = time, status = status)
# cv_big <- cv.coxgpls(
#   data_big,
#   nt = 5,
#   ncores = 1,
#   ind.block.x = c(10, 40)
# )
# cv_big$opt_nt

## ----big-timing---------------------------------------------------------------
if (requireNamespace("bench", quietly = TRUE)) {
  bench::mark(
    streaming = big_pls_cox(X_big, time, status, ncomp = 5, keepX = 0),
    gd = big_pls_cox_gd(X_big, time, status, ncomp = 5, max_iter = 150),
    iterations = 5,
    check = FALSE
  )
}

## ----big-deviance-------------------------------------------------------------
eta_big <- predict(fit_big, type = "link")
dr_cpp <- computeDR(time, status, engine = "cpp", eta = eta_big)
max(abs(dr_cpp - computeDR(time, status)))

## ----cleanup------------------------------------------------------------------
rm(X_big)
file.remove(file.path(big_dir, "X.bin"))
file.remove(file.path(big_dir, "X.desc"))
unlink(big_dir, recursive = TRUE)

