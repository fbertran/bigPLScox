## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/benchmarking-",
  fig.width = 7,
  fig.height = 5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

## ----packages, include=FALSE--------------------------------------------------
needed <- c("plsRcox","boot","survival","glmnet","bigPLScox")
has    <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)

if (!all(has)) {
  missing <- paste(needed[!has], collapse = ", ")
  knitr::opts_chunk$set(eval = FALSE)
  message("Note: skipping code execution because these packages are missing: ", missing)
}

## ----packages_lib-------------------------------------------------------------
library(bigPLScox)
library(survival)
library(bench)
library(bigmemory)
library(plsRcox)

## ----simulate-data------------------------------------------------------------
set.seed(2024)
sim_design <- dataCox(
  n = 2000,
  lambda = 2,
  rho = 1.5,
  x = matrix(rnorm(2000 * 50), ncol = 50),
  beta = c(1, 3, rep(0, 48)),
  cens.rate = 5
)

cox_data <- list(
  x = as.matrix(sim_design[, -(1:3)]),
  time = sim_design$time,
  status = sim_design$status
)
X_big <- bigmemory::as.big.matrix(cox_data$x)

## ----run-benchmarks-----------------------------------------------------------
bench_res <- bench::mark(
  coxgpls = coxgpls(
    cox_data$x,
    cox_data$time,
    cox_data$status,
    ncomp = 5,
    ind.block.x = c(3, 10)
  ),
  big_pls = big_pls_cox(X_big, cox_data$time, cox_data$status, ncomp = 5),
  big_pls_gd = big_pls_cox_gd(X_big, cox_data$time, cox_data$status, ncomp = 5, max_iter = 100),
  survival = coxph(Surv(cox_data$time, cox_data$status) ~ cox_data$x, ties = "breslow"),
  iterations = 100,
  check = FALSE
)
bench_res

bench_summary <- bench_res[, c("expression", "median", "itr/sec")]
bench_summary

## ----bench-plot---------------------------------------------------------------
plot(bench_res, type = "jitter")

## ----accuracy-shared----------------------------------------------------------
set.seed(4242)
split_id <- sample(c("train", "test"), size = nrow(cox_data$x), replace = TRUE, prob = c(0.7, 0.3))
train_idx <- split_id == "train"
test_idx <- !train_idx

train_x <- cox_data$x[train_idx, , drop = FALSE]
test_x <- cox_data$x[test_idx, , drop = FALSE]
train_time <- cox_data$time[train_idx]
test_time <- cox_data$time[test_idx]
train_status <- cox_data$status[train_idx]
test_status <- cox_data$status[test_idx]

big_fit <- big_pls_cox(bigmemory::as.big.matrix(train_x), train_time, train_status, ncomp = 5)
gd_fit <- big_pls_cox_gd(bigmemory::as.big.matrix(train_x), train_time, train_status, ncomp = 5, max_iter = 120)
pls_fit <- plsRcox::plsRcox(train_x, train_time, train_status, nt = 5, verbose = FALSE)

lp_big <- predict(big_fit, newdata = test_x, type = "link")
lp_gd <- predict(gd_fit, newdata = test_x, type = "link")
lp_pls <- predict(pls_fit, newdata = test_x)

surv_test <- survival::Surv(test_time, test_status)
acc_results <- data.frame(
  method = c("big_pls_cox", "big_pls_cox_gd", "plsRcox"),
  c_index = c(
    survival::concordance(surv_test ~ lp_big)$concordance,
    survival::concordance(surv_test ~ lp_gd)$concordance,
    survival::concordance(surv_test ~ lp_pls)$concordance
  )
)
acc_results

## ----export-benchmark, eval = FALSE-------------------------------------------
# if (!dir.exists("inst/benchmarks/results")) {
#   dir.create("inst/benchmarks/results", recursive = TRUE)
# }
# write.csv(bench_res[,1:9], file = "inst/benchmarks/results/benchmarking-demo.csv", row.names = FALSE)

