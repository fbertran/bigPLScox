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

## ----packages-----------------------------------------------------------------
library(bigPLScox)
library(survival)
library(bench)

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

## ----run-benchmarks-----------------------------------------------------------
bench_res <- bench::mark(
  bigPLScox = coxgpls(
    cox_data$x,
    cox_data$time,
    cox_data$status,
    ncomp = 5,
    ind.block.x = c(3, 10)
  ),
  survival = coxph(Surv(cox_data$time, cox_data$status) ~ cox_data$x, ties = "breslow"),
  iterations = 100,
  check = FALSE
)
bench_res

## ----bench-plot---------------------------------------------------------------
plot(bench_res, type = "jitter")

## ----export-benchmark, eval = FALSE-------------------------------------------
# if (!dir.exists("inst/benchmarks/results")) {
#   dir.create("inst/benchmarks/results", recursive = TRUE)
# }
# write.csv(bench_res, file = "inst/benchmarks/results/benchmarking-demo.csv", row.names = FALSE)

