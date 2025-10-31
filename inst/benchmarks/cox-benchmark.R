#!/usr/bin/env Rscript

# Benchmark coxgpls() against survival::coxph() on simulated data. Run from the
# package root with:
#
#   Rscript inst/benchmarks/cox-benchmark.R
#
# Results are stored under inst/benchmarks/results/ with a time-stamped filename.
# Customise the simulation through the options bigPLScox.benchmark.n,
# bigPLScox.benchmark.p, bigPLScox.benchmark.ncomp, and
# bigPLScox.benchmark.iterations.

suppressPackageStartupMessages({
  library(bigPLScox)
  library(survival)
  library(bench)
})

set.seed(2024)

n_obs <- as.integer(getOption("bigPLScox.benchmark.n", 2000))
n_pred <- as.integer(getOption("bigPLScox.benchmark.p", 50))
n_comp <- as.integer(getOption("bigPLScox.benchmark.ncomp", 5))
iterations <- as.integer(getOption("bigPLScox.benchmark.iterations", 50))

message("Simulating dataset with ", n_obs, " observations and ", n_pred, " predictors")

sim_design <- dataCox(
  n = n_obs,
  lambda = 2,
  rho = 1.5,
  x = matrix(rnorm(n_obs * n_pred), ncol = n_pred),
  beta = c(1, 3, rep(0, n_pred - 2)),
  cens.rate = 5
)

cox_data <- list(
  x = as.matrix(sim_design[, -(1:3)]),
  time = sim_design$time,
  status = sim_design$status
)

message("Running benchmarks with ", iterations, " iterations per estimator")

bench_res <- bench::mark(
  bigPLScox = coxgpls(
    cox_data$x,
    cox_data$time,
    cox_data$status,
    ncomp = n_comp
  ),
  survival = coxph(Surv(time, status) ~ x, data = cox_data, ties = "breslow"),
  iterations = iterations,
  check = FALSE
)

print(bench_res)

out_dir <- file.path("inst", "benchmarks", "results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_file <- file.path(out_dir, paste0("cox-benchmark-", timestamp, ".csv"))

message("Writing benchmark results to ", out_file)
write.csv(as.data.frame(bench_res[, 1:9]), file = out_file, row.names = FALSE)

message("Done")
