#!/usr/bin/env Rscript

# Benchmark variants of the bigPLScox solvers against each other and the
# baseline survival::coxph() implementation. Run from the package root with:
#
#   Rscript inst/benchmarks/benchmark_bigPLScox.R
#
# Results are stored under inst/benchmarks/results/ with a time-stamped
# filename. Customise the simulation size via the options
# bigPLScox.benchmark.n, bigPLScox.benchmark.p, and
# bigPLScox.benchmark.ncomp.

suppressPackageStartupMessages({
  library(bench)
  library(bigPLScox)
  library(survival)  
})

set.seed(2024)

n_obs <- as.integer(getOption("bigPLScox.benchmark.n", 20000))
n_pred <- as.integer(getOption("bigPLScox.benchmark.p", 20))
n_comp <- as.integer(getOption("bigPLScox.benchmark.ncomp", 5))
iterations <- as.integer(getOption("bigPLScox.benchmark.iterations", 25))

message("Simulating dataset with ", n_obs, " observations and ", n_pred, " predictors")

sim_data <- dataCox(
  n = n_obs,
  lambda = 2,
  rho = 1.5,
  x = matrix(sample(0:1, size = n_obs * n_pred, replace = TRUE), ncol = n_pred),
  beta = c(1, 0.5, -0.75, 0, 1.25, -1, 0, 0.5, -0.5, 0, runif(max(0, n_pred - 10), -1, 1)),
  cens.rate = 3
)

x_matrix <- bigmemory::as.big.matrix(sim_data[, -(1:3)])

design <- list(
  x = x_matrix,
  time = sim_data$time,
  status = sim_data$status
)

data_frame <- list(
  x = as.matrix(sim_data[, -(1:3)]),
  time = sim_data$time,
  status = sim_data$status
)

message("Running benchmarks with ", iterations, " iterations per estimator")

res <- bench::mark(
  iterations = iterations,
  check = FALSE,
  survival = {
    fit <- coxph(Surv(time, status) ~ x, data = data_frame, ties = "breslow")
    invisible(fit)
  },
  coxgpls = {
    fit <- coxgpls(
      Xplan = data_frame$x,
      time = data_frame$time,
      status = data_frame$status,
      ncomp = n_comp
      )
    invisible(fit)
  },
  big_pls_cox = {
    fit <- big_pls_cox(
      design$x,
      design$time,
      design$status,
      ncomp = n_comp
    )
    invisible(fit)
  },
  big_pls_cox_gd = {
    fit <- big_pls_cox_gd(
      design$x,
      design$time,
      design$status,
      ncomp = n_comp
    )
    invisible(fit)
  },
  min_time = 0.5
)

print(res)

summary_res <- summary(res)
summary_res$speedup_vs_survival <- as.numeric(summary_res$median[1] / summary_res$median)
summary_res$speedup_vs_coxgpls <- as.numeric(summary_res$median[2] / summary_res$median)

ggplot2::autoplot(res)
ggplot2::autoplot(res, "jitter")
plot(res$mem_alloc~factor(attr(res$expression,"description")),xlab="Algorithm",ylab="Memory")

print(summary_res[, c("expression", "median", "itr/sec", "mem_alloc", "gc/sec", "n_itr", "n_gc", "speedup_vs_survival", "speedup_vs_coxgpls")])

out_dir <- file.path("inst", "benchmarks", "results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

readr::write_csv(
  x = summary_res,
  file = file.path(out_dir, sprintf("benchmark_bigPLScox-%s.csv", format(Sys.time(), "%Y%m%d-%H%M%S")))
)

message("Benchmark summary written to ", out_dir)
