# Benchmark script for bigPLScox
#
# Run from the package root with:
#   Rscript inst/benchmarks/cox-benchmark.R
# The script records timings comparing bigPLScox to survival::coxph and writes the
# results to inst/benchmarks/results/.

suppressPackageStartupMessages({
  library(bigPLScox)
  library(survival)
  library(bench)
})

set.seed(2024)
n_obs <- as.integer(getOption("bigPLScox.benchmark.n", 2000))
p <- as.integer(getOption("bigPLScox.benchmark.p", 50))
n_comp <- as.integer(getOption("bigPLScox.benchmark.ncomp", 5))

message("Simulating dataset with ", n_obs, " observations and ", p, " predictors")

sim_design <- dataCox(
  n = n_obs,
  lambda = 2,
  rho = 1.5,
  x = matrix(rnorm(n_obs * p), ncol = p),
  beta = c(1, 3, rep(0, p - 2)),
  cens.rate = 5
)

cox_data <- list(
  x = as.matrix(sim_design[,-(1:3)]),
  time = sim_design$time,
  status = sim_design$status
)

bench_res <- bench::mark(
  bigPLScox = coxgpls(
    cox_data$x,
    cox_data$time,
    cox_data$status,
    ncomp = n_comp
  ),
  survival = coxph(Surv(time, status) ~ x, data=cox_data, ties = "breslow"),
  iterations = 50,
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
write.csv(as.data.frame(bench_res[,1:9]), file = out_file, row.names = FALSE)

message("Done")
