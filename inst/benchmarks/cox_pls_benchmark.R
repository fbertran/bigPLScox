#!/usr/bin/env Rscript

# Benchmark big_pls_cox() against plsRcox::plsRcox() on dense and bigmemory
# matrices. Run from the package root with:
#
#   Rscript inst/benchmarks/cox_pls_benchmark.R
#
# Results are stored under inst/benchmarks/results/ with time-stamped filenames.
# Customise the problem size through the options bigPLScox.benchmark.n,
# bigPLScox.benchmark.p, and bigPLScox.benchmark.iterations.

suppressPackageStartupMessages({
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("Package 'bench' is required for benchmarking. Install it with install.packages('bench').")
  }
  if (!requireNamespace("plsRcox", quietly = TRUE)) {
    stop("Package 'plsRcox' is required for the baseline comparison. Install it with install.packages('plsRcox').")
  }
  library(bigPLScox)
  library(bigmemory)
  library(bench)
})

set.seed(20240513)

n_obs <- as.integer(getOption("bigPLScox.benchmark.n", 800))
n_pred <- as.integer(getOption("bigPLScox.benchmark.p", 60))
n_comp <- as.integer(getOption("bigPLScox.benchmark.ncomp", 5))
iterations <- as.integer(getOption("bigPLScox.benchmark.iterations", 5))

message("Simulating dataset with ", n_obs, " observations and ", n_pred, " predictors")

X_dense <- matrix(rnorm(n_obs * n_pred), nrow = n_obs)
time <- rexp(n_obs, rate = 0.1)
status <- rbinom(n_obs, 1, 0.6)

X_big <- bigmemory::as.big.matrix(X_dense)

run_big_pls_cox <- function(x_mat) {
  big_pls_cox(
    X = x_mat,
    time = time,
    status = status,
    ncomp = n_comp
  )
}

run_plsRcox <- function(x_mat) {
  df <- data.frame(time = time, status = status, x_mat)
  plsRcox::plsRcox(
    ~ x_mat,
    time = time,
    event = status,
    nt = n_comp,
    scaleX = TRUE
  )
}

message("== Timing comparison on in-memory matrices ==")
bench_in_memory <- bench::mark(
  big_pls_cox = run_big_pls_cox(X_dense),
  plsRcox = run_plsRcox(X_dense),
  iterations = iterations,
  check = FALSE
)
print(bench_in_memory[, c("expression", "min", "median", "itr/sec")])

message("== Timing comparison on bigmemory-backed matrices ==")
bench_bigmemory <- bench::mark(
  big_pls_cox = run_big_pls_cox(X_big),
  iterations = iterations,
  check = FALSE
)
print(bench_bigmemory[, c("expression", "min", "median", "itr/sec")])

message("== Column summary statistics computed via C++ helper ==")
col_stats <- bigPLScox:::big_pls_cox_col_stats_cpp(methods::slot(X_big, "address"))
str(col_stats)

message("== Memory footprint (MB) ==")
print(rbind(
  dense = bench::as_bench_bytes(object.size(X_dense)),
  bigmatrix = bench::as_bench_bytes(bigmemory::GetMatrixSize(X_big) + object.size(X_big))
))

if (requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)) {
  message("== Parallel scaling experiment ==")
  cores_to_try <- c(1, 2, parallel::detectCores(logical = FALSE))
  cores_to_try <- unique(pmax(1, cores_to_try))
  res <- lapply(cores_to_try, function(k) {
    doParallel::registerDoParallel(cores = k)
    elapsed <- system.time(run_big_pls_cox(X_big))["elapsed"]
    foreach::registerDoSEQ()
    data.frame(cores = k, elapsed = elapsed)
  })
  scaling <- do.call(rbind, res)
  print(scaling)
} else {
  message("Parallel packages not installed; skipping scaling experiment.")
}

out_dir <- file.path("inst", "benchmarks", "results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
readr::write_csv(bench_in_memory, file.path(out_dir, paste0("cox_pls-benchmark-dense-", timestamp, ".csv")))
readr::write_csv(bench_bigmemory, file.path(out_dir, paste0("cox_pls-benchmark-bigmemory-", timestamp, ".csv")))

message("Benchmark complete. Results written to ", out_dir)
