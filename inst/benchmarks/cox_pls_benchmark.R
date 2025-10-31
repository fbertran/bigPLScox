#!/usr/bin/env Rscript

# Benchmarking script for bigPLScox
#
# The goal is to contrast the native C++-accelerated implementation exposed via
# big_pls_cox() against a reference implementation from the plsRcox package on
# both in-memory and bigmemory-backed matrices. The script also reports the
# memory footprint and scaling behaviour when multiple cores are used.

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
n <- 800
p <- 60

X_dense <- matrix(rnorm(n * p), nrow = n)
time <- rexp(n, rate = 0.1)
status <- rbinom(n, 1, 0.6)

X_big <- bigmemory::as.big.matrix(X_dense)

run_big_pls_cox <- function(x_mat) {
  big_pls_cox(
    X = x_mat,
    time = time,
    status = status,
    ncomp = 5
  )
}

run_plsRcox <- function(x_mat) {
  df <- data.frame(time = time, status = status, x_mat)
  plsRcox::plsRcox(
    ~x_mat, 
    time=time,
    event=status,
    nt = 5,
    scaleX = TRUE
  )
}

cat("\n== Timing comparison on in-memory matrices ==\n")
bench_in_memory <- bench::mark(
  big_pls_cox = run_big_pls_cox(X_dense),
  plsRcox = run_plsRcox(X_dense),
  iterations = 5,
  check = FALSE
)
print(bench_in_memory[, c("expression", "min", "median", "itr/sec")])

cat("\n== Timing comparison on bigmemory-backed matrices ==\n")
bench_bigmemory <- bench::mark(
  big_pls_cox = run_big_pls_cox(X_big),
  iterations = 5,
  check = FALSE
)
print(bench_bigmemory[, c("expression", "min", "median", "itr/sec")])

# Highlight the use of the underlying C++ accelerators directly
cat("\n== Column summary statistics computed via C++ helper ==\n")
col_stats <- bigPLScox:::big_pls_cox_col_stats_cpp(methods::slot(X_big, "address"))
print(str(col_stats))

cat("\n== Memory footprint (MB) ==\n")
print(rbind(
  dense = bench::as_bench_bytes(object.size(X_dense)),
  bigmatrix = bench::as_bench_bytes(bigmemory::GetMatrixSize(X_big) + object.size(X_big))
))

if (requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)) {
  cat("\n== Parallel scaling experiment ==\n")
  cores_to_try <- c(1, 2, parallel::detectCores(logical = FALSE))
  cores_to_try <- unique(pmax(1, cores_to_try))
  res <- lapply(cores_to_try, function(k) {
    doParallel::registerDoParallel(cores = k)
    tictoc <- system.time(run_big_pls_cox(X_big))
    foreach::registerDoSEQ()
    data.frame(cores = k, elapsed = tictoc[["elapsed"]])
  })
  scaling <- do.call(rbind, res)
  print(scaling)
} else {
  cat("\nParallel packages not installed; skipping scaling experiment.\n")
}

cat("\nBenchmark complete.\n")
