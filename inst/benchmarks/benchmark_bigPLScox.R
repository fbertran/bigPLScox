# Benchmarks for bigPLScox
#
# These scripts are intended to be run manually via `bench::press()` or
# `bench::mark()` to compare algorithmic variants provided by the package.
# They are not executed automatically during package checks.

suppressPackageStartupMessages({
  library(bench)
  library(bigPLScox)
})

set.seed(2024)

# Generate a moderately sized simulated dataset -----------------------------
sim_data <- dataCox(
  n = 20000,
  lambda = 2,
  rho = 1.5,
  x = matrix(sample(0:1, size = 20000 * 20, replace = TRUE), ncol = 20),
  beta = c(1, 0.5, -0.75, 0, 1.25, -1, 0, 0.5, -0.5, 0, runif(10,-1,1)),
  cens.rate = 3
)

x_matrix <- bigmemory::as.big.matrix(sim_data[,-c(1,2,3)])

cox_data <- list(
  x = as.matrix(sim_data[,-c(1,2,3)]),
  time = sim_data$time,
  status = sim_data$status
)

design <- list(
  x = x_matrix,
  time = sim_data$time,
  status = sim_data$status
)

# Compare the different Cox-PLS solvers --------------------------------------
res <- bench::mark(
  iterations = 25,
  check = FALSE,
  survival = {
    fit <- coxph(Surv(time, status) ~ x, data = cox_data, ties = "breslow")
    invisible(fit)
  },
  coxgpls = {
    fit <- coxgpls(
      Xplan=cox_data$x,
    time=cox_data$time,
    event=cox_data$status,
    ncomp = 5
    )
    invisible(fit)
  },  
  big_pls_cox = {
    fit <- big_pls_cox(
      design$x,
      design$time,
      design$status,
      ncomp = 5
    )
    invisible(fit)
  },
  big_pls_cox_gd = {
    fit <- big_pls_cox_gd(
      design$x,
      design$time,
      design$status,
      ncomp = 5
    )
    invisible(fit)
  },
  min_time = 0.5
)

print(res)

# Summarise the speedup of the deviance residual variant ---------------------
summary_res <- summary(res)
summary_res$speedup_vs_survival <- as.numeric(summary_res$median[1] /summary_res$median)
summary_res$speedup_vs_coxgpls <- as.numeric(summary_res$median[2] /summary_res$median)

ggplot2::autoplot(res)
ggplot2::autoplot(res, "ridge")
plot(res$mem_alloc~factor(attr(res$expression,"description")),xlab="Algorithm",ylab="Memory")

print(summary_res[, c("expression", "median", "itr/sec", "mem_alloc", "gc/sec", "n_itr", "n_gc", "speedup_vs_survival", "speedup_vs_coxgpls")])

# Save the benchmark results for later inspection ---------------------------
if (!dir.exists("inst/benchmarks/results")) {
  dir.create("inst/benchmarks/results", recursive = TRUE)
}

readr::write_csv(
  x = summary_res,
  file = file.path("inst", "benchmarks", "results", sprintf("benchmark-%s.csv", format(Sys.Date(), "%Y%m%d")))
)
