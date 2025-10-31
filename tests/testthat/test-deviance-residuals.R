skip_if_not_installed("survival")

test_that("C++ deviance residuals match survival implementation", {
  set.seed(123)
  time <- rexp(100, rate = 0.2)
  status <- rbinom(100, 1, 0.6)
  dr_surv <- residuals(survival::coxph(survival::Surv(time, status) ~ 1), type = "deviance")
  dr_surv_cdr <- suppressWarnings(bigPLScox::computeDR(time, status, engine = "survival"))
  dr_cpp_cdr <- bigPLScox::computeDR(time, status, engine = "cpp", eta=rep(0 ,100))
  dr_qcpp_cdr <- bigPLScox::computeDR(time, status, engine = "qcpp")
  expect_equal(unname(dr_surv_cdr), unname(dr_surv), tolerance = 1e-6)
  expect_equal(as.numeric(unname(dr_cpp_cdr)), unname(dr_surv), tolerance = 1e-6)
  expect_equal(as.numeric(unname(dr_qcpp_cdr)), unname(dr_surv), tolerance = 1e-6)
})

test_that("computeDR uses the C++ backend in the simple case", {
  set.seed(42)
  time <- rexp(50)
  status <- rbinom(50, 1, 0.7)
  res <- computeDR(time, status, scaleY = FALSE)
  dr_qcpp <- cox_deviance_details(time, status)
  expect_equal(unname(res), unname(dr_qcpp$deviance), tolerance = 1e-6)
})

test_that("deviance residuals from big.matrix agree with base implementation", {
  skip_if_not_installed("bigmemory")
  set.seed(321)
  time <- rexp(40)
  status <- rbinom(40, 1, 0.5)
  X <- cbind(time, status)
  bm <- bigmemory::as.big.matrix(X)
  dr_big <- cox_deviance_residuals_big(bm, time_col = 1, status_col = 2)
  dr_cpp <- cox_deviance_residuals(time, status)
  expect_equal(dr_big, dr_cpp, tolerance = 1e-6)
})

test_that("partial deviance from big matrices is consistent", {
  skip_if_not_installed("bigmemory")
  set.seed(999)
  X <- matrix(rnorm(60), nrow = 20)
  coef <- rnorm(ncol(X))
  time <- rexp(20)
  status <- rbinom(20, 1, 0.4)
  bm <- bigmemory::as.big.matrix(X)
  stats_big <- cox_partial_deviance_big(bm, coef, time, status)
  lp <- as.numeric(X %*% coef)
  unique_times <- sort(unique(time), decreasing = TRUE)
  loglik_ref <- 0
  for (t in unique_times) {
    idx_time <- which(time == t)
    event_idx <- idx_time[status[idx_time] > 0]
    if (length(event_idx) == 0) {
      next
    }
    risk_set <- which(time >= t)
    loglik_ref <- loglik_ref + sum(lp[event_idx]) - length(event_idx) * log(sum(exp(lp[risk_set])))
  }
  expect_length(stats_big$linear_predictor, nrow(X))
  expect_equal(stats_big$loglik, loglik_ref, tolerance = 1e-6)
})

