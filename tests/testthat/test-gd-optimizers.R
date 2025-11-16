# tests/testthat/test-gd-optimizers.R

test_that("GD diagnostics have consistent lengths", {
  skip_if_not_installed("survival")
  skip_if_not_installed("bigmemory")
  
  set.seed(789)
  
  n <- 60
  p <- 15
  X <- matrix(rnorm(n * p), n, p)
  
  time   <- rexp(n)
  status <- rbinom(n, 1, 0.8)
  
  fit <- big_pls_cox_gd(
    X             = X,
    time          = time,
    status        = status,
    ncomp         = 2L,
    max_iter      = 500L,
    tol           = 1e-7,
    learning_rate = 0.05,
    method        = "bb",
    diag          = TRUE
  )
  
  iter <- fit$iterations
  
  expect_true(is.numeric(iter))
  expect_true(iter >= 1)
  
  expect_equal(length(fit$loglik_trace),   iter)
  expect_equal(length(fit$step_trace),     iter)
  expect_equal(length(fit$gradnorm_trace), iter)
  
  expect_true(all(is.finite(fit$loglik_trace)))
  expect_true(all(is.finite(fit$gradnorm_trace)))
})
