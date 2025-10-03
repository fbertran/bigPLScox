skip_on_os("windows")

library(bigmemory)

test_that("gradient descent returns consistent coefficients", {
  set.seed(42)
  n <- 40
  p <- 6
  base_matrix <- matrix(rnorm(n * p), nrow = n)
  X <- bigmemory::as.big.matrix(base_matrix)
  time <- rexp(n, rate = 0.2)
  status <- rbinom(n, 1, 0.75)
  
  fit1 <- big_pls_cox_gd(X, time, status, ncomp = 3, max_iter = 100, tol = 1e-6)
  fit2 <- big_pls_cox_gd(X, time, status, ncomp = 3, max_iter = 100, tol = 1e-6)
  
  X2 <- bigmemory::as.big.matrix(base_matrix)
  fit3 <- big_pls_cox_gd(X2, time, status, ncomp = 3, max_iter = 100, tol = 1e-6)
  
  expect_equal(length(fit1$coefficients), 3L)
  expect_equal(fit1$coefficients, fit2$coefficients)
  expect_equal(fit1$coefficients, fit3$coefficients)
  expect_true(is.numeric(fit1$loglik))
  expect_true(fit1$iterations >= 1)
})
