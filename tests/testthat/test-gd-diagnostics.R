test_that("GD diagnostics have correct length", {
  skip_if_not(requireNamespace("survival", quietly = TRUE))
  
  set.seed(123)
  n <- 200
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n)
  status <- rbinom(n, 1, 0.7)
  
  bm <- bigmemory::as.big.matrix(X)
  
  fit <- big_pls_cox_gd(
    X = bm,
    time = time,
    status = status,
    ncomp = 3,
    max_iter = 100
  )
  
  d <- gd_diagnostics(fit)
  expect_true(is.list(d))
  expect_true(all(c("iterations", "loglik", "step_sizes", "gradient_norm", "coef_trace", "eta_trace") %in% names(d)))
  
  # lengths should match the actual number of iterations
  iters <- length(d$iterations)
  expect_length(d$iterations, iters)
  expect_length(d$loglik, iters)
  expect_length(d$step_sizes, iters)
  expect_length(d$gradient_norm, iters)
  expect_length(d$coef_trace, iters)
  expect_length(d$eta_trace, iters)
})
