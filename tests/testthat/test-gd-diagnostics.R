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
})

test_that("gradient based methods give consistent latent Cox fits", {
  skip_if_not_installed("survival")
  skip_if_not_installed("bigmemory")
  
  set.seed(123)
  n <- 200
  p <- 40
  
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.7)
  
  X_big <- bigmemory::as.big.matrix(X)

  gd_bb <- big_pls_cox_gd(
    X        = X_big,
    time     = time,
    status   = status,
    ncomp    = 4,
    max_iter = 2000,
    method   = "bb"
  )
  
  gd_bfgs <- big_pls_cox_gd(
    X        = X_big,
    time     = time,
    status   = status,
    ncomp    = 4,
    max_iter = 2000,
    method   = "bfgs"
  )
  
  # same scores for all GD variants, up to sign and small numerical noise
  expect_equal(dim(gd_bb$scores), dim(gd_bfgs$scores))
  cors <- abs(cor(gd_bb$scores, gd_bfgs$scores))
  expect_true(all(diag(cors) > 0.99))
  
  # Cox refits on scores are close, even if raw GD coefficients differ
  lp_gd_bb  <- as.numeric(gd_bb$cox_fit$linear.predictors)
  lp_gd_bfgs <- as.numeric(gd_bfgs$cox_fit$linear.predictors)
  
  expect_true(cor(lp_gd_bfgs, lp_gd_bb) > 0.9)
})

test_that("predict.big_pls_cox_gd handles types and components correctly", {
  skip_if_not_installed("survival")
  skip_if_not_installed("bigmemory")
  
  set.seed(456)
  n <- 150
  p <- 30
  
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, rate = 0.2)
  status <- rbinom(n, 1, 0.6)
  
  X_big <- bigmemory::as.big.matrix(X)
  
  fit <- big_pls_cox_gd(
    X        = X_big,
    time     = time,
    status   = status,
    ncomp    = 3,
    method   = "bb",
    max_iter = 1000
  )
  
  # in sample scores from predict match stored scores
  scores_pred <- predict(fit, type = "components")
  expect_equal(dim(scores_pred), dim(fit$scores))
  expect_equal(scores_pred, fit$scores, tolerance = 1e-6)
  
  # use the Cox refit coefficients as reference
  beta_cox <- stats::coef(fit$cox_fit)
  lp_pred  <- predict(fit, type = "link")
  lp_ref   <- as.numeric(fit$scores %*% beta_cox)
  
  expect_equal(lp_pred, lp_ref, tolerance = 1e-6)
  
  # subset of components
  lp_12 <- predict(fit, type = "link", comps = 1:2, coef = beta_cox)
  expect_length(lp_12, n)
  # reference computation with only first two components
  lp_12_ref <- as.numeric(fit$scores[, 1:2, drop = FALSE] %*% beta_cox[1:2])
  expect_equal(lp_12, lp_12_ref, tolerance = 1e-6)
  
  # risk and response types are positive
  risk <- predict(fit, type = "risk")
  resp <- predict(fit, type = "response")
  expect_true(all(risk > 0))
  expect_true(all(resp > 0))
})

test_that("gradient diagnostics have consistent lengths", {
  skip_if_not_installed("survival")
  skip_if_not_installed("bigmemory")
  
  set.seed(789)
  n <- 120
  p <- 20
  
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, rate = 0.15)
  status <- rbinom(n, 1, 0.7)
  
  X_big <- bigmemory::as.big.matrix(X)
  
  fit <- big_pls_cox_gd(
    X        = X_big,
    time     = time,
    status   = status,
    ncomp    = 2,
    method   = "bb",
    max_iter = 200
  )
  
  lt <- fit$loglik_trace
  st <- fit$step_trace
  gnt <- fit$gradnorm_trace
  expect_true(fit$iterations >= 1L)
  
  same_len <- vapply(
    list(lt, gnt, st),
    function(x) length(x) == fit$iterations,
    logical(1)
  )
  expect_true(all(same_len))
})

