test_that("function returns latent components and consistent fits", {
  set.seed(42)
  n <- 40
  p <- 6
  base_matrix <- matrix(rnorm(n * p), nrow = n)
  X <- bigmemory::as.big.matrix(base_matrix)
  time <- rexp(n, rate = 0.2)
  status <- rbinom(n, 1, 0.75)
  
  fit1 <- big_pls_cox(X, time, status, ncomp = 3)
  fit2 <- big_pls_cox(X, time, status, ncomp = 3)
  
  expect_equal(length(fit1$cox_fit$coefficients), 3L)
  expect_equal(dim(fit1$scores), c(n, 3))
  expect_equal(dim(fit1$loadings), c(p, 3))
  expect_equal(dim(fit1$weights), c(p, 3))
  expect_equal(length(fit1$center), p)
  expect_equal(length(fit1$scale), p)
  expect_equal(fit1$cox_fit$coefficients, fit2$cox_fit$coefficients)
  expect_equal(fit1$cox_fit$loglik, fit2$cox_fit$loglik)
  expect_true(is.numeric(fit1$cox_fit$loglik))
  expect_equal(fit1$scores, fit2$scores)
})

test_that("reordering columns leaves the fit unchanged", {
  set.seed(101)
  n <- 35
  p <- 5
  base_matrix <- matrix(rnorm(n * p), nrow = n)
  time <- rexp(n, rate = 0.3)
  status <- rbinom(n, 1, 0.6)
  
  X <- bigmemory::as.big.matrix(base_matrix)
  fit_orig <- big_pls_cox(X, time, status, ncomp = 3)
  
  perm <- sample.int(p)
  X_perm <- bigmemory::as.big.matrix(base_matrix[, perm])
  fit_perm <- big_pls_cox(X_perm, time, status, ncomp = 3)
  
  expect_equal(fit_orig$cox_fit$coefficients, fit_perm$cox_fit$coefficients, tolerance = 1e-6)
  expect_equal(fit_orig$loglik, fit_perm$loglik, tolerance = 1e-6)
  expect_equal(fit_orig$scores, fit_perm$scores, tolerance = 1e-6)
  
  inv_perm <- match(seq_len(p), perm)
  expect_equal(fit_orig$loadings, fit_perm$loadings[inv_perm, , drop = FALSE],
               tolerance = 1e-6)
  expect_equal(fit_orig$weights, fit_perm$weights[inv_perm, , drop = FALSE],
               tolerance = 1e-6)
  expect_equal(fit_orig$center, fit_perm$center[inv_perm], tolerance = 1e-6)
  expect_equal(fit_orig$scale, fit_perm$scale[inv_perm], tolerance = 1e-6)
})

test_that("big_pls_cox approximates plsRcox", {
  skip_if_not_installed("survival")
  skip_if_not_installed("plsRcox")
  skip_if_not_installed("bigmemory")
  
  set.seed(123)
  n <- 30
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n)
  status <- rbinom(n, 1, 0.6)
  
  my <- big_pls_cox(X, time, status, ncomp = 2)
  myother <- plsRcox::plsRcox(X, time, status, nt = 3, verbose = FALSE)
  
  expect_equal(ncol(my$scores), 2)
  expect_equal(nrow(my$loadings), p)
  
  # Compare component spaces via correlations
  cors <- cor(my$scores, myother$tt[, 1:2])
  expect_true(all(abs(cors) <= 1 + 1e-8))
  expect_true(mean(abs(diag(cors))) > 0.95)
})


test_that("prediction helpers return expected shapes", {
  skip_if_not_installed("survival")
  skip_if_not_installed("bigmemory")
  set.seed(2025)
  n <- 40
  p <- 6
  X <- bigmemory::as.big.matrix(matrix(rnorm(n * p), nrow = n))
  time <- rexp(n)
  status <- rbinom(n, 1, 0.6)
  fit <- big_pls_cox(X, time, status, ncomp = 3)
  risks <- predict(fit, type = "risk")
  expect_length(risks, n)
  comps <- predict(fit, type = "components")
  expect_equal(dim(comps), c(n, 3))
  info <- select_ncomp(fit)
  expect_s3_class(info$summary, "data.frame")
  expect_equal(info$summary$ncomp, seq_len(3))
})
