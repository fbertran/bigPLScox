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
  
  ours <- big_pls_cox(X, time, status, ncomp = 2)
  theirs <- plsRcox::plsRcox(time, status, X, nt = 2)
  
  expect_equal(ncol(ours$scores), 2)
  expect_equal(nrow(ours$loadings), p)
  
  # Compare component spaces via correlations
  cors <- cor(ours$scores, theirs$tt_comp[, 1:2])
  expect_true(all(abs(cors) <= 1 + 1e-8))
  expect_true(mean(abs(diag(cors))) > 0.5)
})
