test_that("Test big features", {
skip_if_not_installed("bigmemory")
skip_if_not_installed("survival")

set.seed(42)
n <- 40
p <- 6
X <- matrix(rnorm(n * p), n, p)
time <- rexp(n)
status <- rbinom(n, 1, 0.7)

X_big <- bigmemory::as.big.matrix(X)

# big_pls_cox supports prediction and sparse selection
fit_stream <- big_pls_cox(X_big, time, status, ncomp = 2, keepX = c(3, 2))
fit_gd <- big_pls_cox_gd(X_big, time, status, ncomp = 2, max_iter = 50, keepX = 2)

scores_stream <- predict(fit_stream, type = "components")
expect_equal(nrow(scores_stream), n)
expect_equal(ncol(scores_stream), 2)

link_stream <- predict(fit_stream, type = "link")
link_gd <- predict(fit_gd, type = "link")
expect_length(link_stream, n)
expect_length(link_gd, n)

# Deviance residuals from C++ engine match survival::coxph
res_surv <- suppressWarnings(computeDR(time, status, engine = "survival"))
res_cpp <- computeDR(time, status, engine = "cpp", eta = rep(1,length(status)))
expect_equal(as.numeric(res_cpp), as.numeric(res_surv), tolerance = 1e-7)

# Information criteria helpers
info <- component_information(fit_stream, max_comp = 2)
expect_s3_class(info, "data.frame")
expect_equal(info$ncomp, 1:2)
sel <- select_ncomp(fit_stream, criterion = "AIC", max_comp = 2)
expect_true(is.list(sel))
expect_true(all(c("summary", "criterion", "opt_ncomp") %in% names(sel)))

# Cross-validation handles big.matrix inputs
cv_fit <- suppressWarnings(cv.coxgpls(list(x = X_big, time = time, status = status), 
                                      nt = 2,plot.it = FALSE))
expect_true(is.list(cv_fit))
})
