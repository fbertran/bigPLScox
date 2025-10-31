## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/benchmark-",
  out.width = "100%"
)

## ----message=FALSE------------------------------------------------------------
library(bigPLScox)
library(survival)
library(bench)

## -----------------------------------------------------------------------------
set.seed(2024)
sim_design <- dataCox(
  n = 2000,
  lambda = 2,
  rho = 1.5,
  x = matrix(rnorm(2000 * 50), ncol = 50),
  beta = c(1, 3, rep(0, 48)),
  cens.rate = 5
)

## -----------------------------------------------------------------------------
cox_data <- list(
  x = as.matrix(sim_design[,-(1:3)]),
  time = sim_design$time,
  status = sim_design$status
)

## -----------------------------------------------------------------------------
res <- bench::mark(
  bigPLScox = coxgpls(
    cox_data$x,
    cox_data$time,
    cox_data$status,
    ncomp = 5,
    ind.block.x = c(3, 10)
  ),
  survival = coxph(Surv(cox_data$time, cox_data$status) ~ cox_data$x, ties = "breslow"),
  iterations = 100,
  check = FALSE
)
res

## -----------------------------------------------------------------------------
plot(res, type = "jitter")

