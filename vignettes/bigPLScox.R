## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dpi = 150,
  message = FALSE,
  warning = FALSE
)

## ----eval = FALSE-------------------------------------------------------------
# # From CRAN (once released)
# install.packages("bigPLScox")
# 
# # Development version
# # remotes::install_github("fbertran/bigPLScox")

## -----------------------------------------------------------------------------
library(bigPLScox)

data(micro.censure)
data(Xmicro.censure_compl_imp)

Y_train <- micro.censure$survyear[1:80]
status_train <- micro.censure$DC[1:80]
X_train <- Xmicro.censure_compl_imp[1:80, -40]

## -----------------------------------------------------------------------------
residuals_plot <- computeDR(Y_train, status_train, plot = TRUE)
head(residuals_plot)

## -----------------------------------------------------------------------------
set.seed(123)
fit <- coxgpls(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 20)
)

fit

## -----------------------------------------------------------------------------
set.seed(123)
cv_results <- cv.coxgpls(
  data = list(x = X_train, time = Y_train, status = status_train),
  nt = 6,
  ind.block.x = c(3, 10, 20),
)

cv_results$opt_nt

## -----------------------------------------------------------------------------
data(micro.censure)
data(Xmicro.censure_compl_imp)

X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]


fit_dr=coxgplsDR(X_train_micro,Y_train_micro,C_train_micro,
ncomp=6,ind.block.x=c(3,10,15),keepX=rep(4,6))

fit_dr

