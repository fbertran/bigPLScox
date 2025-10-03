<!-- README.md is generated from README.Rmd. Please edit that file -->




# bigPLScox <img src="man/figures/logo_bigPLS.svg" align="right" width="200"/>

# bigPLScox, PLS models and their extension for big data in R
## Frédéric Bertrand and Myriam Maumy-Bertrand

<!-- badges: start -->
[![R-CMD-check](https://github.com/fbertran/bigPLScox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/bigPLScox/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/fbertran/bigPLScox/actions/workflows/rhub.yaml/badge.svg)](https://github.com/fbertran/bigPLScox/actions/workflows/rhub.yaml)
<!-- badges: end -->


The goal of bigPLScox is provide Cox models in a high dimensional setting in R.



Support for parallel computation and GPU is being developed.


This website and these examples were created by F. Bertrand and M. Maumy-Bertrand.

## Installation

You can install the released version of bigPLScox from [CRAN](https://CRAN.R-project.org) with:


``` r
install.packages("bigPLS")
```

You can install the development version of bigPLScox from [github](https://github.com) with:


``` r
devtools::install_github("fbertran/bigPLS")
```

# Example

## Allelotyping real dataset

### The dataset


``` r
library(plsRcox)
data(micro.censure)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]
Y_test_micro <- micro.censure$survyear[81:117]
C_test_micro <- micro.censure$DC[81:117]

data(Xmicro.censure_compl_imp)
X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)

X_train_micro_orig <- Xmicro.censure_compl_imp[1:80,]
X_train_micro_orig_df <- data.frame(X_train_micro_orig)
```

### Compute deviance residuals with some options.


``` r
head(computeDR(Y_train_micro,C_train_micro,plot=TRUE))
```

## Simulated data

### Generate dataset


``` r
set.seed(4669)
library(bigPLScox)
x_sim <- matrix(sample(0:1, size = 20000, replace = TRUE), ncol = 2)
dCox_sim <- dataCox(10^4, lambda = 3, rho = 2, x_sim,
beta = c(1,3), cens.rate = 5)
data(dCox_sim)
```

### Compute deviance residuals with some options.


``` r
with(dCox_sim,head(computeDR(time,status,plot=TRUE)))
```

### Model Matrix


``` r
coxgpls(~.,Y_train_micro,C_train_micro,ncomp=6,trace=TRUE,model_matrix=TRUE,dataXplan = X_train_micro_orig_df,ind.block.x=c(3,10,20))
```

### coxgpls


``` r
(cox_gpls_fit=coxgpls(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_gpls_fit2=coxgpls(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_gpls_fit3=coxgpls(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,ind.block.x=c(3,10,15)))

rm(cox_gpls_fit,cox_gpls_fit2,cox_gpls_fit3)
```

### cv.coxgpls










































