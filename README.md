<!-- README.md is generated from README.Rmd. Please edit that file -->



# bigPLScox <img src="man/figures/logo_bigPLScox.svg" align="right" width="200"/>

## PLS models for Cox regression with big data in R  
### Frédéric Bertrand and Myriam Maumy-Bertrand

<https://doi.org/10.32614/CRAN.package.bigPLScox>

<!-- badges: start -->
[![DOI](https://img.shields.io/badge/doi-10.32614/CRAN.package.bigPLScox-blue.svg)](https://doi.org/10.32614/CRAN.package.bigPLScox)
[![R-CMD-check](https://github.com/fbertran/bigPLScox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/bigPLScox/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/fbertran/bigPLScox/actions/workflows/rhub.yaml/badge.svg)](https://github.com/fbertran/bigPLScox/actions/workflows/rhub.yaml)
<!-- badges: end -->

`bigPLScox` provides Partial Least Squares (PLS) methods for Cox proportional
hazards models, with a particular focus on high dimensional and big memory
settings. The package supports classical PLS Cox methods together with
accelerated C++ backends that operate directly on `bigmemory::big.matrix`
objects.

The main design goals are:

* Efficient PLS based Cox models for large p and large n  
* First class support for file backed big matrices  
* Unified prediction, cross validation and diagnostic tools

Standalone benchmarking scripts that complement the vignette live under
`inst/benchmarks/`.

The documentation website and examples are maintained by Frédéric
Bertrand and Myriam Maumy.

> **Conference highlight.** Maumy, M. and Bertrand, F. (2023). 
*"PLS models and their extension for big data"*. Conference presentation at the 
Joint Statistical Meetings (JSM 2023), Toronto, Ontario, Canada, Aug 5–10, 2023.

> **Conference highlight.** Maumy, M. and Bertrand, F. (2023). 
*"bigPLS: Fitting and cross-validating PLS-based Cox models to censored big data"*. 
Poster at BioC2023: The Bioconductor Annual Conference, Dana-Farber Cancer Institute, 
Boston, MA, USA, Aug 2–4, 2023. doi:10.7490/f1000research.1119546.1.


## Core modelling functions

The following families of PLS Cox estimators are available.

* `coxgpls()` and `coxgplsDR()`  
  Generalised PLS Cox regression based on partial likelihood, with an optional
  deviance residual based variant (`coxgplsDR`).

* `coxsgpls()` and `coxsgplsDR()`  
  Sparse PLS Cox estimators that encourage variable selection at the latent
  component level.

* `coxspls_sgpls()` and `coxspls_sgplsDR()`  
  Structured sparse PLS Cox versions that support group information.

* DK style estimators  
  `coxDKgplsDR()`, `coxDKsgplsDR()` and `coxDKspls_sgplsDR()` implement
  deviance residual based variants following the DK strategy.

All these functions come in both `default` and `formula` interfaces and have
matching `predict()` methods with support for `type = "link"`, `"risk"` and
other standard Cox outputs.

Cross validation helpers are provided through:

* `cv.coxgpls()`, `cv.coxgplsDR()`  
* `cv.coxsgpls()`, `cv.coxsgplsDR()`  
* `cv.coxspls_sgpls()` and `cv.coxspls_sgplsDR()`  
* `cv.coxDKgplsDR()`, `cv.coxDKsgplsDR()`, `cv.coxDKspls_sgplsDR()`

These mirror the criteria used in `plsRcox` and include time dependent
survival metrics.

## Big memory PLS Cox backends

The package offers dedicated functions for Cox PLS fits on large matrices,
including file backed `bigmemory::big.matrix` objects.

* `big_pls_cox()`  
  Iterative construction of PLS components for Cox models using big matrices,
  with optional naive sparsity through `keepX`.

* `big_pls_cox_fast()`  
  High performance exact PLS Cox backend. It operates on both standard dense
  matrices and `big.matrix` inputs and is implemented entirely in C++ for
  speed.

* `big_pls_cox_gd()`  
  Gradient based optimisation of the Cox partial likelihood in the latent PLS
  space. The `method` argument selects the optimisation scheme:

  - `"gd"` for a basic fixed step gradient descent  
  - `"bb"` for a Barzilai Borwein step size  
  - `"nesterov"` for Nesterov style acceleration  
  - `"bfgs"` for a quasi Newton type update

  All optimisation methods share the same PLS scores and differ only in how
  the Cox coefficients are updated.

* `big_pls_cox_transform()`  
  Low level interface that applies a trained PLS Cox transformation to new
  data, used internally by the prediction helpers and also exported for
  advanced workflows.

Cross validation for the big memory backends is provided by:

* `cv.big_pls_cox()`  
* `cv.big_pls_cox_gd()`

These functions help select the number of components and compare the exact and
gradient based backends.

## Prediction, plots and summaries

The following S3 methods are provided for PLS Cox fits.

* `predict.big_pls_cox()`  
  Prediction method for the original big memory PLS Cox solver.

* `predict.big_pls_cox_fast()`  
  Unified prediction interface for exact PLS Cox fits on both dense and big
  matrices. Supports:
  - `type = "link"`, `"risk"`, `"response"`  
  - `type = "components"` to return PLS scores  
  - `comps` to select a subset of components  
  - `coef` to supply custom Cox coefficients

* `predict.big_pls_cox_gd()`  
  Prediction for gradient based fits that supports the same `type`, `comps`
  and `coef` arguments and uses the stored Cox fit by default.

* `plot.big_pls_cox()` and `plot.big_pls_cox_gd()`  
  Simple visual summaries of component effects, often used together with
  deviance residual plots.

* `summary.big_pls_cox()`, `summary.big_pls_cox_fast()` and
  `summary.big_pls_cox_gd()`  
  Text summaries that expose the PLS structure, number of components, and the
  embedded Cox fit.

* `print.big_pls_cox()`, `print.big_pls_cox_gd()` and
  `print.summary.big_pls_cox_fast()`  
  Compact console output for quick inspection.

Several internal PLS models from `plsRcox` (for example `gPLS`, `sPLS`,
`sgPLS`, `pls.cox`) also have `stats::predict()` methods registered in the
namespace so that standard `predict()` calls continue to work.

## Diagnostics and model selection

`bigPLScox` provides a range of tools for residual diagnostics, component
selection and inspection of gradient based fits.

* Deviance residual tools  
  - `computeDR()` carries out deviance residual computation and can use a
    pure R or C++ engine, with optional support for big matrices.  
  - `cox_deviance_residuals()` and `cox_deviance_residuals_big()` implement
    low level deviance residuals for dense and big memory data.  
  - `cox_partial_deviance_big()` and `cox_deviance_details()` expose partial
    deviance and internal calculations.  
  - `benchmark_deviance_residuals()` provides a simple wrapper to compare
    different implementations on synthetic data.

* Component summaries  
  - `component_information()` extracts per component information such as
    variance explained and effective variable usage from both `big_pls_cox`
    and `big_pls_cox_gd` fits.  
  - `select_ncomp()` offers information criteria based choices for the number
    of components, for example AIC or BIC like rules.

* Gradient based diagnostics  
  - `gd_diagnostics()` returns optimisation diagnostics for gradient based
    backends, including iteration counts, log likelihood progression, gradient
    norms and step sizes.

These tools are intended to complement classic survival model diagnostics such
as `survival::coxph()` residual plots.

## Utilities, data and scaling

A small number of helper functions and data objects round out the package.

* `bigscale`  
  Scaling of big matrices that is compatible with the big memory PLS Cox
  workflow.

* `bigSurvSGD.na.omit()` and `partialbigSurvSGDv0()`  
  Interfaces for survival stochastic gradient methods provided by the
  companion `bigSurvSGD` package.

* `dataCox`  
  Example survival dataset used in documentation and unit tests.

The package also re exports the `%*%` and `Arith` methods used with some
big matrix types.

## Vignettes and learning material

Several vignettes ship with the package and are accessible once it is
installed.

* Getting started with `bigPLScox`  
* Overview of the main modelling functions and their extensions  
* Big memory workflows with `bigmemory` matrices  
* Benchmarking `bigPLScox` against baseline Cox implementations

Refer to the pkgdown site for rendered versions of these documents and a
complete function reference:

<https://fbertran.github.io/bigPLScox/>

## Installation

You can install the released version of bigPLScox from
[CRAN](https://CRAN.R-project.org) with:


``` r
install.packages("bigPLScox")
```

You can install the development version of bigPLScox from
[GitHub](https://github.com/fbertran/bigPLScox) with:


``` r
# install.packages("devtools")
devtools::install_github("fbertran/bigPLScox")
```

## Minimal example

The following minimal example uses the micro array data bundled with the
package.


``` r
library(bigPLScox)
data(micro.censure)
data(Xmicro.censure_compl_imp)

Y <- micro.censure$survyear
status <- micro.censure$DC
X <- Xmicro.censure_compl_imp

set.seed(123)
fit <- coxgpls(
  Xplan = X,
  time = Y,
  status = status,
  ncomp = 4,
  ind.block.x = c(3, 10, 20)
)
#> Error in colMeans(x, na.rm = TRUE): 'x' must be numeric

summary(fit)
#> Error: object 'fit' not found
```

A big memory workflow uses `bigmemory::big.matrix` objects.


``` r
library(bigmemory)

X_big <- bigmemory::as.big.matrix(X)

fast_fit <- big_pls_cox_fast(
  X = X_big,
  time = Y,
  status = status,
  ncomp = 4
)

lp <- predict(fast_fit, newdata = X_big, type = "link")
head(lp)
#> [1] -0.4296294 -0.7809034  1.6411946 -1.3885315  1.2299486 -1.7144312
```

For more elaborate examples, including cross validation and comparisons between
the exact and gradient based backends, see the vignettes and the scripts under
`inst/benchmarks`.

## Citation

If you use `bigPLScox` in scientific work, please cite the package and the
associated conference material.

Maumy, M. and Bertrand, F. (2023). PLS models and their extension for big
data. Joint Statistical Meetings, Toronto, Ontario, Canada.

Maumy, M. and Bertrand, F. (2023). bigPLS: Fitting and cross validating PLS
based Cox models to censored big data. BioC2023, Dana Farber Cancer Institute,
Boston, MA, poster contribution. doi:10.7490/f1000research.1119546.1.
