<!-- README.md is generated from README.Rmd. Please edit that file -->



# bigPLScox <img src="man/figures/logo_bigPLScox.svg" align="right" width="200"/>

# bigPLScox, PLS models and their extension for big data in R
## Frédéric Bertrand and Myriam Maumy-Bertrand

<!-- badges: start -->
[![R-CMD-check](https://github.com/fbertran/bigPLScox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/bigPLScox/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/fbertran/bigPLScox/actions/workflows/rhub.yaml/badge.svg)](https://github.com/fbertran/bigPLScox/actions/workflows/rhub.yaml)
<!-- badges: end -->


`bigPLScox` provides Partial Least Squares (PLS) methods tailored for Cox
proportional hazards models with large, high-dimensional feature matrices. The
package works directly with [`bigmemory`](https://cran.r-project.org/package=bigmemory)
objects, enabling native C++ accelerators and iterative algorithms to run without loading the full dataset
into memory. In addition to the classical `coxgpls()` solver, the package
contains accelerated variants, cross-validation helpers, and model diagnostics.

- Generalised PLS Cox regression via `coxgpls()` with support for grouped predictors.
- Sparse and structured sparse extensions (`coxsgpls()`, `coxspls_sgpls()`).
- Deviance-residual estimators (`coxgplsDR()`, `coxsgplsDR()`) for robust fits.
- Cross-validation helpers (`cv.coxgpls()`, `cv.coxsgpls()`, …) to select the number of latent components.
- Dataset generators, and diagnostics such as `computeDR()` for quick residual exploration.
- Sparse, group-sparse, and stochastic gradient variants able to consume
  file-backed `big.matrix` objects while leveraging `foreach` parallelism.
- Interfaces for big-memory data through `big_pls_cox()` and `big_pls_cox_gd()`.

GPU support is **not** available in the current release; ongoing development
focuses on improving the multi-core CPU back-end instead.

Additional articles are available in the `vignettes/` directory:

* *Getting started with bigPLScox* — a walkthrough of core modelling,
  cross-validation, and diagnostic helpers.
* *Overview of bigPLScox* — a tour of the main modelling functions with
  practical guidance on choosing estimators.
* *Big-memory workflows with bigPLScox* — guidance on using `bigmemory`
  matrices and parallel back-ends.
* *Benchmarking bigPLScox* — reproducible performance comparisons using the
  **bench** package.

Standalone benchmarking scripts that complement the vignette live under
`inst/benchmarks/`.

The documentation website and examples are maintained by Frédéric
Bertrand and Myriam Maumy.

## Key features

* **Scalable Cox-PLS solvers** (`coxgpls()`, `coxgplsDR()`) that operate on big
  matrices stored on disk.
* **Cross-validation tooling** to select the optimal number of PLS components
  with time-dependent performance metrics.
* **Model diagnostics** such as deviance residual visualisation through
  `computeDR()`.
* **Benchmark scripts** in `inst/benchmarks/` to quantify runtime trade-offs
  between the available solvers.
* **Comprehensive vignette** (`vignettes/bigPLScox.Rmd`) showing a complete
  modelling workflow.

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

## Learning materials

* Browse the **Getting started** vignette with `vignette("getting-started",
  package = "bigPLScox")` for a worked example.
* Explore `vignette("bigPLScox", package = "bigPLScox")` for big-memory
  workflows and streaming solvers.
* Consult the function reference at <https://fbertran.github.io/bigPLScox/>.
* Run the benchmarking scripts in `inst/benchmarks/` to compare solver
  performance on simulated data.

# Quick start

The following example demonstrates the typical workflow on a subset of the
allelotyping dataset bundled with the package. Chunks are evaluated by default
when the README is rendered locally, but they can be toggled with
`knitr::opts_chunk$set(eval = FALSE)` for faster builds.


``` r
library(bigPLScox)
data(micro.censure)
data(Xmicro.censure_compl_imp)
Y_train <- micro.censure$survyear[1:80]
status_train <- micro.censure$DC[1:80]
X_train <- Xmicro.censure_compl_imp[1:80, ]
```

Visualise deviance residuals to assess the baseline model fit:

``` r
residuals_overview <- computeDR(Y_train, status_train, plot = TRUE)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-4-1.png" alt="plot of chunk unnamed-chunk-4" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-4</p>
</div>

``` r
head(residuals_overview)
#>          1          2          3          4          5 
#> -1.4843296 -0.5469540 -0.2314550 -0.3400301 -0.9763372 
#>          6 
#> -0.3866766
```

Fit a Cox-PLS model with six components and inspect the fit summary:

``` r
set.seed(123)
cox_pls_fit <- coxgpls(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = 6,
  ind.block.x = c(3, 10, 20)
)
#> Error in colMeans(x, na.rm = TRUE): 'x' must be numeric
cox_pls_fit
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_gpls)
#> 
#>           coef exp(coef) se(coef)      z        p
#> dim.1 -0.53932   0.58314  0.08723 -6.183  6.3e-10
#> dim.2 -0.39532   0.67347  0.10387 -3.806 0.000141
#> dim.3 -0.29623   0.74362  0.10763 -2.752 0.005918
#> dim.4 -0.29523   0.74436  0.11762 -2.510 0.012074
#> dim.5 -0.11801   0.88869  0.09157 -1.289 0.197498
#> dim.6 -0.09332   0.91091  0.10717 -0.871 0.383910
#> 
#> Likelihood ratio test=56.09  on 6 df, p=2.792e-10
#> n= 80, number of events= 80
```

Cross-validate the number of components and re-fit using the deviance residual
solver for comparison:

``` r
set.seed(123)
cv_results <- cv.coxgpls(
  dataY = list(x = X_train, time = Y_train, status = status_train),
  nt = 6,
  ind.block.x = c(3, 10, 20)
)
#> Error in cv.coxgpls(dataY = list(x = X_train, time = Y_train, status = status_train), : argument "data" is missing, with no default
cv_results$opt_nt
#> NULL
cox_pls_dr <- coxgplsDR(
  Xplan = X_train,
  time = Y_train,
  status = status_train,
  ncomp = cv_results$opt_nt,
  ind.block.x = c(3, 10, 20)
)
#> Error in colMeans(x, na.rm = TRUE): 'x' must be numeric
cox_pls_dr
#> Call:
#> coxph(formula = YCsurv ~ ., data = tt_gplsDR)
#> 
#>          coef exp(coef) se(coef)     z        p
#> dim.1 0.53154   1.70155  0.08694 6.114 9.72e-10
#> dim.2 0.39604   1.48593  0.10932 3.623 0.000292
#> dim.3 0.36946   1.44696  0.12793 2.888 0.003876
#> dim.4 0.21361   1.23814  0.09522 2.243 0.024872
#> dim.5 0.12640   1.13474  0.08615 1.467 0.142322
#> dim.6 0.05705   1.05871  0.10539 0.541 0.588302
#> 
#> Likelihood ratio test=54.01  on 6 df, p=7.348e-10
#> n= 80, number of events= 80
```

Explore alternative estimators such as `coxgplsDR()` for deviance-residual fitting or `coxsgpls()` for sparse component selection. Refer to the package reference for the full list of available models and helper functions.

## Benchmarking

We provide reproducible benchmarks that compare `coxgpls()` and the big-memory
solvers against `survival::coxph()`. Start with the **Benchmarking bigPLScox**
vignette for an interactive tour.

For command-line experiments, execute the scripts in `inst/benchmarks/` after
installing the optional dependencies listed under `Suggests` in the
`DESCRIPTION` file. Each script accepts environment variables (for example,
`bigPLScox.benchmark.n`, `bigPLScox.benchmark.p`, and
`bigPLScox.benchmark.ncomp`) to control the simulation size.

```bash
Rscript inst/benchmarks/cox-benchmark.R
Rscript inst/benchmarks/cox_pls_benchmark.R
Rscript inst/benchmarks/benchmark_bigPLScox.R
```

Results are stored under `inst/benchmarks/results/` with time-stamped filenames
for traceability.

## Vignettes and documentation

Four vignettes ship with the package:

1. **Getting started with bigPLScox** – an end-to-end introduction covering data
   preparation, fitting, and validation workflows.
2. **Overview of bigPLScox** – a high-level description of the modelling
   functions and their typical use cases.
3. **Big-memory workflows with bigPLScox** – instructions for working with
   `bigmemory` matrices and the streaming solvers.
4. **Benchmarking bigPLScox** – guidance for evaluating performance against
   baseline Cox implementations using the **bench** package.
   
The full reference documentation and pkgdown website are available at
<https://fbertran.github.io/bigPLScox/>.

## Bug reports and feature requests

Bug reports and feature requests can be
filed on the [issue tracker](https://github.com/fbertran/bigPLScox/issues/). Please make
sure that new code comes with unit tests or reproducible examples when applicable.
