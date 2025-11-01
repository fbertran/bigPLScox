# bigPLScox 0.6.0

See the "Release highlights" section of the README for a condensed overview of
these changes.

* Added C++ implementations for Cox deviance residuals with streaming support
  for \pkg{bigmemory} matrices together with benchmarking utilities.
* Introduced prediction wrappers and component selection helpers (AIC/BIC) for
  `big_pls_cox()` and `big_pls_cox_gd()`.
* Enabled naive sparsity control in `big_pls_cox()` and exposed survival model
  objects for downstream predictions.
* Added cross-validation helpers `cv.big_pls_cox()` and `cv.big_pls_cox_gd()`
  mirroring the `plsRcox` criteria, including the recommended survivalROC
  iAUC metric by default.
* Documented the legacy and big-memory prediction helpers with runnable
  examples and cross references to diagnostic utilities.
* Extended unit test coverage for the new deviance and prediction features.
* Fixed `cv.coxgpls()` to accept `big.matrix` predictors without coercion
  errors.

# bigPLScox 0.5.0

* Added reproducible benchmarking utilities under `inst/benchmarks` comparing
  `big_pls_cox()` against `plsRcox::plsRcox()` on in-memory and file-backed
  matrices.
* Published two package vignettes that cover introductory workflows and
  large-scale analyses with `bigmemory`.
* Added an introductory vignette covering the core Cox-PLS workflow.
* Refreshed the README and website copy to highlight core functionality and to
  demonstrate working examples without warnings, including guidance on learning
  materials and benchmarking resources.
* Refreshed the README with guidance on learning materials and benchmarking.
* Completed package-level documentation with bibliographic references.
* Updated package metadata to list optional dependencies used in docs and
  benchmarks.

# bigPLScox 0.4.0

* Updated maintainer contact details in `DESCRIPTION`.
* Added unit tests for `big_pls_cox()` and `big_pls_cox_gd()` stability checks.
* Added unit tests covering the new C++-accelerated Cox PLS implementation and
  cross-validation utilities.

# bigPLScox 0.3.0 

* Improved `big_pls_cox()` numerical stability and added support for additional
  convergence diagnostics in the gradient-descent solver.
* Refactored stochastic gradient solvers to better integrate with
  `bigmemory` file-backed matrices.
* Improved numerical stability of the deviance residual computations.

# bigPLScox 0.2.0 

* Expanded documentation examples for deviance residuals and Cox model
  utilities.
* Added dataset documentation for `micro.censure` and simulated Cox examples.
* Added pkgdown site configuration and continuous integration workflows.

# bigPLScox 0.1.0

* Introduced gPLS and sgPLS model families with support for grouped predictors
  and deviance residual pipelines with cross-validation support.

# bigPLScox 0.0.1 

* Initial package skeleton with core data objects and helper routines.

