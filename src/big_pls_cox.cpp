#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp, bigmemory)]]

// Compute column-wise means and standard deviations without loading the whole matrix
// [[Rcpp::export(name = "big_pls_cox_col_stats_cpp")]]
List big_pls_cox_col_stats_cpp(SEXP xpMat) {
  XPtr<BigMatrix> pMat(xpMat);
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  MatrixAccessor<double> accessor(*pMat);
  
  NumericVector means(p);
  NumericVector sds(p);
  
  for (std::size_t j = 0; j < p; ++j) {
    double sum = 0.0;
    double sumsq = 0.0;
    double* col = accessor[j];
    for (std::size_t i = 0; i < n; ++i) {
      double value = col[i];
      sum += value;
      sumsq += value * value;
    }
    double mean = sum / static_cast<double>(n);
    double var = 0.0;
    if (n > 1) {
      var = (sumsq - static_cast<double>(n) * mean * mean) /
        static_cast<double>(n - 1);
    }
    if (!std::isfinite(var) || var <= 0.0) {
      var = 0.0;
    }
    means[j] = mean;
    sds[j] = (var > 0.0) ? std::sqrt(var) : 1.0;
  }
  
  return List::create(Named("mean") = means,
                      Named("sd") = sds);
}

// Compute the next PLS component given martingale residuals
// [[Rcpp::export(name = "big_pls_cox_component_cpp")]]
List big_pls_cox_component_cpp(SEXP xpMat,
                               NumericVector residuals,
                               NumericMatrix scores_prev,
                               NumericMatrix loadings_prev,
                               NumericVector means,
                               NumericVector sds) {
  XPtr<BigMatrix> pMat(xpMat);
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  const int n_prev = scores_prev.ncol();
  
  if (static_cast<std::size_t>(residuals.size()) != n) {
    stop("Residual vector length does not match number of rows");
  }
  if (loadings_prev.nrow() != static_cast<int>(p)) {
    stop("Loadings matrix must have one row per predictor");
  }
  
  MatrixAccessor<double> accessor(*pMat);
  
  NumericVector weights(p);
  NumericVector score(n);
  NumericVector loading(p);
  
  // Compute weights
  double norm_sq = 0.0;
  for (std::size_t j = 0; j < p; ++j) {
    double* col = accessor[j];
    const double mean = means[j];
    const double sd = sds[j];
    double accum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double value = (col[i] - mean) / sd;
      for (int h = 0; h < n_prev; ++h) {
        value -= scores_prev(i, h) * loadings_prev(j, h);
      }
      accum += residuals[i] * value;
    }
    weights[j] = accum;
    norm_sq += accum * accum;
  }
  
  if (norm_sq <= 0.0) {
    stop("Unable to compute weight vector; residuals may be zero");
  }
  
  const double norm = std::sqrt(norm_sq);
  for (std::size_t j = 0; j < p; ++j) {
    weights[j] /= norm;
  }
  
  // Compute the new score vector
  for (std::size_t i = 0; i < n; ++i) {
    double accum = 0.0;
    for (std::size_t j = 0; j < p; ++j) {
      double* col = accessor[j];
      const double mean = means[j];
      const double sd = sds[j];
      double value = (col[i] - mean) / sd;
      for (int h = 0; h < n_prev; ++h) {
        value -= scores_prev(i, h) * loadings_prev(j, h);
      }
      accum += value * weights[j];
    }
    score[i] = accum;
  }
  
  double score_norm_sq = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    score_norm_sq += score[i] * score[i];
  }
  if (score_norm_sq <= 0.0) {
    stop("Computed score vector has zero variance");
  }
  
  // Compute the loading vector
  for (std::size_t j = 0; j < p; ++j) {
    double* col = accessor[j];
    const double mean = means[j];
    const double sd = sds[j];
    double accum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double value = (col[i] - mean) / sd;
      for (int h = 0; h < n_prev; ++h) {
        value -= scores_prev(i, h) * loadings_prev(j, h);
      }
      accum += value * score[i];
    }
    loading[j] = accum / score_norm_sq;
  }
  
  return List::create(Named("weights") = weights,
                      Named("scores") = score,
                      Named("loadings") = loading);
}
