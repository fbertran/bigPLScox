#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace Rcpp;

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

namespace {

std::vector<bool> build_keep_mask(const NumericVector& weights, int keep) {
  const std::size_t p = weights.size();
  std::vector<bool> mask(p, true);
  if (keep <= 0 || keep >= static_cast<int>(p)) {
    return mask;
  }
  std::vector<std::size_t> indices(p);
  std::iota(indices.begin(), indices.end(), 0);
  auto comp = [&](std::size_t lhs, std::size_t rhs) {
    return std::fabs(weights[lhs]) > std::fabs(weights[rhs]);
  };
  std::partial_sort(indices.begin(), indices.begin() + keep, indices.end(), comp);
  std::fill(mask.begin(), mask.end(), false);
  for (int k = 0; k < keep; ++k) {
    mask[indices[k]] = true;
  }
  return mask;
}
  
  template <typename Accessor>
  NumericMatrix compute_scores_impl(Accessor& accessor,
                                    std::size_t n,
                                    std::size_t p,
                                    const NumericVector& means,
                                    const NumericVector& sds,
                                    const NumericMatrix& weights,
                                    const NumericMatrix& loadings,
                                    const IntegerVector& comps) {
    const int total_components = weights.ncol();
    if (loadings.ncol() != total_components) {
      stop("weights and loadings must have the same number of columns");
    }
    if (means.size() != p || sds.size() != p) {
      stop("Length of 'means' and 'sds' must match the number of predictors");
    }
    
    IntegerVector comps_sorted = clone(comps);
    std::sort(comps_sorted.begin(), comps_sorted.end());
    for (int idx = 0; idx < comps_sorted.size(); ++idx) {
      if (comps_sorted[idx] < 1 || comps_sorted[idx] > total_components) {
        stop("Component indices are out of bounds");
      }
    }
    
    const int k = comps.size();
    NumericMatrix scores(n, k);
    
    for (int order_idx = 0; order_idx < comps_sorted.size(); ++order_idx) {
      const int comp_id = comps_sorted[order_idx] - 1;
      NumericVector score(n);
      
      for (std::size_t j = 0; j < p; ++j) {
        const double w = weights(j, comp_id);
        if (w == 0.0) {
          continue;
        }
        const double mean = means[j];
        const double sd = sds[j];
        
        auto column = accessor.column(j);
        for (std::size_t i = 0; i < n; ++i) {
          double value = (column[i] - mean) / sd;
          for (int prev_idx = 0; prev_idx < order_idx; ++prev_idx) {
            const int prev_comp = comps_sorted[prev_idx] - 1;
            value -= scores(i, prev_idx) * loadings(j, prev_comp);
          }
          score[i] += value * w;
        }
      }
      
      for (std::size_t i = 0; i < n; ++i) {
        scores(i, order_idx) = score[i];
      }
    }
    
    if (!std::is_sorted(comps.begin(), comps.end())) {
      NumericMatrix reordered(n, k);
      for (int idx = 0; idx < k; ++idx) {
        int target = -1;
        for (int order_idx = 0; order_idx < comps_sorted.size(); ++order_idx) {
          if (comps_sorted[order_idx] == comps[idx]) {
            target = order_idx;
            break;
          }
        }
        if (target < 0) {
          stop("Internal error while reordering components");
        }
        for (std::size_t i = 0; i < n; ++i) {
          reordered(i, idx) = scores(i, target);
        }
      }
      return reordered;
    }
    
    return scores;
  }
  
  class BigMatrixAccessorWrapper {
  public:
    explicit BigMatrixAccessorWrapper(BigMatrix& mat) : accessor_(mat) {}
    
    double* column(std::size_t j) {
      return accessor_[j];
    }
    
  private:
    MatrixAccessor<double> accessor_;
  };
  
  class NumericMatrixAccessorWrapper {
  public:
    explicit NumericMatrixAccessorWrapper(NumericMatrix mat) : mat_(mat) {}
    
    double* column(std::size_t j) {
      return mat_.begin() + static_cast<R_xlen_t>(j) * mat_.nrow();
    }
    
  private:
    NumericMatrix mat_;
  };
  
} // namespace

// Compute the next PLS component given martingale residuals
// [[Rcpp::export(name = "big_pls_cox_component_cpp")]]
List big_pls_cox_component_cpp(SEXP xpMat,
                               NumericVector residuals,
                               NumericMatrix scores_prev,
                               NumericMatrix loadings_prev,
                               NumericVector means,
                               NumericVector sds,
                               int keepX) {
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
  }
  
  std::vector<bool> mask = build_keep_mask(weights, keepX);
  double norm_sq = 0.0;
  for (std::size_t j = 0; j < p; ++j) {
    if (!mask[j]) {
      weights[j] = 0.0;
      continue;
    }
    norm_sq += weights[j] * weights[j];
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

// [[Rcpp::export(name = "big_pls_cox_transform_cpp")]]
NumericMatrix big_pls_cox_transform_cpp(SEXP xpMat,
                                        NumericVector means,
                                        NumericVector sds,
                                        NumericMatrix weights,
                                        NumericMatrix loadings,
                                        IntegerVector comps) {
  XPtr<BigMatrix> pMat(xpMat);
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  BigMatrixAccessorWrapper accessor(*pMat);
  return compute_scores_impl(accessor, n, p, means, sds, weights, loadings, comps);
}

// [[Rcpp::export(name = "matrix_pls_cox_transform_cpp")]]
NumericMatrix matrix_pls_cox_transform_cpp(NumericMatrix X,
                                           NumericVector means,
                                           NumericVector sds,
                                           NumericMatrix weights,
                                           NumericMatrix loadings,
                                           IntegerVector comps) {
  const std::size_t n = X.nrow();
  const std::size_t p = X.ncol();
  NumericMatrixAccessorWrapper accessor(X);
  return compute_scores_impl(accessor, n, p, means, sds, weights, loadings, comps);
}
