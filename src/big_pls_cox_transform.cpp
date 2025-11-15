#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::depends(Rcpp, bigmemory)]]

namespace {

// core implementation for both dense and big backends
template <typename ColumnGetter>
NumericMatrix transform_impl(
    ColumnGetter get_x,        // function: (std::size_t i, std::size_t j) -> double
    std::size_t n,             // number of rows
    std::size_t p,             // number of columns
    const NumericVector& means,
    const NumericVector& sds,
    const NumericMatrix& weights,
    const NumericMatrix& loadings,
    const IntegerVector& comps)
{
  const int total_components = weights.ncol();
  if (loadings.ncol() != total_components) {
    stop("weights and loadings must have the same number of columns");
  }
  if (means.size() != static_cast<int>(p) || sds.size() != static_cast<int>(p)) {
    stop("Length of 'means' and 'sds' must match the number of predictors");
  }
  
  // sort requested components but remember original order
  IntegerVector comps_sorted = clone(comps);
  std::sort(comps_sorted.begin(), comps_sorted.end());
  for (int idx = 0; idx < comps_sorted.size(); ++idx) {
    if (comps_sorted[idx] < 1 || comps_sorted[idx] > total_components) {
      stop("Component indices are out of bounds");
    }
  }
  
  const int K = comps.size();
  NumericMatrix scores(n, K);
  
  // work buffer for the *ordered* scores (by comps_sorted)
  NumericMatrix scores_ord(n, K);
  
  // loop over components in sorted order (for correct deflation)
  for (int ord_idx = 0; ord_idx < comps_sorted.size(); ++ord_idx) {
    const int comp_id = comps_sorted[ord_idx] - 1; // 0-based
    
    // temporary score vector for this component
    std::vector<double> t_raw(n, 0.0);
    
    // 1) build raw scores with deflation in predictor space
    for (std::size_t j = 0; j < p; ++j) {
      const double w_j = weights(j, comp_id);
      if (w_j == 0.0)
        continue;
      
      const double mu_j = means[j];
      const double sd_j = sds[j];
      
      for (std::size_t i = 0; i < n; ++i) {
        double value = (get_x(i, j) - mu_j) / sd_j;
        
        // subtract contributions of previous components (already in scores_ord)
        for (int prev_idx = 0; prev_idx < ord_idx; ++prev_idx) {
          const int prev_comp = comps_sorted[prev_idx] - 1;
          value -= scores_ord(i, prev_idx) * loadings(j, prev_comp);
        }
        
        t_raw[i] += value * w_j;
      }
    }
    
    // 2) center + variance-normalise (sample variance == 1)
    double mean_t = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      mean_t += t_raw[i];
    }
    mean_t /= static_cast<double>(n);
    
    for (std::size_t i = 0; i < n; ++i) {
      t_raw[i] -= mean_t;
    }
    
    double ss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      ss += t_raw[i] * t_raw[i];
    }
    
    double denom = 1.0;
    if (n > 1) {
      denom = std::sqrt(ss / static_cast<double>(n - 1));
      if (!std::isfinite(denom) || denom <= 1e-12) {
        denom = 1.0;
      }
    }
    
    for (std::size_t i = 0; i < n; ++i) {
      const double t_std = t_raw[i] / denom;
      scores_ord(i, ord_idx) = t_std;
    }
  }
  
  // 3) if comps were already sorted, return scores_ord directly
  bool already_sorted = true;
  for (int k = 0; k < K; ++k) {
    if (comps[k] != comps_sorted[k]) {
      already_sorted = false;
      break;
    }
  }
  
  if (already_sorted) {
    return scores_ord;
  }
  
  // 4) otherwise, reorder the columns back to the user-requested order
  for (int k = 0; k < K; ++k) {
    const int target_comp = comps[k];
    int src_idx = -1;
    for (int ord_idx = 0; ord_idx < K; ++ord_idx) {
      if (comps_sorted[ord_idx] == target_comp) {
        src_idx = ord_idx;
        break;
      }
    }
    if (src_idx < 0) {
      stop("Internal error while reordering components");
    }
    for (std::size_t i = 0; i < n; ++i) {
      scores(i, k) = scores_ord(i, src_idx);
    }
  }
  
  return scores;
}

} // end anonymous namespace


//------------------------------------------------------------------------------
// Dense (base R matrix) backend
//------------------------------------------------------------------------------

// [[Rcpp::export(name = "big_pls_cox_transform_dense_cpp")]]
NumericMatrix big_pls_cox_transform_dense_cpp(
    NumericMatrix X,
    NumericVector means,
    NumericVector sds,
    NumericMatrix weights,
    NumericMatrix loadings,
    IntegerVector comps)
{
  const std::size_t n = X.nrow();
  const std::size_t p = X.ncol();
  
  // getter: value (i,j) from dense matrix
  auto get_x = [&X](std::size_t i, std::size_t j) -> double {
    return X(i, j);
  };
  
  return transform_impl(get_x, n, p, means, sds, weights, loadings, comps);
}


//------------------------------------------------------------------------------
// big.matrix backend
//------------------------------------------------------------------------------

// [[Rcpp::export(name = "big_pls_cox_transform_big_cpp")]]
NumericMatrix big_pls_cox_transform_big_cpp(
    SEXP xpMat,
    NumericVector means,
    NumericVector sds,
    NumericMatrix weights,
    NumericMatrix loadings,
    IntegerVector comps)
{
  XPtr<BigMatrix> pMat(xpMat);
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  
  MatrixAccessor<double> acc(*pMat);
  
  // getter: value (i,j) from big.matrix
  auto get_x = [&acc](std::size_t i, std::size_t j) -> double {
    return acc[j][i];
  };
  
  return transform_impl(get_x, n, p, means, sds, weights, loadings, comps);
}