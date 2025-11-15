#ifndef FAST_COMMON_H
#define FAST_COMMON_H

#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// =========================================================
// Cox score residuals for a given eta
// r_i = status_i – exp(eta_i)/sum_{j: t_j >= t_i} exp(eta_j)
// =========================================================
inline arma::vec cox_score_residuals(
    const arma::vec& time,
    const arma::vec& status,
    const arma::vec& eta)
{
  const std::size_t n = time.n_elem;
  arma::uvec ord = arma::sort_index(time, "descend");
  
  arma::vec exp_eta = arma::exp(eta);
  arma::vec r(n, arma::fill::zeros);
  
  double running = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = ord[k];
    running += exp_eta[i];
    r[i] = status[i] - exp_eta[i] / running;
  }
  return r;
}

// =========================================================
// Build keepX mask
// =========================================================
inline std::vector<bool> build_mask(const arma::vec& w, int keepX) {
  const std::size_t p = w.n_elem;
  std::vector<bool> mask(p, true);
  if (keepX <= 0 || keepX >= static_cast<int>(p)) return mask;
  
  arma::uvec ord = arma::sort_index(arma::abs(w), "descend");
  for (std::size_t k = keepX; k < p; ++k) {
    mask[ord[k]] = false;
  }
  return mask;
}

// =========================================================
// Compute one PLS component: weights, scores, loadings
// ColFun x(j,i) must return standardised X[i,j] as double
// =========================================================
template <class ColFun>
inline void compute_component(
    std::size_t h,
    std::size_t n,
    std::size_t p,
    const arma::vec& r,
    const arma::mat& scores_prev,
    const arma::mat& load_prev,
    int keepX,
    ColFun x,
    arma::vec& w,
    arma::vec& t,
    arma::vec& load)
{
  // ---- weights ----
  for (std::size_t j = 0; j < p; ++j) {
    double acc = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double v = x(j, i);
      for (std::size_t k = 0; k < h; ++k) {
        v -= scores_prev(i, k) * load_prev(j, k);
      }
      acc += r[i] * v;
    }
    w[j] = acc;
  }
  
  // keepX masking
  std::vector<bool> mask = build_mask(w, keepX);
  double ss = 0.0;
  for (std::size_t j = 0; j < p; ++j) {
    if (!mask[j]) {
      w[j] = 0.0;
      continue;
    }
    ss += w[j] * w[j];
  }
  
  double wn = std::sqrt(ss);
  if (!std::isfinite(wn) || wn < 1e-15) wn = 1.0;
  w /= wn;
  
  // deterministic sign: force sum(w) >= 0
  if (arma::sum(w) < 0.0) w = -w;
  
  // ---- scores t = X_deflated * w ----
  for (std::size_t i = 0; i < n; ++i) {
    double acc = 0.0;
    for (std::size_t j = 0; j < p; ++j) {
      double v = x(j, i);
      for (std::size_t k = 0; k < h; ++k) {
        v -= scores_prev(i, k) * load_prev(j, k);
      }
      acc += v * w[j];
    }
    t[i] = acc;
  }
  
  // center + variance-normalize (sample var = 1)
  t -= arma::mean(t);
  double var_t = arma::var(t);
  if (!std::isfinite(var_t) || var_t < 1e-15) var_t = 1.0;
  t /= std::sqrt(var_t);
  
  // ---- loadings: load = X_deflated^T t / (t^T t) ----
  const double tt = arma::dot(t, t); // ~ n - 1
  for (std::size_t j = 0; j < p; ++j) {
    double acc = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      double v = x(j, i);
      for (std::size_t k = 0; k < h; ++k) {
        v -= scores_prev(i, k) * load_prev(j, k);
      }
      acc += v * t[i];
    }
    load[j] = acc / tt;
  }
}

// =========================================================
// Dense driver: X is an arma::mat
// =========================================================
inline Rcpp::List run_pls_cox_dense(
    const arma::mat& X,
    const arma::vec& time,
    const arma::vec& status,
    int ncomp,
    const arma::vec& means,
    const arma::vec& sds,
    const arma::ivec& keepX)
{
  const std::size_t n = X.n_rows;
  const std::size_t p = X.n_cols;
  
  arma::mat scores(n, ncomp, arma::fill::zeros);
  arma::mat loadings(p, ncomp, arma::fill::zeros);
  arma::mat weights(p, ncomp, arma::fill::zeros);
  
  arma::mat scores_prev(n, 0);
  arma::mat load_prev(p, 0);
  
  // single set of Cox score residuals (eta = 0)
  arma::vec eta0(n, arma::fill::zeros);
  arma::vec r = cox_score_residuals(time, status, eta0);
  
  auto colfun = [&](std::size_t j, std::size_t i) -> double {
    return (X(i, j) - means[j]) / sds[j];
  };
  
  for (int h = 0; h < ncomp; ++h) {
    arma::vec w(p, arma::fill::zeros);
    arma::vec t(n, arma::fill::zeros);
    arma::vec load(p, arma::fill::zeros);
    
    compute_component(
      static_cast<std::size_t>(h),
      n, p, r,
      scores_prev, load_prev,
      keepX[h],
           colfun,
           w, t, load
    );
    
    weights.col(h)  = w;
    scores.col(h)   = t;
    loadings.col(h) = load;
    
    scores_prev = scores.cols(0, h);
    load_prev   = loadings.cols(0, h);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("scores")   = scores,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("weights")  = weights,
    Rcpp::Named("center")   = means,
    Rcpp::Named("scale")    = sds
  );
}

// =========================================================
// Big driver: X is a big.matrix
// =========================================================
inline Rcpp::List run_pls_cox_big(
    SEXP xp,
    const arma::vec& time,
    const arma::vec& status,
    int ncomp,
    const arma::vec& means,
    const arma::vec& sds,
    const arma::ivec& keepX)
{
  Rcpp::XPtr<BigMatrix> pMat(xp);
  MatrixAccessor<double> acc(*pMat);
  
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  
  arma::mat scores(n, ncomp, arma::fill::zeros);
  arma::mat loadings(p, ncomp, arma::fill::zeros);
  arma::mat weights(p, ncomp, arma::fill::zeros);
  
  arma::mat scores_prev(n, 0);
  arma::mat load_prev(p, 0);
  
  // single set of Cox score residuals (eta = 0)
  arma::vec eta0(n, arma::fill::zeros);
  arma::vec r = cox_score_residuals(time, status, eta0);
  
  auto colfun = [&](std::size_t j, std::size_t i) -> double {
    return (acc[j][i] - means[j]) / sds[j];
  };
  
  for (int h = 0; h < ncomp; ++h) {
    arma::vec w(p, arma::fill::zeros);
    arma::vec t(n, arma::fill::zeros);
    arma::vec load(p, arma::fill::zeros);
    
    compute_component(
      static_cast<std::size_t>(h),
      n, p, r,
      scores_prev, load_prev,
      keepX[h],
           colfun,
           w, t, load
    );
    
    weights.col(h)  = w;
    scores.col(h)   = t;
    loadings.col(h) = load;
    
    scores_prev = scores.cols(0, h);
    load_prev   = loadings.cols(0, h);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("scores")   = scores,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("weights")  = weights,
    Rcpp::Named("center")   = means,
    Rcpp::Named("scale")    = sds
  );
}

#endif // FAST_COMMON_H