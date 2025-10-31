#include <Rcpp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp, bigmemory)]]

namespace {

struct DevianceResult {
  NumericVector deviance;
  NumericVector martingale;
  NumericVector cumhaz;
  NumericVector eta;
};

DevianceResult compute_deviance_core(const NumericVector& time,
                                     const NumericVector& status,
                                     const NumericVector& eta,
                                     const std::string& method) {
  const std::size_t n = static_cast<std::size_t>(time.size());
  if (status.size() != static_cast<R_xlen_t>(n) ||
      eta.size() != static_cast<R_xlen_t>(n)) {
    stop("time, status and eta must have the same length");
  }
  std::vector<double> exp_eta(n);
  for (std::size_t i = 0; i < n; ++i) {
    exp_eta[i] = std::exp(eta[i]);
  }
  
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    if (time[lhs] == time[rhs]) {
      return status[lhs] > status[rhs];
    }
    return time[lhs] < time[rhs];
  });
  
  double risk_sum = std::accumulate(exp_eta.begin(), exp_eta.end(), 0.0);
  NumericVector cum_base(n);
  double baseline_cum = 0.0;
  std::size_t idx = 0;
  
  while (idx < n) {
    const double current_time = time[order[idx]];
    double event_count = 0.0;
    double event_risk = 0.0;
    std::vector<std::size_t> block;
    
    while (idx < n && time[order[idx]] == current_time) {
      const std::size_t id = order[idx];
      block.push_back(id);
      if (status[id] > 0.0) {
        event_count += status[id];
        event_risk += exp_eta[id];
      }
      ++idx;
    }
    
    double hazard_increment = 0.0;
    if (event_count > 0.0) {
      if (method == "breslow") {
        if (risk_sum <= 0.0) {
          stop("Risk set sum became non-positive under Breslow approximation");
        }
        hazard_increment = event_count / risk_sum;
      } else {
        // Efron tie handling (default)
        const int tied_events = static_cast<int>(std::round(event_count));
        for (int l = 0; l < tied_events; ++l) {
          const double denom = risk_sum - (static_cast<double>(l) / event_count) * event_risk;
          if (denom <= 0.0) {
            stop("Risk set sum became non-positive under Efron approximation");
          }
          hazard_increment += 1.0 / denom;
        }
      }
    }
    
    baseline_cum += hazard_increment;
    for (std::size_t id : block) {
      cum_base[id] = baseline_cum;
    }
    
    for (std::size_t id : block) {
      risk_sum -= exp_eta[id];
    }
  }
  
  NumericVector martingale(n);
  NumericVector deviance(n);
  NumericVector cumhaz(n);
  
  for (std::size_t i = 0; i < n; ++i) {
    const double cumhaz_i = cum_base[i] * exp_eta[i];
    cumhaz[i] = cumhaz_i;
    const double m = status[i] - cumhaz_i;
    martingale[i] = m;
    double dev = 0.0;
    if (status[i] > 0.0) {
      const double arg = std::max(1e-12, status[i] - m);
      dev = -2.0 * (m + status[i] * std::log(arg));
    } else {
      dev = -2.0 * m;
    }
    dev = std::sqrt(std::max(0.0, dev));
    if (m < 0.0) {
      dev = -dev;
    }
    deviance[i] = dev;
  }
  
  DevianceResult result{deviance, martingale, cumhaz, eta};
  return result;
}

NumericVector compute_eta_matrix(const NumericMatrix& X,
                                 const NumericVector& coef,
                                 const NumericVector& center,
                                 const NumericVector& scale) {
  const std::size_t n = X.nrow();
  const std::size_t p = X.ncol();
  if (coef.size() != static_cast<R_xlen_t>(p)) {
    stop("Length of coef must match number of columns in X");
  }
  if (!center.isNULL() && center.size() != static_cast<R_xlen_t>(p)) {
    stop("Length of center must match number of columns in X");
  }
  if (!scale.isNULL() && scale.size() != static_cast<R_xlen_t>(p)) {
    stop("Length of scale must match number of columns in X");
  }
  NumericVector eta(n, 0.0);
  for (std::size_t j = 0; j < p; ++j) {
    const double coef_j = coef[j];
    if (coef_j == 0.0) {
      continue;
    }
    const double mean = center.isNULL() ? 0.0 : center[j];
    double sd = scale.isNULL() ? 1.0 : scale[j];
    if (sd <= 0.0 || !std::isfinite(sd)) {
      sd = 1.0;
    }
    for (std::size_t i = 0; i < n; ++i) {
      eta[i] += ((X(i, j) - mean) / sd) * coef_j;
    }
  }
  return eta;
}

NumericVector compute_eta_bigmatrix(SEXP xpMat,
                                    const NumericVector& coef,
                                    const NumericVector& center,
                                    const NumericVector& scale) {
  XPtr<BigMatrix> pMat(xpMat);
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  if (coef.size() != static_cast<R_xlen_t>(p)) {
    stop("Length of coef must match number of columns in X");
  }
  if (!center.isNULL() && center.size() != static_cast<R_xlen_t>(p)) {
    stop("Length of center must match number of columns in X");
  }
  if (!scale.isNULL() && scale.size() != static_cast<R_xlen_t>(p)) {
    stop("Length of scale must match number of columns in X");
  }
  MatrixAccessor<double> accessor(*pMat);
  NumericVector eta(n, 0.0);
  for (std::size_t j = 0; j < p; ++j) {
    const double coef_j = coef[j];
    if (coef_j == 0.0) {
      continue;
    }
    const double mean = center.isNULL() ? 0.0 : center[j];
    double sd = scale.isNULL() ? 1.0 : scale[j];
    if (sd <= 0.0 || !std::isfinite(sd)) {
      sd = 1.0;
    }
    double* column = accessor[j];
    for (std::size_t i = 0; i < n; ++i) {
      eta[i] += ((column[i] - mean) / sd) * coef_j;
    }
  }
  return eta;
}



struct CoxDevianceResult {
  std::vector<double> cumhaz;
  std::vector<double> martingale;
  std::vector<double> deviance;
};

std::vector<std::size_t> order_desc(const NumericVector& time) {
  const std::size_t n = static_cast<std::size_t>(time.size());
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), static_cast<std::size_t>(0));
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    if (time[a] == time[b]) {
      return a > b;
    }
    return time[a] > time[b];
  });
  return order;
}

CoxDevianceResult compute_deviance_impl(const NumericVector& time,
                                        const NumericVector& status,
                                        const NumericVector& weights) {
  const std::size_t n = static_cast<std::size_t>(time.size());
  if (status.size() != static_cast<int>(n) || weights.size() != static_cast<int>(n)) {
    stop("`time`, `status`, and `weights` must have the same length");
  }
  
  std::vector<double> cumhaz(n, 0.0);
  std::vector<double> martingale(n, 0.0);
  std::vector<double> deviance(n, 0.0);
  
  const std::vector<std::size_t> order = order_desc(time);
  
  double risk_sum = 0.0;
  double cumulative_hazard = 0.0;
  
  std::size_t idx = 0;
  while (idx < n) {
    const double current_time = time[order[idx]];
    double block_weight_sum = 0.0;
    double block_event_sum = 0.0;
    std::vector<std::size_t> block_indices;
    while (idx < n && time[order[idx]] == current_time) {
      const std::size_t obs = order[idx];
      const double w = weights[obs];
      block_weight_sum += w;
      if (status[obs] > 0.0) {
        block_event_sum += w;
      }
      block_indices.push_back(obs);
      ++idx;
    }
    
    const double risk_with_block = risk_sum + block_weight_sum;
    const double hazard_increment = (risk_with_block > 0.0) ? block_event_sum / risk_with_block : 0.0;
    cumulative_hazard += hazard_increment;
    for (const auto obs : block_indices) {
      cumhaz[obs] = cumulative_hazard;
    }
    risk_sum += block_weight_sum;
  }
  
  const double eps = 1e-12;
  for (std::size_t i = 0; i < n; ++i) {
    const double delta = status[i];
    const double w = weights[i];
    const double hz = cumhaz[i];
    const double mart = w * (delta - hz);
    martingale[i] = mart;
    if (delta > 0.0) {
      const double term = std::max(hz, eps);
      deviance[i] = std::copysign(std::sqrt(std::max(0.0, -2.0 * (mart + w * std::log(term)))), mart);
    } else {
      deviance[i] = std::copysign(std::sqrt(std::max(0.0, -2.0 * mart)), mart);
    }
  }
  
  return {cumhaz, martingale, deviance};
}

NumericVector normalise_weights(Nullable<NumericVector> weights, std::size_t n) {
  if (weights.isNotNull()) {
    NumericVector w(weights);
    if (static_cast<std::size_t>(w.size()) != n) {
      stop("`weights` must have the same length as the input vectors");
    }
    for (double value : w) {
      if (!std::isfinite(value) || value <= 0.0) {
        stop("`weights` must contain strictly positive finite numbers");
      }
    }
    return w;
  }
  NumericVector ones(n, 1.0);
  return ones;
}

} // namespace





// [[Rcpp::export(name = "deviance_residuals_cpp")]]
List deviance_residuals_cpp(NumericVector time,
                            NumericVector status,
                            NumericVector eta,
                            std::string method = "efron") {
  method = Rcpp::as<std::string>(CharacterVector::create(method)[0]);
  std::transform(method.begin(), method.end(), method.begin(), ::tolower);
  if (method != "efron" && method != "breslow") {
    stop("method must be 'efron' or 'breslow'");
  }
  DevianceResult result = compute_deviance_core(time, status, eta, method);
  return List::create(Named("deviance") = result.deviance,
                      Named("martingale") = result.martingale,
                      Named("cumhaz") = result.cumhaz,
                      Named("linear_predictor") = result.eta);
}

// [[Rcpp::export(name = "matrix_deviance_residuals_cpp")]]
List matrix_deviance_residuals_cpp(NumericMatrix X,
                                   NumericVector coef,
                                   NumericVector time,
                                   NumericVector status,
                                   Nullable<NumericVector> center = R_NilValue,
                                   Nullable<NumericVector> scale = R_NilValue,
                                   std::string method = "efron") {
  NumericVector center_vec = center.isNotNull() ? NumericVector(center) : NumericVector();
  NumericVector scale_vec = scale.isNotNull() ? NumericVector(scale) : NumericVector();
  NumericVector eta = compute_eta_matrix(X, coef, center_vec, scale_vec);
  DevianceResult result = compute_deviance_core(time, status, eta, method);
  return List::create(Named("deviance") = result.deviance,
                      Named("martingale") = result.martingale,
                      Named("cumhaz") = result.cumhaz,
                      Named("linear_predictor") = result.eta);
}

// [[Rcpp::export(name = "big_deviance_residuals_cpp")]]
List big_deviance_residuals_cpp(SEXP xpMat,
                                NumericVector coef,
                                NumericVector time,
                                NumericVector status,
                                Nullable<NumericVector> center = R_NilValue,
                                Nullable<NumericVector> scale = R_NilValue,
                                std::string method = "efron") {
  NumericVector center_vec = center.isNotNull() ? NumericVector(center) : NumericVector();
  NumericVector scale_vec = scale.isNotNull() ? NumericVector(scale) : NumericVector();
  NumericVector eta = compute_eta_bigmatrix(xpMat, coef, center_vec, scale_vec);
  DevianceResult result = compute_deviance_core(time, status, eta, method);
  return List::create(Named("deviance") = result.deviance,
                      Named("martingale") = result.martingale,
                      Named("cumhaz") = result.cumhaz,
                      Named("linear_predictor") = result.eta);
}


// [[Rcpp::export]]
NumericVector cox_deviance_residuals_cpp(NumericVector time,
                                         NumericVector status,
                                         Nullable<NumericVector> weights = R_NilValue) {
  if (time.size() == 0) {
    stop("`time` must have positive length");
  }
  NumericVector w = normalise_weights(weights, static_cast<std::size_t>(time.size()));
  CoxDevianceResult res = compute_deviance_impl(time, status, w);
  NumericVector out(res.deviance.begin(), res.deviance.end());
  return out;
}

// [[Rcpp::export]]
List cox_deviance_details_cpp(NumericVector time,
                              NumericVector status,
                              Nullable<NumericVector> weights = R_NilValue) {
  if (time.size() == 0) {
    stop("`time` must have positive length");
  }
  NumericVector w = normalise_weights(weights, static_cast<std::size_t>(time.size()));
  CoxDevianceResult res = compute_deviance_impl(time, status, w);
  return List::create(
    Named("cumulative_hazard") = NumericVector(res.cumhaz.begin(), res.cumhaz.end()),
    Named("martingale") = NumericVector(res.martingale.begin(), res.martingale.end()),
    Named("deviance") = NumericVector(res.deviance.begin(), res.deviance.end())
  );
}

// [[Rcpp::export]]
NumericVector cox_deviance_residuals_big_cpp(SEXP xpMat,
                                             int time_col,
                                             int status_col,
                                             Nullable<NumericVector> weights = R_NilValue) {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  if (time_col < 1 || status_col < 1 ||
      time_col > static_cast<int>(pMat->ncol()) ||
      status_col > static_cast<int>(pMat->ncol())) {
    stop("Column indices out of range");
  }
  const std::size_t n = pMat->nrow();
  MatrixAccessor<double> accessor(*pMat);
  NumericVector time(n), status(n);
  double* time_ptr = accessor[time_col - 1];
  double* status_ptr = accessor[status_col - 1];
  for (std::size_t i = 0; i < n; ++i) {
    time[i] = time_ptr[i];
    status[i] = status_ptr[i];
  }
  NumericVector w = normalise_weights(weights, n);
  CoxDevianceResult res = compute_deviance_impl(time, status, w);
  return NumericVector(res.deviance.begin(), res.deviance.end());
}

// [[Rcpp::export]]
List cox_partial_deviance_big_cpp(SEXP xpMat,
                                  NumericVector coef,
                                  NumericVector time,
                                  NumericVector status) {
  Rcpp::XPtr<BigMatrix> pMat(xpMat);
  const std::size_t n = pMat->nrow();
  const std::size_t p = pMat->ncol();
  if (coef.size() != static_cast<int>(p)) {
    stop("Length of `coef` must match number of columns in the design matrix");
  }
  if (time.size() != static_cast<int>(n) || status.size() != static_cast<int>(n)) {
    stop("`time` and `status` must have length equal to the number of rows of the design matrix");
  }
  
  MatrixAccessor<double> accessor(*pMat);
  NumericVector eta(n);
  for (std::size_t j = 0; j < p; ++j) {
    double beta = coef[j];
    if (beta == 0.0) {
      continue;
    }
    double* col = accessor[j];
    for (std::size_t i = 0; i < n; ++i) {
      eta[i] += beta * col[i];
    }
  }
  
  std::vector<std::size_t> order = order_desc(time);
  double risk_sum = 0.0;
  double loglik = 0.0;
  
  std::size_t pos = 0;
  while (pos < n) {
    const double current_time = time[order[pos]];
    double block_sum_exp = 0.0;
    double block_event_sum_eta = 0.0;
    double block_event_count = 0.0;
    std::vector<std::size_t> block_indices;
    while (pos < n && time[order[pos]] == current_time) {
      const std::size_t idx = order[pos];
      const double e = std::exp(eta[idx]);
      block_sum_exp += e;
      if (status[idx] > 0.0) {
        block_event_sum_eta += eta[idx];
        block_event_count += status[idx];
      }
      block_indices.push_back(idx);
      ++pos;
    }
    const double risk_with_block = risk_sum + block_sum_exp;
    if (block_event_count > 0.0) {
      if (risk_with_block <= 0.0) {
        stop("Numerical issue: risk set sum non-positive");
      }
      loglik += block_event_sum_eta - block_event_count * std::log(risk_with_block);
    }
    risk_sum += block_sum_exp;
  }
  
  const double deviance = -2.0 * loglik;
  return List::create(
    Named("loglik") = loglik,
    Named("deviance") = deviance,
    Named("linear_predictor") = eta
  );
}


