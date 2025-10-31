#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <numeric>

using namespace Rcpp;

namespace {

inline bool is_event(double status_value) {
  return status_value > 0.0;
}

inline bool almost_equal(double a, double b) {
  const double diff = std::fabs(a - b);
  const double scale = std::max(std::fabs(a), std::fabs(b));
  const double eps = std::numeric_limits<double>::epsilon();
  return diff <= eps * (1.0 + scale);
}

inline double regularise_positive(double value, double reference,
                                  const char* message) {
  if (value > 0.0 && std::isfinite(value)) {
    return value;
  }
  if (std::isfinite(reference) && reference > 0.0) {
    const double tol = std::max(1e-12, reference * 1e-12);
    if (std::isfinite(value) && value >= -tol) {
      return std::max(tol, std::numeric_limits<double>::min());
    }
  } else if (std::isfinite(value) && value >= -1e-12) {
    return std::max(1e-12, std::numeric_limits<double>::min());
  }
  stop(message);
}

} // anonymous namespace

// [[Rcpp::export(name = "cox_partial_loglik_cpp")]]
SEXP cox_partial_loglik_cpp(NumericMatrix X,
                            NumericVector time,
                            NumericVector status,
                            NumericMatrix beta,
                            std::string method,
                            bool return_all) {
  const int n = X.nrow();
  const int p = X.ncol();
  if (time.size() != n || status.size() != n) {
    stop("Length of 'time' and 'status' must match the number of rows in 'x'.");
  }
  if (beta.nrow() != p) {
    stop("Number of rows in 'b' must match the number of columns in 'x'.");
  }
  const int k = beta.ncol();
  if (k == 0) {
    NumericVector loglik_empty;
    if (!return_all) {
      return loglik_empty;
    }
    NumericMatrix w_empty(n, 0);
    NumericMatrix eta_empty(n, 0);
    NumericMatrix dmat_empty(n, 0);
    IntegerVector oo(n);
    for (int i = 0; i < n; ++i) {
      oo[i] = i + 1;
    }
    return List::create(Named("loglik") = loglik_empty,
                        Named("w") = w_empty,
                        Named("eta") = eta_empty,
                        Named("dmat") = dmat_empty,
                        Named("oo") = oo);
  }
  
  bool method_is_efron;
  if (method == "breslow") {
    method_is_efron = false;
  } else if (method == "efron") {
    method_is_efron = true;
  } else {
    stop("Unknown partial likelihood method: %s", method);
  }
  
  std::vector<int> order_idx(n);
  std::iota(order_idx.begin(), order_idx.end(), 0);
  std::stable_sort(order_idx.begin(), order_idx.end(), [&](int lhs, int rhs) {
    const double time_lhs = time[lhs];
    const double time_rhs = time[rhs];
    if (time_lhs < time_rhs - std::numeric_limits<double>::epsilon()) {
      return true;
    }
    if (time_lhs > time_rhs + std::numeric_limits<double>::epsilon()) {
      return false;
    }
    const double status_lhs = status[lhs];
    const double status_rhs = status[rhs];
    if (status_lhs > status_rhs) {
      return true;
    }
    if (status_lhs < status_rhs) {
      return false;
    }
    return lhs < rhs;
  });
  
  NumericVector time_sorted(n);
  NumericVector status_sorted(n);
  for (int i = 0; i < n; ++i) {
    const int idx = order_idx[i];
    time_sorted[i] = time[idx];
    status_sorted[i] = status[idx];
  }
  
  std::vector<int> group_start;
  std::vector< std::vector<int> > events_by_group;
  std::vector<int> event_counts;
  NumericVector rept(n);
  
  int pos = 0;
  while (pos < n) {
    const int start = pos;
    const double current_time = time_sorted[pos];
    while (pos < n && almost_equal(time_sorted[pos], current_time)) {
      ++pos;
    }
    group_start.push_back(start);
    std::vector<int> group_events;
    for (int i = start; i < pos; ++i) {
      if (is_event(status_sorted[i])) {
        group_events.push_back(i);
      }
    }
    const int event_count = static_cast<int>(group_events.size());
    event_counts.push_back(event_count);
    int counter = event_count;
    for (int idx : group_events) {
      rept[idx] = counter;
      --counter;
    }
    events_by_group.push_back(group_events);
  }
  
  std::vector<int> event_rows;
  std::vector<int> event_col_for_row(n, -1);
  int event_col = 0;
  for (std::size_t g = 0; g < events_by_group.size(); ++g) {
    for (int row : events_by_group[g]) {
      event_rows.push_back(row);
      event_col_for_row[row] = event_col;
      ++event_col;
    }
  }
  const int m = static_cast<int>(event_rows.size());
  if (m == 0) {
    stop("No complete observation. Failed to compute partial likelihood.");
  }
  
  NumericMatrix eta_orig(n, k);
  for (int col = 0; col < k; ++col) {
    NumericVector beta_col = beta(_, col);
    for (int row = 0; row < n; ++row) {
      double sum = 0.0;
      for (int j = 0; j < p; ++j) {
        sum += X(row, j) * beta_col[j];
      }
      eta_orig(row, col) = sum;
    }
  }
  
  NumericVector loglik(k);
  NumericMatrix w_last;
  NumericMatrix dmat_last;
  NumericVector wsum_last;
  if (return_all) {
    w_last = NumericMatrix(n, m);
    dmat_last = NumericMatrix(n, m);
    wsum_last = NumericVector(m);
  }
  
  for (int col = 0; col < k; ++col) {
    std::vector<double> eta_sorted(n);
    std::vector<double> suffix_scale(n);
    std::vector<double> suffix_sum(n);
    double current_scale = -std::numeric_limits<double>::infinity();
    double current_sum = 0.0;
    for (int i = n - 1; i >= 0; --i) {
      const int idx = order_idx[i];
      const double eta_val = eta_orig(idx, col);
      eta_sorted[i] = eta_val;
      if (!std::isfinite(current_scale)) {
        current_scale = eta_val;
        current_sum = 1.0;
      } else {
        const double new_scale = std::max(current_scale, eta_val);
        const double scaled_existing = current_sum * std::exp(current_scale - new_scale);
        const double scaled_new = std::exp(eta_val - new_scale);
        current_scale = new_scale;
        current_sum = scaled_existing + scaled_new;
      }
      suffix_scale[i] = current_scale;
      suffix_sum[i] = current_sum;
    }
    
    double loglik_col = 0.0;
    for (std::size_t g = 0; g < events_by_group.size(); ++g) {
      const int start = group_start[g];
      const int event_count = event_counts[g];
      if (event_count == 0) {
        continue;
      }
      const double risk_scale = suffix_scale[start];
      const double risk_sum_scaled = suffix_sum[start];
      double event_sum_eta = 0.0;
      double event_sum_scaled = 0.0;
      const std::vector<int>& group_events = events_by_group[g];
      for (int row : group_events) {
        event_sum_eta += eta_sorted[row];
        event_sum_scaled += std::exp(eta_sorted[row] - risk_scale);
      }
      loglik_col += event_sum_eta;
      for (int idx_ev = 0; idx_ev < event_count; ++idx_ev) {
        const int row = group_events[idx_ev];
        const double di = static_cast<double>(event_count);
        const double adjust = method_is_efron ? (di - rept[row]) / di : 0.0;
        double denom_scaled = risk_sum_scaled - adjust * event_sum_scaled;
        denom_scaled = regularise_positive(denom_scaled, risk_sum_scaled,
                                           "Non-positive risk set sum encountered while computing partial log-likelihood.");
        loglik_col -= risk_scale;
        loglik_col -= std::log(denom_scaled);
        if (return_all && col == k - 1) {
          const int event_col_idx = event_col_for_row[row];
          for (int i = start; i < n; ++i) {
            dmat_last(i, event_col_idx) = 1.0;
          }
          if (method_is_efron && adjust != 0.0) {
            for (int event_row : group_events) {
              dmat_last(event_row, event_col_idx) -= adjust;
            }
          }
          double wsum = 0.0;
          for (int i = start; i < n; ++i) {
            const double val = dmat_last(i, event_col_idx) * std::exp(eta_sorted[i] - risk_scale);
            w_last(i, event_col_idx) = val;
            wsum += val;
          }
          wsum_last[event_col_idx] = wsum;
        }
      }
    }
    loglik[col] = loglik_col;
  }
  
  if (!return_all) {
    return loglik;
  }
  
  for (int col_idx = 0; col_idx < m; ++col_idx) {
    const double wsum_raw = wsum_last[col_idx];
    const double wsum = regularise_positive(wsum_raw, wsum_raw,
                                            "Non-positive risk set sum encountered while scaling weights.");
    for (int i = 0; i < n; ++i) {
      w_last(i, col_idx) /= wsum;
    }
  }
  
  IntegerVector oo(n);
  for (int i = 0; i < n; ++i) {
    oo[i] = order_idx[i] + 1;
  }
  
  return List::create(Named("loglik") = loglik,
                      Named("w") = w_last,
                      Named("eta") = eta_orig,
                      Named("dmat") = dmat_last,
                      Named("oo") = oo);
}
