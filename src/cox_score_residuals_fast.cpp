// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// -------------------------------------------------------------
// Cox score residuals (fast, O(n), fully in C++)
// r_i = d_i - exp(eta_i) / sum_{j in R_i} exp(eta_j)
// -------------------------------------------------------------

// [[Rcpp::export]]
arma::vec cox_score_residuals_fast(
    const arma::vec &time,
    const arma::vec &status,
    const arma::vec &eta
) {
  int n = time.n_elem;
  
  // Sort indices decreasing by time
  arma::uvec ord = arma::sort_index(time, "descend");
  
  arma::vec exp_eta = arma::exp(eta);
  arma::vec r(n);
  double risk_sum = 0.0;
  
  for (int k = 0; k < n; k++) {
    int i = ord[k];
    risk_sum += exp_eta[i];
    double haz = exp_eta[i] / risk_sum;
    r[i] = status[i] - haz;
  }
  
  // Center for numerical stability
  r -= arma::mean(r);
  
  return r;
}
