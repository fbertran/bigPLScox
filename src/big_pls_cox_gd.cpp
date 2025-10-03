#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <limits>

using namespace Rcpp;

namespace {

arma::mat extract_matrix(const Rcpp::XPtr<BigMatrix>& mat_ptr) {
  const std::size_t n = mat_ptr->nrow();
  const std::size_t p = mat_ptr->ncol();
  arma::mat result(n, p);
  if (mat_ptr->matrix_type() != 8) {
    throw std::runtime_error("big.matrix must be of type double");
  }
  MatrixAccessor<double> accessor(*mat_ptr);
  for (std::size_t j = 0; j < p; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      result(i, j) = accessor[j][i];
    }
  }
  return result;
}

arma::uvec order_desc(const arma::vec& time) {
  arma::uvec idx = arma::sort_index(time, "descend");
  return idx;
}

double compute_loglik(const arma::mat& X, const arma::vec& status,
                      const arma::vec& eta) {
  const std::size_t n = X.n_rows;
  arma::vec exp_eta = arma::exp(eta);
  arma::vec cum_exp = arma::cumsum(exp_eta);
  arma::vec running_sum = arma::zeros<arma::vec>(X.n_cols);
  double loglik = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    running_sum += X.row(i).t() * exp_eta[i];
    if (status[i] > 0.5) {
      double denom = cum_exp[i];
      if (denom <= 0) {
        throw std::runtime_error("Numerical issue: risk set sum is non-positive");
      }
      loglik += eta[i] - std::log(denom);
    }
  }
  return loglik;
}

arma::vec compute_gradient(const arma::mat& X, const arma::vec& status,
                           const arma::vec& eta) {
  const std::size_t n = X.n_rows;
  const std::size_t p = X.n_cols;
  arma::vec grad(p, arma::fill::zeros);
  arma::vec exp_eta = arma::exp(eta);
  arma::vec cum_exp = arma::cumsum(exp_eta);
  arma::vec running_sum = arma::zeros<arma::vec>(p);
  for (std::size_t i = 0; i < n; ++i) {
    running_sum += X.row(i).t() * exp_eta[i];
    if (status[i] > 0.5) {
      const double denom = cum_exp[i];
      grad += X.row(i).t();
      grad -= running_sum / denom;
    }
  }
  return grad;
}

} // namespace

// [[Rcpp::export(name = "big_pls_cox_gd_cpp")]]
Rcpp::List big_pls_cox_gd_cpp(SEXP X_ptr, Rcpp::NumericVector time,
                              Rcpp::NumericVector status, int ncomp,
                              int max_iter, double tol, double learning_rate) {
  if (max_iter <= 0) {
    Rcpp::stop("`max_iter` must be positive");
  }
  if (tol <= 0) {
    Rcpp::stop("`tol` must be positive");
  }
  if (learning_rate <= 0) {
    Rcpp::stop("`learning_rate` must be positive");
  }
  
  Rcpp::XPtr<BigMatrix> mat_ptr(X_ptr);
  arma::mat Xfull = extract_matrix(mat_ptr);
  const std::size_t n = Xfull.n_rows;
  const std::size_t p = Xfull.n_cols;
  if (n == 0 || p == 0) {
    Rcpp::stop("`X` must have positive dimensions");
  }
  if (time.size() != static_cast<int>(n) || status.size() != static_cast<int>(n)) {
    Rcpp::stop("Length of `time` and `status` must equal number of rows of `X`");
  }
  
  if (ncomp < 1 || ncomp > static_cast<int>(p)) {
    Rcpp::stop("`ncomp` must be between 1 and ncol(X)");
  }
  const arma::uword k = static_cast<arma::uword>(ncomp);
  
  arma::vec time_vec(time.begin(), time.size(), false);
  arma::vec status_vec(status.begin(), status.size(), false);
  
  arma::uvec ord = order_desc(time_vec);
  arma::mat X = Xfull.rows(ord);
  arma::vec status_ord = status_vec.elem(ord);
  
  arma::mat Xuse = X.cols(0, k - 1);
  
  arma::vec beta = arma::zeros<arma::vec>(k);
  double prev_loglik = -std::numeric_limits<double>::infinity();
  bool converged = false;
  arma::vec grad(k);
  double loglik = prev_loglik;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::vec eta = Xuse * beta;
    grad = compute_gradient(Xuse, status_ord, eta);
    arma::vec beta_new = beta + learning_rate * grad;
    loglik = compute_loglik(Xuse, status_ord, Xuse * beta_new);
    if (iter > 0 && std::abs(loglik - prev_loglik) < tol) {
      converged = true;
      beta = beta_new;
      prev_loglik = loglik;
      return Rcpp::List::create(Rcpp::Named("coefficients") = beta,
                                Rcpp::Named("loglik") = loglik,
                                Rcpp::Named("iterations") = iter + 1,
                                Rcpp::Named("converged") = converged);
    }
    if (arma::norm(beta_new - beta, 2) < tol) {
      beta = beta_new;
      loglik = compute_loglik(Xuse, status_ord, Xuse * beta);
      converged = true;
      prev_loglik = loglik;
      return Rcpp::List::create(Rcpp::Named("coefficients") = beta,
                                Rcpp::Named("loglik") = loglik,
                                Rcpp::Named("iterations") = iter + 1,
                                Rcpp::Named("converged") = converged);
    }
    beta = beta_new;
    prev_loglik = loglik;
  }
  arma::vec eta_final = Xuse * beta;
  loglik = compute_loglik(Xuse, status_ord, eta_final);
  return Rcpp::List::create(Rcpp::Named("coefficients") = beta,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("iterations") = max_iter,
                            Rcpp::Named("converged") = converged);
}
