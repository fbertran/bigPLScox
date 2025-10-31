#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <limits>
#include <cmath>
#include <string>
#include <algorithm>

using namespace Rcpp;

namespace {

arma::Mat<double> map_bigmatrix(const Rcpp::XPtr<BigMatrix>& mat_ptr) {
  if (mat_ptr->matrix_type() != 8) {
    throw std::runtime_error("big.matrix must be of type double");
  }
  const std::size_t n = mat_ptr->nrow();
  const std::size_t p = mat_ptr->ncol();
  double* ptr = static_cast<double*>(mat_ptr->matrix());
  if (ptr == nullptr) {
    throw std::runtime_error("Null pointer encountered when mapping big.matrix");
  }
  return arma::Mat<double>(ptr, n, p, /*copy_aux_mem =*/ false, /*strict =*/ true);
}

struct PLSDecomposition {
  arma::mat scores;
  arma::mat loadings;
  arma::mat weights;
  arma::vec center;
  arma::vec scale;
};

arma::vec compute_column_means(const arma::mat& X) {
  return arma::mean(X, 0).t();
}

arma::vec compute_column_sds(const arma::mat& X, const arma::vec& means) {
  const std::size_t n = X.n_rows;
  const std::size_t p = X.n_cols;
  arma::vec sds(p, arma::fill::ones);
  for (std::size_t j = 0; j < p; ++j) {
    double sumsq = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double diff = X(i, j) - means[j];
      sumsq += diff * diff;
    }
    double var = 0.0;
    if (n > 1) {
      var = sumsq / static_cast<double>(n - 1);
    }
    if (!std::isfinite(var) || var <= 0.0) {
      sds[j] = 1.0;
    } else {
      sds[j] = std::sqrt(var);
    }
  }
  return sds;
}

PLSDecomposition compute_pls_components(const arma::mat& X,
                                        const Rcpp::NumericVector& time,
                                        const Rcpp::NumericVector& status,
                                        int ncomp,
                                        const Rcpp::IntegerVector& keepX) {
  const std::size_t n = X.n_rows;
  const std::size_t p = X.n_cols;
  if (n == 0 || p == 0) {
    Rcpp::stop("`X` must have positive dimensions");
  }

  arma::vec means = compute_column_means(X);
  arma::vec sds = compute_column_sds(X, means);

  arma::mat centered = X;
  for (std::size_t j = 0; j < p; ++j) {
    centered.col(j) -= means[j];
    centered.col(j) /= sds[j];
  }

  arma::mat deflated = centered;
  arma::mat scores(n, ncomp, arma::fill::zeros);
  arma::mat loadings(p, ncomp, arma::fill::zeros);
  arma::mat weights(p, ncomp, arma::fill::zeros);

  Rcpp::Environment survival_env = Rcpp::Environment::namespace_env("survival");
  Rcpp::Function coxph_fit = survival_env["coxph.fit"];
  Rcpp::Function coxph_control = survival_env["coxph.control"];
  Rcpp::List control = coxph_control();
  if (control.containsElementNamed("iter.max")) {
    int iter_max = Rcpp::as<int>(control["iter.max"]);
    if (iter_max < 50) {
      control["iter.max"] = 50;
    }
  } else {
    control["iter.max"] = 50;
  }

  Rcpp::NumericMatrix y(n, 2);
  for (std::size_t i = 0; i < n; ++i) {
    y(i, 0) = time[i];
    y(i, 1) = status[i];
  }
  Rcpp::NumericVector strata(n, 1.0);
  Rcpp::NumericVector offset(n, 0.0);
  Rcpp::NumericVector weights_vec(n, 1.0);
  Rcpp::CharacterVector rownms(n);
  for (std::size_t i = 0; i < n; ++i) {
    rownms[i] = std::to_string(i + 1);
  }

  for (int h = 0; h < ncomp; ++h) {
    Rcpp::NumericMatrix current_scores(n, h);
    for (int col = 0; col < h; ++col) {
      for (std::size_t row = 0; row < n; ++row) {
        current_scores(row, col) = scores(row, col);
      }
    }

    Rcpp::NumericVector init(h);
    for (int col = 0; col < h; ++col) {
      init[col] = 0.0;
    }

    Rcpp::List fit = coxph_fit(Rcpp::_["x"] = current_scores,
                               Rcpp::_["y"] = y,
                               Rcpp::_["strata"] = strata,
                               Rcpp::_["offset"] = offset,
                               Rcpp::_["init"] = init,
                               Rcpp::_["control"] = control,
                               Rcpp::_["weights"] = weights_vec,
                               Rcpp::_["method"] = "efron",
                               Rcpp::_["rownames"] = rownms);

    Rcpp::NumericVector residuals_r = fit["residuals"];
    arma::vec residuals = Rcpp::as<arma::vec>(residuals_r);

    arma::vec weight_vec = deflated.t() * residuals;
    if (keepX.size() > 0) {
      const int keep = (keepX.size() == 1) ? keepX[0] : keepX[h];
      if (keep > 0 && keep < static_cast<int>(p)) {
        arma::uvec order = arma::sort_index(arma::abs(weight_vec), "descend");
        arma::vec mask(p, arma::fill::zeros);
        for (int idx = 0; idx < keep; ++idx) {
          mask[order[idx]] = 1.0;
        }
        weight_vec %= mask;
      }
    }
    const double weight_norm_sq = arma::dot(weight_vec, weight_vec);
    if (weight_norm_sq <= 0.0 || !std::isfinite(weight_norm_sq)) {
      Rcpp::stop("Unable to compute weight vector; residuals may be zero");
    }
    const double weight_norm = std::sqrt(weight_norm_sq);
    weight_vec /= weight_norm;
    weights.col(h) = weight_vec;

    arma::vec score_vec = deflated * weight_vec;
    const double score_norm_sq = arma::dot(score_vec, score_vec);
    if (score_norm_sq <= 0.0 || !std::isfinite(score_norm_sq)) {
      Rcpp::stop("Computed score vector has zero variance");
    }
    scores.col(h) = score_vec;

    arma::vec loading_vec = (deflated.t() * score_vec) / score_norm_sq;
    loadings.col(h) = loading_vec;

    deflated -= score_vec * loading_vec.t();
  }

  PLSDecomposition result;
  result.scores = scores;
  result.loadings = loadings;
  result.weights = weights;
  result.center = means;
  result.scale = sds;
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
  arma::vec denom(n, arma::fill::zeros);
  double running_denom = 0.0;
  for (arma::uword idx = n; idx-- > 0;) {
    running_denom += exp_eta[idx];
    denom[idx] = running_denom;
  }

  double loglik = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    if (status[i] > 0.5) {
      const double risk_sum = denom[i];
      if (risk_sum <= 0) {
        throw std::runtime_error("Numerical issue: risk set sum is non-positive");
      }
      loglik += eta[i] - std::log(risk_sum);
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
  arma::vec denom(n, arma::fill::zeros);
  arma::mat cum_x(n, p, arma::fill::zeros);

  double running_denom = 0.0;
  arma::vec running_x(p, arma::fill::zeros);
  for (arma::uword idx = n; idx-- > 0;) {
    const double w = exp_eta[idx];
    running_denom += w;
    running_x += X.row(idx).t() * w;
    denom[idx] = running_denom;
    cum_x.row(idx) = running_x.t();
  }

  for (std::size_t i = 0; i < n; ++i) {
    if (status[i] > 0.5) {
      grad += X.row(i).t();
      grad -= cum_x.row(i).t() / denom[i];
    }
  }
  return grad;
}

} // namespace

// [[Rcpp::export(name = "big_pls_cox_gd_cpp")]]
Rcpp::List big_pls_cox_gd_cpp(SEXP X_ptr, Rcpp::NumericVector time,
                              Rcpp::NumericVector status, int ncomp,
                              int max_iter, double tol, double learning_rate,
                              Rcpp::IntegerVector keepX) {
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
  arma::Mat<double> Xfull = map_bigmatrix(mat_ptr);
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
  PLSDecomposition pls = compute_pls_components(Xfull, time, status, ncomp, keepX);

  arma::mat scores_ord = pls.scores.rows(ord);
  arma::vec status_ord = status_vec.elem(ord);

  arma::vec beta = arma::zeros<arma::vec>(k);
  double prev_loglik = -std::numeric_limits<double>::infinity();
  bool converged = false;
  arma::vec grad(k);
  double loglik = prev_loglik;

  for (int iter = 0; iter < max_iter; ++iter) {
    arma::vec eta = scores_ord * beta;
    grad = compute_gradient(scores_ord, status_ord, eta);
    arma::vec beta_new = beta + learning_rate * grad;
    loglik = compute_loglik(scores_ord, status_ord, scores_ord * beta_new);
    if (iter > 0 && std::abs(loglik - prev_loglik) < tol) {
      converged = true;
      beta = beta_new;
      prev_loglik = loglik;
      return Rcpp::List::create(
          Rcpp::Named("coefficients") = beta,
          Rcpp::Named("loglik") = loglik,
          Rcpp::Named("iterations") = iter + 1,
          Rcpp::Named("converged") = converged,
          Rcpp::Named("scores") = pls.scores,
          Rcpp::Named("loadings") = pls.loadings,
          Rcpp::Named("weights") = pls.weights,
          Rcpp::Named("center") = pls.center,
          Rcpp::Named("scale") = pls.scale);
    }
    if (arma::norm(beta_new - beta, 2) < tol) {
      beta = beta_new;
      loglik = compute_loglik(scores_ord, status_ord, scores_ord * beta);
      converged = true;
      prev_loglik = loglik;
      return Rcpp::List::create(
          Rcpp::Named("coefficients") = beta,
          Rcpp::Named("loglik") = loglik,
          Rcpp::Named("iterations") = iter + 1,
          Rcpp::Named("converged") = converged,
          Rcpp::Named("scores") = pls.scores,
          Rcpp::Named("loadings") = pls.loadings,
          Rcpp::Named("weights") = pls.weights,
          Rcpp::Named("center") = pls.center,
          Rcpp::Named("scale") = pls.scale);
    }
    beta = beta_new;
    prev_loglik = loglik;
  }
  arma::vec eta_final = scores_ord * beta;
  loglik = compute_loglik(scores_ord, status_ord, eta_final);
  return Rcpp::List::create(Rcpp::Named("coefficients") = beta,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("iterations") = max_iter,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("scores") = pls.scores,
                            Rcpp::Named("loadings") = pls.loadings,
                            Rcpp::Named("weights") = pls.weights,
                            Rcpp::Named("center") = pls.center,
                            Rcpp::Named("scale") = pls.scale);
}
