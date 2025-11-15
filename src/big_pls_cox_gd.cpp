#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include <limits>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// ----------------------------------------------------------------------------
// Zero-copy view of big.matrix using Armadillo
// ----------------------------------------------------------------------------
arma::mat map_bigmatrix(const Rcpp::XPtr<BigMatrix>& bm) {
  if (bm->matrix_type() != 8)
    Rcpp::stop("big.matrix must be double type");
  
  double* ptr = static_cast<double*>(bm->matrix());
  if (!ptr) Rcpp::stop("Null pointer in big.matrix");
  
  std::size_t n = bm->nrow();
  std::size_t p = bm->ncol();
  
  return arma::mat(ptr, n, p, /*copy_aux_mem =*/ false, /*strict =*/ true);
}

// ----------------------------------------------------------------------------
// Column means and sds
// ----------------------------------------------------------------------------
arma::vec col_means(const arma::mat& X) {
  return arma::mean(X, 0).t();
}

arma::vec col_sds(const arma::mat& X, const arma::vec& means) {
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  arma::vec sds(p);
  
  for (arma::uword j = 0; j < p; ++j) {
    double ss = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      double d = X(i,j) - means(j);
      ss += d*d;
    }
    double var = (n > 1 ? ss / double(n - 1) : 0.0);
    sds(j) = (var <= 0 || !std::isfinite(var)) ? 1.0 : std::sqrt(var);
  }
  return sds;
}

// ----------------------------------------------------------------------------
// PLS EXACT (Armardillo-matching)
// ----------------------------------------------------------------------------
struct PLSDecomposition {
  arma::mat scores;
  arma::mat loadings;
  arma::mat weights;
  arma::vec center;
  arma::vec scale;
};

PLSDecomposition compute_pls_components(
    const arma::mat& X,
    const NumericVector& time,
    const NumericVector& status,
    int ncomp,
    const IntegerVector& keepX)
{
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  
  arma::vec means = col_means(X);
  arma::vec sds   = col_sds(X, means);
  
  arma::mat Xs = X;
  for (arma::uword j = 0; j < p; ++j) {
    Xs.col(j) -= means(j);
    Xs.col(j) /= sds(j);
  }
  
  arma::mat deflated = Xs;
  arma::mat T(n, ncomp, arma::fill::zeros);
  arma::mat P(p, ncomp, arma::fill::zeros);
  arma::mat W(p, ncomp, arma::fill::zeros);
  
  // Prepare Cox residual engine
  Environment surv = Environment::namespace_env("survival");
  Function coxph_fit = surv["coxph.fit"];
  Function coxph_control = surv["coxph.control"];
  List control = coxph_control();
  control["iter.max"] = 50;
  
  NumericMatrix y(n, 2);
  for (arma::uword i = 0; i < n; ++i) {
    y(i,0) = time[i];
    y(i,1) = status[i];
  }
  
  NumericVector strata(n, 1.0);
  NumericVector offset(n, 0.0);
  NumericVector weights_vec(n, 1.0);
  CharacterVector rownames(n);
  for (arma::uword i=0; i<n; ++i)
    rownames[i] = std::to_string(i+1);
  
  arma::vec residuals(n);
  
  for (int h = 0; h < ncomp; ++h) {
    
    // -------------------------------------
    // Step 1: Compute residuals
    // -------------------------------------
    if (h == 0) {
      residuals = arma::vec(reinterpret_cast<double*>(const_cast<double*>(status.begin())),
                            n, true, true);
      residuals -= arma::mean(residuals);
      
    } else {
      NumericMatrix Tprev(n, h);
      for (int j = 0; j < h; ++j)
        for (arma::uword i=0; i<n; ++i)
          Tprev(i,j) = T(i,j);
      
      NumericVector init(h);
      for (int j=0;j<h;++j) init[j] = 0.0;
      
      List fit = coxph_fit(
        _["x"] = Tprev,
        _["y"] = y,
        _["strata"] = strata,
        _["offset"] = offset,
        _["init"] = init,
        _["control"] = control,
        _["weights"] = weights_vec,
        _["method"] = "efron",
        _["rownames"] = rownames
      );
      
      NumericVector r = fit["residuals"];
      residuals = arma::vec(reinterpret_cast<double*>(r.begin()), n, true, true);
      residuals -= arma::mean(residuals);
    }
    
    if (!residuals.is_finite() || arma::dot(residuals,residuals) < 1e-14) {
      residuals = arma::vec(reinterpret_cast<double*>(const_cast<double*>(status.begin())),
                            n, true, true);
      residuals -= arma::mean(residuals);
    }
    
    // -------------------------------------
    // Step 2: weight w = deflatedᵀ residuals
    // -------------------------------------
    arma::vec w = deflated.t() * residuals;
    
    if (keepX.size() > 0) {
      int keep = (keepX.size()==1 ? keepX[0] : keepX[h]);
      if (keep > 0 && keep < (int)p) {
        arma::uvec srt = arma::sort_index(arma::abs(w), "descend");
        arma::vec mask(p, arma::fill::zeros);
        for (int i = 0; i < keep; ++i) mask(srt(i)) = 1.0;
        w %= mask;
      }
    }
    
    double wnorm2 = arma::dot(w,w);
    if (wnorm2 <= 1e-20 || !std::isfinite(wnorm2)) {
      arma::uword best=0;
      double bestv=-1.0;
      for (arma::uword j=0;j<p;++j) {
        double v = std::abs(arma::dot(deflated.col(j), residuals));
        if (v>bestv) { bestv=v; best=j; }
      }
      w.zeros();
      w(best)=1.0;
      wnorm2=1.0;
    }
    
    w /= std::sqrt(wnorm2);
    W.col(h) = w;
    
    // -------------------------------------
    // Step 3: Score t = deflated w
    // -------------------------------------
    arma::vec t = deflated * w;
    
    // Armadillo PLS: center only, no variance normalization
    t -= arma::mean(t);
    
    // -------------------------------------
    // Step 4: Loading p = (deflatedᵀ t) / (tᵀ t)
    // -------------------------------------
    double tTt = arma::dot(t,t);
    if (tTt <= 1e-20) tTt = 1.0;
    
    arma::vec pvec = (deflated.t() * t) / tTt;
    
    T.col(h) = t;
    P.col(h) = pvec;
    
    // -------------------------------------
    // Step 5: Deterministic sign convention
    // -------------------------------------
    arma::uword jmax = arma::index_max(arma::abs(pvec));
    if (pvec(jmax) < 0.0) {
      T.col(h) *= -1.0;
      W.col(h) *= -1.0;
      P.col(h) *= -1.0;
      t       *= -1.0;
      pvec    *= -1.0;
    }
    
    // -------------------------------------
    // Step 6: Deflation
    // -------------------------------------
    deflated -= t * pvec.t();
  }
  
  PLSDecomposition res;
  res.scores = T;
  res.loadings = P;
  res.weights = W;
  res.center = means;
  res.scale = sds;
  return res;
}

// ----------------------------------------------------------------------------
// Cox loglik and gradient
// ----------------------------------------------------------------------------
double compute_loglik(const arma::mat& X, const arma::vec& status, const arma::vec& eta)
{
  arma::uword n = X.n_rows;
  arma::vec exp_eta = arma::exp(eta);
  
  arma::vec denom(n);
  double run = 0.0;
  for (arma::sword i = n - 1; i >= 0; --i) {
    run += exp_eta(i);
    denom(i) = run;
  }
  
  double ll = 0.0;
  for (arma::uword i=0;i<n;++i)
    if (status(i) > 0.5)
      ll += eta(i) - std::log(denom(i));
    
    return ll;
}

arma::vec compute_gradient(const arma::mat& X, const arma::vec& status, const arma::vec& eta)
{
  arma::uword n = X.n_rows;
  arma::uword k = X.n_cols;
  
  arma::vec grad(k, arma::fill::zeros);
  arma::vec exp_eta = arma::exp(eta);
  
  arma::vec denom(n);
  arma::mat cumX(n, k);
  
  double run = 0.0;
  arma::vec runx(k, arma::fill::zeros);
  
  for (arma::sword i = n - 1; i >= 0; --i) {
    double w = exp_eta(i);
    run  += w;
    runx += X.row(i).t() * w;
    denom(i) = run;
    cumX.row(i) = runx.t();
  }
  
  for (arma::uword i=0;i<n;++i) {
    if (status(i) > 0.5) {
      grad += X.row(i).t();
      grad -= cumX.row(i).t() / denom(i);
    }
  }
  
  return grad;
}

// ----------------------------------------------------------------------------
// Main function with stable backtracking line search
// ----------------------------------------------------------------------------

// [[Rcpp::export(name="big_pls_cox_gd_cpp")]]
List big_pls_cox_gd_cpp(SEXP X_ptr,
                        NumericVector time,
                        NumericVector status,
                        int ncomp,
                        int max_iter,
                        double learning_rate,
                        double tol,
                        IntegerVector keepX)
{
  Rcpp::XPtr<BigMatrix> xp(X_ptr);
  arma::mat X = map_bigmatrix(xp);
  
  arma::uword n = X.n_rows;
  if (time.size() != n || status.size() != n)
    stop("time/status length mismatch with X");
  
  // Extract EXACT PLS components
  PLSDecomposition pls = compute_pls_components(X, time, status, ncomp, keepX);
  
  // Order rows by decreasing time
  arma::vec time_v(reinterpret_cast<double*>(time.begin()), n, true, true);
  arma::vec status_v(reinterpret_cast<double*>(status.begin()), n, true, true);
  
  arma::uvec ord = arma::sort_index(time_v, "descend");
  arma::mat T_ord = pls.scores.rows(ord);
  arma::vec status_ord = status_v.elem(ord);
  
  // Gradient descent with backtracking
  arma::uword k = pls.scores.n_cols;
  arma::vec beta = arma::zeros(k);
  double prev_ll = -std::numeric_limits<double>::infinity();
  bool conv = false;
  
  std::vector<double> loglik_trace;
  std::vector<double> gradnorm_trace;
  std::vector<double> step_trace;
  std::vector< std::vector<double> > coef_trace;
  std::vector< std::vector<double> > eta_trace;
  
  for (int it = 0; it < max_iter; ++it) {
    
    arma::vec eta = T_ord * beta;
    arma::vec grad = compute_gradient(T_ord, status_ord, eta);
    
    double step = learning_rate;
    arma::vec beta_new;
    double ll_new = -std::numeric_limits<double>::infinity();
    
    while (step > 1e-20) {
      beta_new = beta + step * grad;
      ll_new = compute_loglik(T_ord, status_ord, T_ord * beta_new);
      if (ll_new > prev_ll || it == 0) break;
      step *= 0.5;
    }
    
    if (step <= 1e-20 || std::abs(ll_new - prev_ll) < tol) {
      conv = true;
      beta = beta_new;
      prev_ll = ll_new;
      break;
    }
    
    beta = beta_new;
    prev_ll = ll_new;
    
    loglik_trace.push_back(ll_new);
    gradnorm_trace.push_back(arma::norm(grad, 2));
    step_trace.push_back(step);
    
    coef_trace.push_back( arma::conv_to<std::vector<double>>::from(beta_new) );
    eta_trace.push_back( arma::conv_to<std::vector<double>>::from(T_ord * beta_new) );
  }
  
  return List::create(
    _["coefficients"] = beta,
    _["loglik"] = prev_ll,
    _["iterations"] = max_iter,
    _["converged"] = conv,
    _["scores"] = pls.scores,
    _["loadings"] = pls.loadings,
    _["weights"] = pls.weights,
    _["center"] = pls.center,
    _["scale"] = pls.scale,
    _["keepX"] = keepX,
    _["time"] = time,
    _["status"] = status,
    _["loglik_trace"] = loglik_trace,
    _["gradnorm_trace"] = gradnorm_trace,
    _["step_trace"] = step_trace,
    _["coef_trace"] = coef_trace,
    _["eta_trace"] = eta_trace
  );
}

