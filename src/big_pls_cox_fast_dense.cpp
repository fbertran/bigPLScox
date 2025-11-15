#include <RcppArmadillo.h>
#include "fast_common.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
SEXP big_pls_cox_fast_dense_cpp(
    const arma::mat& X,
    const arma::vec& time,
    const arma::vec& status,
    int ncomp,
    const arma::vec& means,
    const arma::vec& sds,
    const arma::ivec& keepX)
{
  return run_pls_cox_dense(X, time, status, ncomp, means, sds, keepX);
}