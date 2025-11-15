#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>
#include "fast_common.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, bigmemory)]]

// [[Rcpp::export]]
SEXP big_pls_cox_fast_big_cpp(
    SEXP xpMat,
    const arma::vec& time,
    const arma::vec& status,
    int ncomp,
    const arma::vec& means,
    const arma::vec& sds,
    const arma::ivec& keepX)
{
  return run_pls_cox_big(xpMat, time, status, ncomp, means, sds, keepX);
}