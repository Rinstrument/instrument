#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec validate_proposal(arma::vec proposal, arma::vec indexes, arma::vec lower) {
  int n = indexes.n_rows;
  for(int i = 0; i < n; i++) {
    if(proposal(indexes(i)) < lower(i)) {
      proposal(indexes(i)) = 0.0001;
    }
  }
  return proposal;
}