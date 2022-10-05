#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
// arma::vec validate_proposal(arma::vec proposal, arma::vec indexes, arma::vec lower) {
//   int n = indexes.n_rows;
//   for(int i = 0; i < n; i++) {
//     if(proposal(indexes(i)) < lower(i)) {
//       proposal(indexes(i)) = 0.0001;
//     }
//   }
//   return proposal;
// }

// [[Rcpp::export]]
arma::vec validate_proposal(arma::vec proposal, arma::vec indexes, arma::vec lower, arma::vec upper) {
  int n = indexes.n_rows;
  for(int i = 0; i < n; i++) {
    if(proposal(indexes(i)) < lower(i)) {
      proposal(indexes(i)) = lower(i);
    }
    if(proposal(indexes(i)) > upper(i)) {
      proposal(indexes(i)) = upper(i);
    }
  }
  return proposal;
}