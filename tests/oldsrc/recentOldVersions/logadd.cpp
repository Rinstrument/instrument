#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double logadd(arma::vec x) {
  arma::vec x_sorted = sort(x, "descent");
//z <- x[1] + log(1 + sum(exp(x[-1] - x[1])))
  double x_1 = x_sorted(0);
  double z = x_1 + std::log(1.0 + arma::accu(arma::exp(x_sorted.rows(1, x_sorted.n_rows - 1) - x_1)));
  return z;
}