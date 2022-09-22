#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::export]]
double min(double a, double b) {
  return 0.5 * (a + b - std::abs(a - b));
}

// [[Rcpp::export]]
arma::uvec get(std::list<arma::uvec> list, int l){
  std::list<arma::uvec>::iterator it = list.begin();
  for(int i = 0; i < l; i++){
    it++;
  }
  return *it;
}

// [[Rcpp::export]]
double sum_list(std::list<arma::uvec> a) {
  double b = 0;
  for(int j = 0; j < 3; j++) {
    b += arma::accu(get(a, j));
  }
  return b;
}

// [[Rcpp::export]]
int invoke(int x, int y,
           int (*func)(int*, int*))
{
  return func(x, y);
}

// // [[Rcpp::export]]
// double evaluate_log_likelihood(arma::vec x, arma::mat data, std::list<arma::uvec> parameter_indexes, double(*log_likelihood)(arma::vec, arma::mat, std::list<arma::uvec>)) {
//   return log_likelihood(x, data, parameter_indexes);
// }


/*** R
min(2, 5.5)
get(list(1:4, 5:7, 9:20), 0)
get(list(1:4, 5:7, 9:20), 1)
get(list(1:4, 5:7, 9:20), 2)

sum_list(list(1:4, 5:7, 9:20))
sum(unlist(list(1:4, 5:7, 9:20)))


*/
