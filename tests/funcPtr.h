#include <RcppArmadillo.h>
//https://stackoverflow.com/questions/43616778/passing-user-created-c-functions-in-rcpp-as-arguments
typedef double (*funcPtr)(arma::vec &, arma::mat &, int);
typedef Rcpp::XPtr<funcPtr> XPtrFuncPtr_t;