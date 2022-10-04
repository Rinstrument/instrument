#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]

// [[Rcpp::export]]
double min(double a, double b) {
  if(std::isinf(a) | std::isinf(b)) { return 1.0; }
  return 0.5 * (a + b - std::abs(a - b));
}
