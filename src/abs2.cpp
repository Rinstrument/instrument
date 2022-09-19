#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]

// [[Rcpp::export]]
double abs2(double x) {
  double nx = 0.001;
  double nx_upper = 9.999;
  if(x < 0.0) {
    return nx;
  } else if(x >= 10) {
    return nx_upper;
  } else {
    return x;
  }
}