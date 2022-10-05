#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]

// [[Rcpp::export]]
double lognorm_dens_dx(double x, double mean, double sd) {
  double res = 0.0;
  res = (1.0 / (sd * x * std::sqrt(2 * 3.14159265358979323846))) * std::exp(-(pow((std::log(x) - mean), 2.0)) / (2.0 * pow(sd, 2.0)));
  if(std::isnan(res)) { res = 0.0; }
  return res;
}

// [[Rcpp::export]]
double logtruncnorm_dens_dx(double x, double mean, double sd, double a, double b) {
  if(x < a | x > b) { return 0.0; }
  double st_norm = (1.0 / std::sqrt(2 * 3.14159265358979323846)) * std::exp(-0.5 * pow((x - mean) / sd, 2.0));
  double phi_b = 0.5 * (1.0 + erf(((b - mean) / sd) / std::sqrt(2)));
  double phi_a = 0.5 * (1.0 + erf(((a - mean) / sd) / std::sqrt(2)));
  return std::log((1 / sd) * (st_norm / (phi_b - phi_a)));
  //return -(std::log(st_norm) - std::log(phi_b - phi_a));
}

// [[Rcpp::export]]
double log_normd_dx(double x, double mean, double sd) {
  return -std::log(sd) - 0.5 * pow((x - mean) / sd, 2.0);
}

// [[Rcpp::export]]
double normd_dx(double x, double mean, double sd) {
  return (1.0 / (sd * std::sqrt(2 * 3.14159265358979323846))) * std::exp(-0.5 * pow((x - mean) / sd, 2.0));
}

// [[Rcpp::export]]
double log_unifd_dx(double x, double l, double u) {
  if(x < l | x > u) { return 0.0; }
  double res = 1.0 / (u - l);
  return std::log(res);
}

// [[Rcpp::export]]
double unifd_dx(double x, double l, double u) {
  if(x < l | x > u) { return 0.0; }
  double res = 1.0 / (u - l);
  return res;
}