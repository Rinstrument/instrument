#include <RcppArmadillo.h>

#ifndef RANDOM_DISTRIBUTIONS
#define RANDOM_DISTRIBUTIONS

double r8poly_value_horner(int m, arma::vec c, double x);
double r8_uniform_01(int &seed);
double normal_01_sample(int & seed);
double normal_01_cdf(double x);
double normal_01_cdf_inv(double p);
double truncated_normal_ab_sample(double mu, double sigma, double a, double b, int &seed);
double normal_cdf_inv(double cdf, double mu, double sigma);
double log_normal_cdf_inv(double cdf, double mu, double sigma);
double log_normal_sample(double mu, double sigma, int &seed);

#endif