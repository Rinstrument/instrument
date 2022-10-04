#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "densities.h"

// [[Rcpp::export]]
double lp_2pl_logit(arma::vec & x, arma::mat & data) {
  double eps = 1e-16;
  int n = data.n_rows;
  int j = data.n_cols;
  double theta_prob = 0.0;
  double lp = 0.0;
  for(int jj = 0; jj < j; jj++) {
    for(int i = 0; i < n; i++) {
      theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*j))) - x(jj + j))));
      lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
      if(jj == 0) {
        lp += log_normd_dx(x(i + (2*j)), 0.0, 1.0);
      }
    }
    lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
    lp += log_normd_dx(x(jj + j), 0.0, 5.0);
  }
  return lp;
}

// [[Rcpp::export]]
double lp_2pl_logit_reg(arma::vec & x, arma::mat & data, int p) {
  int n = data.n_rows;
  int dj = data.n_cols;
  int j = dj - p;
  double lp = 0.0;
  for(int i = j; i < dj; i++) {
    for(int nn = 0; nn < n; nn++) {
      x(2*j + nn) += x(2*j + n + i - j) * data(nn,i);
    }
    lp += log_normd_dx(x(2*j + n + i - j), 0.0, 1.0);
  }
  arma::vec x_ = x.subvec(0, 2*j + n - 1);
  arma::mat data_ = data.cols(0, j - 1);
  lp += lp_2pl_logit(x_, data_);
  return lp;
}
