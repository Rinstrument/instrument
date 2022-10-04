#include <RcppArmadillo.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppArmadilloExtensions/sample.h>
#include <cmath> /* erf */
#include "min.h"

// [[Rcpp::export]]
double gzero(double a) {
  if(a > 0.0) {
    return a;
  } else {
    return 0.01;
  }
}

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
double log_unifd_dx(double x, double l, double u) {
  if(x < l | x > u) { return 0.0; }
  double res = 1.0 / (u - l);
  return std::log(res);
}

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

typedef double (*funcPtr)(arma::vec &, arma::mat &, int);
// [[Rcpp::export]]
Rcpp::XPtr<funcPtr> postFunc(std::string option) {
 if (option == "lp2pl") {
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&lp_2pl_logit_reg)));
  }
  else if (option == "lp2plr") {
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&lp_2pl_logit_reg)));
  }
  else {
    return Rcpp::XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
  }
}

// [[Rcpp::export]]
arma::vec amc(arma::mat & x, arma::vec x_start, int iter, int burn, int greedy_iterations, double a, arma::mat & data,
              SEXP logPostPtr, arma::vec & accept, arma::vec validation_indexes, arma::vec validation_lower, int p_reg) {
  int p = x.n_rows;
  double l_def = 2.38 * 2.38 / p;
  arma::vec x_current = x_start;
  arma::vec prop_mu_t = arma::vec(p, arma::fill::zeros);
  arma::vec prop_sigma_t = arma::vec(p, arma::fill::ones);
  arma::ivec p_update_index = arma::randi(iter, arma::distr_param(1, +(p)));
  arma::vec gam_correct = arma::vec(p, arma::fill::zeros);
  arma::vec gam_correct_iter = arma::vec(p, arma::fill::zeros);
  arma::vec gam_correct_iter_post_burn = arma::vec(p, arma::fill::zeros);
  arma::vec l_scaling_t = arma::vec(p, arma::fill::value(l_def));
  // random normal draws for proposal distribution
  arma::vec ru_prop_mh = arma::randu(iter, arma::distr_param(0, 1));
  arma::vec x_proposal = x_current;
  double R = 1;
  Rcpp::XPtr<funcPtr> logPost(logPostPtr);
  funcPtr lp = *logPost;
  // double (* log_likelihood)(arma::vec &, arma::mat &, int);
  // if(log_likelihood_selector == 0) { // 1-P logit
  //   log_likelihood = lp_2pl_logit_reg;
  // } else if(log_likelihood_selector == 1) { // 2-P logit
  //   log_likelihood = lp_2pl_logit_reg;
  // }
  for(int it = 0; it < iter; it++) {
    arma::uvec updates = arma::randperm(p, p_update_index(it));
    // // 1. Sample candidate value for component p_i
    x_proposal(updates) = l_scaling_t(updates) % (arma::randn(p_update_index(it), arma::distr_param(0, 1)) % prop_sigma_t(updates)) + x_current(updates);
    x_proposal = validate_proposal(x_proposal, validation_indexes, validation_lower);
    if(it >= greedy_iterations) {
      R = min(std::exp(lp(x_proposal, data, p_reg) - lp(x_current, data, p_reg)), 1.0) > ru_prop_mh(it);
    } else {
      R = min(std::exp(lp(x_proposal, data, p_reg) - lp(x_current, data, p_reg)), 1.0) > 0.234;
    }
    // 2. Select the value for x(t+1) according to MH criterion
    if(R == 1) {
      x_current(updates) = x_proposal(updates);
      accept(it) = 1.0;
    } else {
      x_proposal(updates) = x_current(updates);
      accept(it) = 0.0;
    }
    // 3. Adaptation step: update proposal distribution variance in two
    // steps and update component scaling parameter
    arma::vec x_current_diff_prop_mu = arma::vec(p_update_index(it), arma::fill::zeros);
    x_current_diff_prop_mu = x_current(updates) - prop_mu_t(updates);
    if(it >= burn) {
      for(int update_j = 0; update_j < p_update_index(it); update_j++) {
        int to_update = updates(update_j);
        x(to_update, gam_correct_iter_post_burn(to_update)) = x_current(to_update);
        gam_correct_iter_post_burn(to_update) += 1.0;
      }
    }
    if(it >= greedy_iterations) {
      gam_correct(updates) = 1.0 / (gam_correct_iter(updates) + 1.0);
      gam_correct_iter(updates) += 1.0;
      prop_mu_t(updates) = prop_mu_t(updates) + (gam_correct(updates) % x_current_diff_prop_mu);
      prop_sigma_t(updates) = prop_sigma_t(updates) + gam_correct(updates) % ((x_current_diff_prop_mu % x_current_diff_prop_mu) - prop_sigma_t(updates));
      l_scaling_t(updates) = arma::exp(arma::log(l_scaling_t(updates)) + (gam_correct(updates)) * (R - a));
    }
  }
  return gam_correct_iter_post_burn;
}
