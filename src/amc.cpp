#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath> /* erf */
#include "min.h"
#include "validate_proposal.h"
#include "densities.h"
#include "posteriors.h"

// [[Rcpp::export]]
SEXP amc(arma::mat & x, arma::vec x_start, int iter, int burn, int greedy_iterations, double a, arma::mat & data,
         int lp_select, arma::vec & accept, arma::vec validation_indexes, arma::vec validation_lower, arma::vec validation_upper,
         arma::vec & gam_correct_iter_post_burn, int p_reg) {
  int p = x.n_rows;
  double l_def = 2.38 * 2.38 / p;
  arma::vec x_current = x_start;
  arma::vec prop_mu_t = arma::vec(p, arma::fill::zeros);
  arma::vec prop_sigma_t = arma::vec(p, arma::fill::ones);
  arma::ivec p_update_index = arma::randi(iter, arma::distr_param(1, +(p)));
  arma::vec gam_correct = arma::vec(p, arma::fill::zeros);
  arma::vec gam_correct_iter = arma::vec(p, arma::fill::zeros);
  //arma::vec gam_correct_iter_post_burn = arma::vec(p, arma::fill::zeros);
  arma::vec l_scaling_t = arma::vec(p, arma::fill::value(l_def));
  // random normal draws for proposal distribution
  arma::vec ru_prop_mh = arma::randu(iter, arma::distr_param(0, 1));
  arma::vec x_proposal = x_current;
  double R = 1;
  double (* lp)(arma::vec &, arma::mat &, int);
  if(lp_select == 0) { // 1-P logit
    lp = lp_lm; //lp_2pl_logit_reg;
  } else if(lp_select == 1) { // 2-P logit
    lp = lp_lm; //lp_2pl_logit_reg;
  }
  for(int it = 0; it < iter; it++) {
    arma::uvec updates = arma::randperm(p, p_update_index(it));
    // // 1. Sample candidate value for component p_i
    x_proposal(updates) = l_scaling_t(updates) % (arma::randn(p_update_index(it), arma::distr_param(0, 1)) % prop_sigma_t(updates)) + x_current(updates);
    x_proposal = validate_proposal(x_proposal, validation_indexes, validation_lower, validation_upper);
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
  return R_NilValue;
}
