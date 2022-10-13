#include <RcppArmadillo.h>
#include "RcppArmadilloExtensions/sample.h"
#include <RcppArmadilloExtensions/fixprob.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "densities.h"
#include "random_distributions.h"
#include "logadd.h"

// [[Rcpp::export]]
arma::vec sis_theta_model(arma::mat & data, int n) {
  double eps = 1e-16;
  double alpha = 0.2;
  // int n = data.n_rows;
  int p = data.n_rows;
  int j = data.n_cols;
  int n_param = data.n_rows;
  arma::ivec seeds = arma::randi<arma::ivec>(p*n, arma::distr_param(+(1), +(999999999)));
  arma::vec x = arma::vec(p*n, arma::fill::zeros);
  arma::vec x_weighted = arma::vec(p*n, arma::fill::zeros);
  arma::vec w = arma::vec(n, arma::fill::ones);
  arma::vec w_history = arma::vec(p*n, arma::fill::ones);
  arma::vec u = arma::vec(n, arma::fill::zeros);
  int t = 1;
  double theta_prob = 0.0;
  double n_eff = 0.0;
  
  for(int p_ = 0; p_ < n_param; p_++) {
    for(int n_ = 0; n_ < n; n_++) {
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
        double lk = 0.0;
        double lp = 0.0;
        for(int j_ = 0; j_ < j; j_++) {
          theta_prob = 1.0 / (1.0 + std::exp(-(x((t-1)*n + n_))));
          lk += data(p_, j_)*std::log(theta_prob + eps) + (1.0 - data(p_, j_))*std::log(1.0 - theta_prob + eps);
        }
        // lp = log_normd_dx(x((t-1)*n_param + n_), 0.0, 1.0);
        u(n_) = lk;
        w(n_) = w(n_) + u(n_);
    }
    double md = arma::max(w);
    arma::vec lw = w - md;
    arma::vec probs = arma::exp(lw - logadd(lw));
    for(int w_ = 0; w_ < n; w_++) {
      w_history((t-1)*n + w_) = probs(w_);
    }
    // Rcpp::Rcout << ""
    // Rcpp::Rcout << "prob1: " << probs(0) << std::endl;
    // Rcpp::Rcout << "prob20: " << probs(19) << std::endl;
    // Rcpp::Rcout << "prob200: " << probs(199) << std::endl;
    // Rcpp::Rcout << "prob50: " << probs(49) << std::endl;
    // Rcpp::Rcout << "prob400: " << probs(399) << std::endl;
    double norm_const = arma::accu(probs);
    // Rcpp::Rcout << "norm_const: " << norm_const << std::endl;

    for(int n_ = 0; n_ < n; n_++) { 
      x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
    }
    //Rcpp::Rcout << "w: " << w << std::endl;
    n_eff = arma::accu(exp(w) % exp(w)); // 1.0 / 
    Rcpp::Rcout << "n_eff: " << n_eff << std::endl;

    // Rcpp::Rcout << "n_eff: " << n_eff << std::endl;
    if(n_eff < 2.0e-10) { //<  alpha*n
        Rcpp::Rcout << "triggered: " << std::endl;
        // // arma::vec w_temp = log(w);
        // double md = arma::max(w);
        // arma::vec lw = w - md;
        // arma::vec probs = arma::exp(lw - logadd(lw));
        arma::ivec sindex = arma::regspace<arma::ivec>(0, t*n - 1);
        arma::vec prev_probs = arma::vec(t*n, arma::fill::zeros);
        for(int w_ = 0; w_ < t*n; w_++) {
          prev_probs(w_) = w_history(w_);
        }
        arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs);
        for(int x_t = 0; x_t < t*n; x_t++) {
            x(x_t) = x(samp_indices(x_t));
        }
        for(int n_ = 0; n_ < n; n_++) {
            w(n_) = 1.0 / n;
        }
    }
    t += 1;
  }
  return x_weighted;
}













// [[Rcpp::export]]
arma::vec sis_theta_model2(arma::mat & data, int n) {
  double eps = 1e-16;
  double alpha = 0.2;
  int j = data.n_cols;
  int n_theta = data.n_rows;
  int p = j + n_theta;
  int n_param = p;
  arma::ivec seeds = arma::randi<arma::ivec>(p*n, arma::distr_param(+(1), +(999999999)));
  // x stores thetas (the latent vectors)
  arma::vec x = arma::vec(n_theta*n, arma::fill::zeros);
  arma::vec x_weighted = arma::vec(n_theta*n, arma::fill::zeros);
  // y stores deltas (item difficulties)
  arma::vec y = arma::vec(j*n, arma::fill::zeros);
  arma::vec y_weighted = arma::vec(j*n, arma::fill::zeros);

  arma::vec w = arma::vec(n, arma::fill::ones);
  arma::vec w_history = arma::vec(p*n, arma::fill::ones);

  arma::vec u = arma::vec(n, arma::fill::zeros);
  int t = 1;
  double theta_prob = 0.0;
  double n_eff = 0.0;
  
  for(int p_ = 0; p_ < n_param; p_++) {
    for(int n_ = 0; n_ < n; n_++) {
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
        double lk = 0.0;
        double lp = 0.0;
        for(int j_ = 0; j_ < j; j_++) {
          theta_prob = 1.0 / (1.0 + std::exp(-(x((t-1)*n + n_))));
          lk += data(p_, j_)*std::log(theta_prob + eps) + (1.0 - data(p_, j_))*std::log(1.0 - theta_prob + eps);
        }
        // lp = log_normd_dx(x((t-1)*n_param + n_), 0.0, 1.0);
        u(n_) = lk;
        w(n_) = w(n_) + u(n_);
    }
    double md = arma::max(w);
    arma::vec lw = w - md;
    arma::vec probs = arma::exp(lw - logadd(lw));
    for(int w_ = 0; w_ < n; w_++) {
      w_history((t-1)*n + w_) = probs(w_);
    }
    // Rcpp::Rcout << ""
    // Rcpp::Rcout << "prob1: " << probs(0) << std::endl;
    // Rcpp::Rcout << "prob20: " << probs(19) << std::endl;
    // Rcpp::Rcout << "prob200: " << probs(199) << std::endl;
    // Rcpp::Rcout << "prob50: " << probs(49) << std::endl;
    // Rcpp::Rcout << "prob400: " << probs(399) << std::endl;
    double norm_const = arma::accu(probs);
    // Rcpp::Rcout << "norm_const: " << norm_const << std::endl;

    for(int n_ = 0; n_ < n; n_++) { 
      x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
    }
    //Rcpp::Rcout << "w: " << w << std::endl;
    n_eff = arma::accu(exp(w) % exp(w)); // 1.0 / 
    Rcpp::Rcout << "n_eff: " << n_eff << std::endl;

    // Rcpp::Rcout << "n_eff: " << n_eff << std::endl;
    if(n_eff < 2.0e-30) { //<  alpha*n
        Rcpp::Rcout << "triggered: " << std::endl;
        // // arma::vec w_temp = log(w);
        // double md = arma::max(w);
        // arma::vec lw = w - md;
        // arma::vec probs = arma::exp(lw - logadd(lw));
        arma::ivec sindex = arma::regspace<arma::ivec>(0, t*n - 1);
        arma::vec prev_probs = arma::vec(t*n, arma::fill::zeros);
        for(int w_ = 0; w_ < t*n; w_++) {
          prev_probs(w_) = w_history(w_);
        }
        arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs);
        for(int x_t = 0; x_t < t*n; x_t++) {
            x(x_t) = x(samp_indices(x_t));
        }
        for(int n_ = 0; n_ < n; n_++) {
            w(n_) = 1.0 / n;
        }
    }
    t += 1;
  }
  return x_weighted;
}