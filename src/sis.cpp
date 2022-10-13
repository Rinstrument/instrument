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
  // Draw from posterior distribution by sequential update with resampling
  // using the sequential importance sampling with resampling angorithm 
  // Reference: 
  // Arguments: 
  //  data - input matrix of observations
  //  n - number of resampling iterations
  // Value:
  //  vector of weighted estimates for each parameter over each of the n resampling iterations
  double eps = 1e-16;               // epsilon for log()
  int d_nrow = data.n_rows;         // number of rows in data
  int d_ncol = data.n_cols;         // number of cols in data
  int n_theta = d_nrow;             // number of theta parameters (theta = latent ability)
  int n_delta = d_ncol;             // number of deltas, difficult parameters
  int n_param = n_theta + n_delta;  // total number of parameters to estimate
  // pre-allocation of memory for model estimation
  arma::ivec seeds = arma::randi<arma::ivec>(n_param*n, arma::distr_param(+(1), +(999999999))); // random seeds for every draw of every parameter
  arma::vec x = arma::vec(n_param*n, arma::fill::zeros);           // store all draws
  arma::vec x_weighted = arma::vec(n_param*n, arma::fill::zeros);  // store all weighted draws
  arma::vec w = arma::vec(n, arma::fill::ones);                    // store current weights
  arma::vec w_history = arma::vec(n_param*n, arma::fill::ones);    // store all weights
  arma::vec u = arma::vec(n, arma::fill::zeros);                   // store updates for current iteration
  // pre-allocation of useful values
  int t = 1;                 // t indexes parameter during SIS iterations
  double theta_prob = 0.0;   // likelihood propability
  double n_eff = 0.0;        // ???
  // loop over each parameter
  //   for each parameter, resample n times
  //   update likelihood (weights), save history, check for degeneracy
  //   if weights are near degenerate, resample parameters and reset weights
  for(int p_ = 0; p_ < n_param; p_++) {   // p_ indexes current parameter to update
    for(int n_ = 0; n_ < n; n_++) {         // n_ indexes the resampling iterations (sample each "p_" "n_" times)
      if(t > n_theta) {                   // if t > n_theta, then we sample difficuty prior
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*5.0;
      } else {                            // if t <= n_theta, sample theta prior
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
      }
        double lk = 0.0;  // pre-allocate likelihood memory
        if(t > n_theta) {
          for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step
            int d_current_index = t - n_theta;
            theta_prob = 1.0 / (1.0 + std::exp(-(x((d_current_index-1)*n + n_) - x((t-1)*n + n_))));
            lk += data(i_, d_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, d_current_index - 1))*std::log(1.0 - theta_prob + eps);
          }
        } else {
          for(int j_ = 0; j_ < d_ncol; j_++) {  // compute the likelihood for the given step
            int theta_current_index = t;
            theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_current_index-1)*n + n_))));
            lk += data(theta_current_index - 1, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_current_index - 1, j_))*std::log(1.0 - theta_prob + eps);
          }
        }
        u(n_) = lk;             // update using the likelihood at "prior distribution" draw (log scale)
        w(n_) = w(n_) + u(n_);  // update weights on the log scale
    }
    // This section transforms the log weights by exponentiating them.
    // We cannot just exponentiate the weights, so we take the 
    // difference between weight and max weight,
    // and then exponentiate the difference minus logadd.
    // The weight adjustment method was taken from the SIR method of LaplacesDemon R package
    double md = arma::max(w);
    arma::vec lw = w - md;
    arma::vec probs = arma::exp(lw - logadd(lw));
    // End weight scale change
    // Update weight history
    for(int w_ = 0; w_ < n; w_++) {
      w_history((t-1)*n + w_) = probs(w_);
    }
    double norm_const = arma::accu(probs);   // sum of probs should be 1 since probs are probabilities
    // weight observations by their importance. This is the step by which the prior draws
    // are weighted by the likelihood, yielding draws from the joint posterior distribution.
    for(int n_ = 0; n_ < n; n_++) { 
      x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
    }
    n_eff = arma::accu(exp(w) % exp(w));
    Rcpp::Rcout << "n_eff: " << n_eff << std::endl;
    if(n_eff < 2.0e-10) {
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
