#include <RcppArmadillo.h>
#include "RcppArmadilloExtensions/sample.h"
#include <RcppArmadilloExtensions/fixprob.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "densities.h"
#include "random_distributions.h"
#include "logadd.h"

// // [[Rcpp::export]]
// arma::vec sis_theta_model(arma::mat & data, int n) {
//   double eps = 1e-16;
//   double alpha = 0.2;
//   // int n = data.n_rows;
//   int p = data.n_rows;
//   int j = data.n_cols;
//   int n_param = data.n_rows;
//   arma::ivec seeds = arma::randi<arma::ivec>(p*n, arma::distr_param(+(1), +(999999999)));
//   arma::vec x = arma::vec(p*n, arma::fill::zeros);
//   arma::vec x_weighted = arma::vec(p*n, arma::fill::zeros);
//   arma::vec w = arma::vec(n, arma::fill::ones);
//   arma::vec w_history = arma::vec(p*n, arma::fill::ones);
//   arma::vec u = arma::vec(n, arma::fill::zeros);
//   int t = 1;
//   double theta_prob = 0.0;
//   double n_eff = 0.0;
  
//   for(int p_ = 0; p_ < n_param; p_++) {
//     for(int n_ = 0; n_ < n; n_++) {
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
//         double lk = 0.0;
//         double lp = 0.0;
//         for(int j_ = 0; j_ < j; j_++) {
//           theta_prob = 1.0 / (1.0 + std::exp(-(x((t-1)*n + n_))));
//           lk += data(p_, j_)*std::log(theta_prob + eps) + (1.0 - data(p_, j_))*std::log(1.0 - theta_prob + eps);
//         }
//         // lp = log_normd_dx(x((t-1)*n_param + n_), 0.0, 1.0);
//         u(n_) = lk;
//         w(n_) = w(n_) + u(n_);
//     }
//     double md = arma::max(w);
//     arma::vec lw = w - md;
//     arma::vec probs = arma::exp(lw - logadd(lw));
//     for(int w_ = 0; w_ < n; w_++) {
//       w_history((t-1)*n + w_) = probs(w_);
//     }
//     // Rcpp::Rcout << ""
//     // Rcpp::Rcout << "prob1: " << probs(0) << std::endl;
//     // Rcpp::Rcout << "prob20: " << probs(19) << std::endl;
//     // Rcpp::Rcout << "prob200: " << probs(199) << std::endl;
//     // Rcpp::Rcout << "prob50: " << probs(49) << std::endl;
//     // Rcpp::Rcout << "prob400: " << probs(399) << std::endl;
//     double norm_const = arma::accu(probs);
//     // Rcpp::Rcout << "norm_const: " << norm_const << std::endl;

//     for(int n_ = 0; n_ < n; n_++) { 
//       x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
//     }
//     //Rcpp::Rcout << "w: " << w << std::endl;
//     n_eff = arma::accu(exp(w) % exp(w)); // 1.0 / 
//     Rcpp::Rcout << "n_eff: " << n_eff << std::endl;

//     // Rcpp::Rcout << "n_eff: " << n_eff << std::endl;
//     if(n_eff < 2.0e-10) { //<  alpha*n
//         Rcpp::Rcout << "triggered: " << std::endl;
//         // // arma::vec w_temp = log(w);
//         // double md = arma::max(w);
//         // arma::vec lw = w - md;
//         // arma::vec probs = arma::exp(lw - logadd(lw));
//         arma::ivec sindex = arma::regspace<arma::ivec>(0, t*n - 1);
//         arma::vec prev_probs = arma::vec(t*n, arma::fill::zeros);
//         for(int w_ = 0; w_ < t*n; w_++) {
//           prev_probs(w_) = w_history(w_);
//         }
//         arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs);
//         for(int x_t = 0; x_t < t*n; x_t++) {
//             x(x_t) = x(samp_indices(x_t));
//         }
//         for(int n_ = 0; n_ < n; n_++) {
//             w(n_) = 1.0 / n;
//         }
//     }
//     t += 1;
//   }
//   return x_weighted;
// }



// // [[Rcpp::export]]
// arma::vec sis_theta_model2(arma::mat & data, int n) {
//   // Draw from posterior distribution by sequential update with resampling
//   // using the sequential importance sampling with resampling angorithm 
//   // Reference: 
//   // Arguments: 
//   //  data - input matrix of observations
//   //  n - number of resampling iterations
//   // Value:
//   //  vector of weighted estimates for each parameter over each of the n resampling iterations
//   double eps = 1e-16;               // epsilon for log()
//   int d_nrow = data.n_rows;         // number of rows in data
//   int d_ncol = data.n_cols;         // number of cols in data
//   int n_theta = d_nrow;             // number of theta parameters (theta = latent ability)
//   int n_delta = d_ncol;             // number of deltas, difficult parameters
//   int n_param = n_theta + n_delta;  // total number of parameters to estimate
//   // pre-allocation of memory for model estimation
//   arma::ivec seeds = arma::randi<arma::ivec>(n_param*n, arma::distr_param(+(1), +(999999999))); // random seeds for every draw of every parameter
//   arma::vec x = arma::vec(n_param*n, arma::fill::zeros);           // store all draws
//   arma::vec x_weighted = arma::vec(n_param*n, arma::fill::zeros);  // store all weighted draws
//   arma::vec w = arma::vec(n, arma::fill::ones);                    // store current weights
//   arma::vec w_history = arma::vec(n_param*n, arma::fill::ones);    // store all weights
//   arma::vec u = arma::vec(n, arma::fill::zeros);                   // store updates for current iteration
//   // pre-allocation of useful values
//   int t = 1;                 // t indexes parameter during SIS iterations
//   double theta_prob = 0.0;   // likelihood propability
//   double n_eff = 0.0;        // ???
//   // loop over each parameter
//   //   for each parameter, resample n times
//   //   update likelihood (weights), save history, check for degeneracy
//   //   if weights are near degenerate, resample parameters and reset weights
//   for(int p_ = 0; p_ < n_param; p_++) {   // p_ indexes current parameter to update
//     for(int n_ = 0; n_ < n; n_++) {         // n_ indexes the resampling iterations (sample each "p_" "n_" times)
//       if(t > n_theta) {                   // if t > n_theta, then we sample difficuty prior
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*2.0;
//       } else {                            // if t <= n_theta, sample theta prior
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
//       }
//       double lk = 0.0;  // pre-allocate likelihood memory
//       if(t > n_theta) {
//         for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step
//           int d_current_index = t - n_theta;
//           theta_prob = 1.0 / (1.0 + std::exp(-(x((i_)*n + n_) - x((t-1)*n + n_))));
//           lk += data(i_, d_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, d_current_index - 1))*std::log(1.0 - theta_prob + eps);
//         }
//       } else {
//         for(int j_ = 0; j_ < d_ncol; j_++) {  // compute the likelihood for the given step
//           int theta_current_index = t;
//           theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_current_index-1)*n + n_))));
//           lk += data(theta_current_index - 1, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_current_index - 1, j_))*std::log(1.0 - theta_prob + eps);
//         }
//       }
//       u(n_) = lk;             // update using the likelihood at "prior distribution" draw (log scale)
//       w(n_) = w(n_) + u(n_);  // update weights on the log scale
//     }
//     // This section transforms the log weights by exponentiating them.
//     // We cannot just exponentiate the weights, so we take the 
//     // difference between weight and max weight,
//     // and then exponentiate the difference minus logadd.
//     // The weight adjustment method was taken from the SIR method of LaplacesDemon R package
//     double md = arma::max(w);
//     arma::vec lw = w - md;
//     arma::vec probs = arma::exp(lw - logadd(lw));
//     // End weight scale change
//     // Update weight history
//     for(int w_ = 0; w_ < n; w_++) {
//       w_history((t-1)*n + w_) = probs(w_);
//     }
//     double norm_const = arma::accu(probs);   // sum of probs should be 1 since probs are probabilities
//     // weight observations by their importance. This is the step by which the prior draws
//     // are weighted by the likelihood, yielding draws from the joint posterior distribution.
//     for(int n_ = 0; n_ < n; n_++) { 
//       x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
//     }
//     // This is the stage where weights are checked for degeneracy (and rejuvinated).
//     // 1) In short this section does the following: compute some measure of "issues" with the weights.
//     // 2) If "issues" is sufficiently small (indicated by a cutoff value), then we know that the weights are
//     // degenerate or approaching degeneracy.
//     // 3) Fix the degeneracy by resampling all parameter draws with probability proportional to weigths and
//     //    reset all weights to 1 / n. This will be triggered a number of times during runtime to repair
//     //    weights that are overly concentrated near zero, and keep going. Without this step, the algorithm 
//     //    will not work.
//     n_eff = arma::accu(exp(w) % exp(w));              // check weights for issues using sum of squared weights
//     //Rcpp::Rcout << "n_eff: " << n_eff << std::endl;
//     if(n_eff < 2.0e-10) {   // if weights are under the minimum that we decide is tolerable, then do 2) and 3).
//         arma::ivec sindex = arma::regspace<arma::ivec>(0, t*n - 1);  // vector from 0, ..., t*n - 1 (zero up to current)
//         arma::vec prev_probs = arma::vec(t*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//         for(int w_ = 0; w_ < t*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//           prev_probs(w_) = w_history(w_);
//         }
//         // Now, we complete 3) by resampling using weight history
//         arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs);
//         for(int x_t = 0; x_t < t*n; x_t++) {   // update each parameter draw with the resampled values
//             x(x_t) = x(samp_indices(x_t));
//         }
//         for(int n_ = 0; n_ < n; n_++) {   // weight = 1.0 / n for each weight. This is a reset of weights with equal value per
//             w(n_) = 1.0 / n;              // weight. It's like starting the weighting over again.
//         }
//     }   // End of the weight rejuvination step.
//     t += 1;  // That was for one parameter (t = 1), increment t by one and repeat for the rest of the parameters
//   }
//   return x_weighted;  // return vector of weighted estimates (will have to sum across draws for the expected value (posterior mean))
// }

















// // [[Rcpp::export]]
// arma::vec sis3(arma::mat & data, int n, double tol) {
//   // Draw from posterior distribution by sequential update with resampling
//   // using the sequential importance sampling with resampling angorithm 
//   // Reference: 
//   // Arguments: 
//   //  data - input matrix of observations
//   //  n - number of resampling iterations
//   // Value:
//   //  vector of weighted estimates for each parameter over each of the n resampling iterations
//   double eps = 1e-16;               // epsilon for log()
//   int d_nrow = data.n_rows;         // number of rows in data
//   int d_ncol = data.n_cols;         // number of cols in data
//   int n_theta = d_nrow;             // number of theta parameters (theta = latent ability)
//   int n_delta = d_ncol;             // number of deltas, difficult parameters
//   int n_alpha = d_ncol;             // number of thetas, ability parameters
//   int n_param = n_theta + n_delta + n_alpha;  // total number of parameters to estimate
//   // pre-allocation of memory for model estimation
//   arma::ivec seeds = arma::randi<arma::ivec>(n_param*n, arma::distr_param(+(1), +(999999999))); // random seeds for every draw of every parameter
//   arma::vec x = arma::vec(n_param*n, arma::fill::zeros);           // store all draws
//   arma::vec x_weighted = arma::vec(n_param*n, arma::fill::zeros);  // store all weighted draws
//   arma::vec w = arma::vec(n, arma::fill::ones);                    // store current weights
//   arma::vec w_history = arma::vec(n_param*n, arma::fill::ones);    // store all weights
//   arma::vec u = arma::vec(n, arma::fill::zeros);                   // store updates for current iteration
//   // pre-allocation of useful values
//   int t = 1;                 // t indexes parameter during SIS iterations
//   double theta_prob = 0.0;   // likelihood propability
//   double n_eff = 0.0;        // ???
//   // loop over each parameter
//   //   for each parameter, resample n times
//   //   update likelihood (weights), save history, check for degeneracy
//   //   if weights are near degenerate, resample parameters and reset weights
//   for(int p_ = 0; p_ < n_param; p_++) {   // p_ indexes current parameter to update
//     for(int n_ = 0; n_ < n; n_++) {         // n_ indexes the resampling iterations (sample each "p_" "n_" times)
//       if(t <= n_theta) {                  // if t <= n_theta, sample theta prior
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
//       } else if(t <= n_theta + n_delta) { // if t > n_theta and t <= n_theta + n_delta, then we sample difficuty prior
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*2.0;
//       } else {                            // if t > n_theta + n_delta and t <= n_theta + n_delta + n_alpha, then we sample ability prior
//         //x((t-1)*n + n_) = truncated_normal_ab_sample(1.0, 2.0, 0.0, , seeds((t-1)*n + n_));
//         Rcpp::Rcout << "triggered: " << std::endl;
//         double a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0;
//         while (a_samp < 0.0) {
//           Rcpp::Rcout << "a_samp < 0.0: " << a_samp << std::endl;
//           a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0;
//         }
//         x((t-1)*n + n_) = a_samp;
//         Rcpp::Rcout << a_samp << std::endl;
//       }
//       double lk = 0.0;  // pre-allocate likelihood memory
//       // --------------------------
//       if(t <= n_theta) { 
//         for(int j_ = 0; j_ < d_ncol; j_++) {  // compute the likelihood for the given step (theta step)
//           int theta_current_index = t;
//           theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_current_index-1)*n + n_))));
//           lk += data(theta_current_index - 1, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_current_index - 1, j_))*std::log(1.0 - theta_prob + eps);
//         }
//       } else if(t <= n_theta + n_delta) { 
//         for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (delta step)
//           int d_current_index = t - n_theta;
//           theta_prob = 1.0 / (1.0 + std::exp(-(x((i_)*n + n_) - x((t-1)*n + n_))));
//           lk += data(i_, d_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, d_current_index - 1))*std::log(1.0 - theta_prob + eps);
//         }
//       } else {
//         for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (alpha step)
//           int a_current_index = t - n_theta - n_delta; //*x((i_)*n + n_)) - x((t-n_delta-1)*n + n_)
//           theta_prob = 1.0 / (1.0 + std::exp(-(   (x((t-1)*n + n_)*x((i_)*n + n_)  )  - x((t-n_delta-1)*n + n_)     ))); //t-n_delta-1???
//           lk += data(i_, a_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, a_current_index - 1))*std::log(1.0 - theta_prob + eps);
//         }
//         Rcpp::Rcout << "lk: " << lk << std::endl;
//       }
//       u(n_) = std::exp(lk);             // update using the likelihood at "prior distribution" draw (log scale)
//       w(n_) = w(n_) * u(n_);  // update weights on the log scale
//     }
//     // This section transforms the log weights by exponentiating them.
//     // We cannot just exponentiate the weights, so we take the 
//     // difference between weight and max weight,
//     // and then exponentiate the difference minus logadd.
//     // The weight adjustment method was taken from the SIR method of LaplacesDemon R package

//     // -----------------------------------------
//     // double md = arma::max(w);
//     // arma::vec lw = w - md;
//     // arma::vec probs = arma::exp(lw - logadd(lw));
//     // -----------------------------------------

//     // End weight scale change
//     // Update weight history
//     for(int w_ = 0; w_ < n; w_++) {
//       w_history((t-1)*n + w_) = probs(w_);
//     }
//     double norm_const = arma::accu(probs);   // sum of probs should be 1 since probs are probabilities
//     // weight observations by their importance. This is the step by which the prior draws
//     // are weighted by the likelihood, yielding draws from the joint posterior distribution.
//     for(int n_ = 0; n_ < n; n_++) { 
//       x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
//     }
//     // This is the stage where weights are checked for degeneracy (and rejuvinated).
//     // 1) In short this section does the following: compute some measure of "issues" with the weights.
//     // 2) If "issues" is sufficiently small (indicated by a cutoff value), then we know that the weights are
//     // degenerate or approaching degeneracy.
//     // 3) Fix the degeneracy by resampling all parameter draws with probability proportional to weigths and
//     //    reset all weights to 1 / n. This will be triggered a number of times during runtime to repair
//     //    weights that are overly concentrated near zero, and keep going. Without this step, the algorithm 
//     //    will not work.
//     n_eff = arma::accu(exp(w) % exp(w));              // check weights for issues using sum of squared weights
//     //Rcpp::Rcout << "n_eff: " << n_eff << std::endl;
//     if(n_eff < tol) {   // if weights are under the minimum that we decide is tolerable, then do 2) and 3).
//       // update each parameter draw with the resampled values
//       // since the alpha parameters must be positive, resample two separate subsets
//       // if(t <= n_theta + n_delta) {   // non-alpha parameters
//       //   arma::ivec sindex = arma::regspace<arma::ivec>(0, t*n - 1);  // vector from 0, ..., t*n - 1 (zero up to current)
//       //   arma::vec prev_probs = arma::vec(t*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//       //   for(int w_ = 0; w_ < t*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//       //     prev_probs(w_) = w_history(w_);
//       //   }
//       //   // Now, we complete 3) by resampling using weight history
//       //   arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs);
//       //   for(int x_t = 0; x_t < t*n; x_t++) {
//       //     x(x_t) = x(samp_indices(x_t));
//       //   }
//       // } else {                       // alpha parameters & non-alpha parameters
//       //   arma::ivec sindex = arma::regspace<arma::ivec>(0, (n_theta + n_delta)*n - 1);  // vector from 0, ..., (n_theta + n_delta)*n - 1 (zero up to theta, delta)
//       //   arma::ivec sindex2 = arma::regspace<arma::ivec>((n_theta + n_delta)*n, t*n - 1);  // vector from (n_theta + n_delta)*n, ..., t*n - 1 (zero up to current alpha)
//       //   arma::vec prev_probs = arma::vec((n_theta + n_delta)*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//       //   arma::vec prev_probs2 = arma::vec(t*n - (n_theta + n_delta)*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//       //   for(int w_ = 0; w_ < (n_theta + n_delta)*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//       //     prev_probs(w_) = w_history(w_);
//       //   }
//       //   prev_probs = prev_probs / arma::accu(prev_probs);
//       //   for(int w_ = (n_theta + n_delta)*n; w_ < t*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//       //     prev_probs2(w_ - (n_theta + n_delta)*n) = w_history(w_);
//       //   }
//       //   prev_probs2 = prev_probs2 / arma::accu(prev_probs2);
//       //   // Now, we complete 3) by resampling using weight history
//       //   arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs); 
//       //   arma::ivec samp_indices2 = Rcpp::RcppArmadillo::sample(sindex2, t*n, true, prev_probs2); 
//       //   for(int x_t = 0; x_t < (n_theta + n_delta)*n; x_t++) {
//       //     x(x_t) = x(samp_indices(x_t));
//       //   }
//       //   for(int x_t = (n_theta + n_delta)*n; x_t < t*n; x_t++) {
//       //     x(x_t) = x(samp_indices2(x_t));
//       //   }
//       // }

//       arma::ivec sindex = arma::regspace<arma::ivec>((t-1)*n, t*n - 1); //0  // vector from 0, ..., t*n - 1 (zero up to current)
//       arma::vec prev_probs = arma::vec(n, arma::fill::zeros);  //t*n  // get previous probabilities for the resampling step
//       for(int w_ = (t-1)*n; w_ < t*n; w_++) { //w_=0  // for every previous probability, use the weight which was used earlier 
//         prev_probs(w_ - (t-1)*n) = w_history(w_);
//       }
//       prev_probs = prev_probs / arma::accu(prev_probs);
//       // Now, we complete 3) by resampling using weight history
//       arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, n, true, prev_probs); // t*n
//       for(int x_t = (t-1)*n; x_t < t*n; x_t++) { //t*n
//         x(x_t) = x(samp_indices(x_t - (t-1)*n));
//         // if(samp_indices(x_t) >= (n_theta + n_delta)*n) {
//         //   x(x_t) = std::abs(x(x_t));
//         // }
//       }
//       for(int n_ = 0; n_ < n; n_++) {   // weight = 1.0 / n for each weight. This is a reset of weights with equal value per
//         w(n_) = 1.0 / n;              // weight. It's like starting the weighting over again.
//       }
//     }   // End of the weight rejuvination step.
//     t += 1;  // That was for one parameter (t = 1), increment t by one and repeat for the rest of the parameters
//   }
//   return x_weighted;  // return vector of weighted estimates (will have to sum across draws for the expected value (posterior mean))
// }








// // since the alpha parameters must be positive, resample two separate subsets
//       // if(t <= n_theta + n_delta) {   // non-alpha parameters
//       //   arma::ivec sindex = arma::regspace<arma::ivec>(0, t*n - 1);  // vector from 0, ..., t*n - 1 (zero up to current)
//       //   arma::vec prev_probs = arma::vec(t*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//       //   for(int w_ = 0; w_ < t*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//       //     prev_probs(w_) = w_history(w_);
//       //   }
//       //   // Now, we complete 3) by resampling using weight history
//       //   arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs);
//       //   for(int x_t = 0; x_t < t*n; x_t++) {
//       //     x(x_t) = x(samp_indices(x_t));
//       //   }
//       // } else {                       // alpha parameters & non-alpha parameters
//       //   arma::ivec sindex = arma::regspace<arma::ivec>(0, (n_theta + n_delta)*n - 1);  // vector from 0, ..., (n_theta + n_delta)*n - 1 (zero up to theta, delta)
//       //   arma::ivec sindex2 = arma::regspace<arma::ivec>((n_theta + n_delta)*n, t*n - 1);  // vector from (n_theta + n_delta)*n, ..., t*n - 1 (zero up to current alpha)
//       //   arma::vec prev_probs = arma::vec((n_theta + n_delta)*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//       //   arma::vec prev_probs2 = arma::vec(t*n - (n_theta + n_delta)*n, arma::fill::zeros);    // get previous probabilities for the resampling step
//       //   for(int w_ = 0; w_ < (n_theta + n_delta)*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//       //     prev_probs(w_) = w_history(w_);
//       //   }
//       //   prev_probs = prev_probs / arma::accu(prev_probs);
//       //   for(int w_ = (n_theta + n_delta)*n; w_ < t*n; w_++) {   // for every previous probability, use the weight which was used earlier 
//       //     prev_probs2(w_ - (n_theta + n_delta)*n) = w_history(w_);
//       //   }
//       //   prev_probs2 = prev_probs2 / arma::accu(prev_probs2);
//       //   // Now, we complete 3) by resampling using weight history
//       //   arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, t*n, true, prev_probs); 
//       //   arma::ivec samp_indices2 = Rcpp::RcppArmadillo::sample(sindex2, t*n, true, prev_probs2); 
//       //   for(int x_t = 0; x_t < (n_theta + n_delta)*n; x_t++) {
//       //     x(x_t) = x(samp_indices(x_t));
//       //   }
//       //   for(int x_t = (n_theta + n_delta)*n; x_t < t*n; x_t++) {
//       //     x(x_t) = x(samp_indices2(x_t));
//       //   }
//       // }













// [[Rcpp::export]]
arma::vec sis4(arma::mat & data, int n, double tol) {
  // Draw from posterior distribution by sequential update with resampling
  // using the sequential importance sampling with resampling angorithm 
  // Reference: 
  // Arguments: 
  //  data - input matrix of observations
  //  n - number of resampling iterations
  //  tol - tolerance for resampling. If sum of squared weight falls below tol,
  //        resampling is performed and weights reset.
  // Value:
  //  vector of weighted estimates for each parameter over each of the n resampling iterations
  double eps = 1e-16;               // epsilon for log()
  int d_nrow = data.n_rows;         // number of rows in data
  int d_ncol = data.n_cols;         // number of cols in data
  int n_theta = d_nrow;             // number of theta parameters (theta = latent ability)
  int n_delta = d_ncol;             // number of deltas, difficult parameters
  int n_alpha = d_ncol;             // number of thetas, ability parameters
  int n_param = n_theta + n_delta + n_alpha;  // total number of parameters to estimate
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
      if(t <= n_theta) {                  // if t <= n_theta, sample theta prior
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
      } else if(t <= n_theta + n_delta) { // if t > n_theta and t <= n_theta + n_delta, then we sample difficuty prior
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*2.0;
      } else {                            // if t > n_theta + n_delta and t <= n_theta + n_delta + n_alpha, then we sample ability prior
        //x((t-1)*n + n_) = truncated_normal_ab_sample(1.0, 2.0, 0.0, , seeds((t-1)*n + n_));
        double a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0;
        while (a_samp < 0.0) {
          a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0;
        }
        x((t-1)*n + n_) = a_samp;
      }
      double lk = 0.0;  // pre-allocate likelihood memory
      // --------------------------
      if(t <= n_theta) { 
        for(int j_ = 0; j_ < d_ncol; j_++) {  // compute the likelihood for the given step (theta step)
          int theta_current_index = t;
          theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_current_index-1)*n + n_))));
          lk += data(theta_current_index - 1, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_current_index - 1, j_))*std::log(1.0 - theta_prob + eps);
        }
      } else if(t <= n_theta + n_delta) { 
        for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (delta step)
          int d_current_index = t - n_theta;
          theta_prob = 1.0 / (1.0 + std::exp(-(x((i_)*n + n_) - x((t-1)*n + n_))));
          lk += data(i_, d_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, d_current_index - 1))*std::log(1.0 - theta_prob + eps);
        }
      } else {
        for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (alpha step)
          int a_current_index = t - n_theta - n_delta; //*x((i_)*n + n_)) - x((t-n_delta-1)*n + n_)
          theta_prob = 1.0 / (1.0 + std::exp(-(   (x((t-1)*n + n_)*x((i_)*n + n_)  )  - x((t-n_delta-1)*n + n_)     ))); //t-n_delta-1???
          lk += data(i_, a_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, a_current_index - 1))*std::log(1.0 - theta_prob + eps);
        }
      }
      u(n_) = lk;   // update using the likelihood at "prior distribution" draw (log scale)
      w(n_) = w(n_) + u(n_);  // update weights on the log scale
    }
    // This section transforms the log weights by exponentiating them.
    // We cannot just exponentiate the weights, so we take the 
    // difference between weight and max weight,
    // and then exponentiate the difference minus logadd.
    // The weight adjustment method was taken from the SIR method of LaplacesDemon R package
    // -----------------------------------------
    double md = arma::max(w);
    arma::vec lw = w - md;
    arma::vec probs = arma::exp(lw - logadd(lw));
    // -----------------------------------------
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
    // This is the stage where weights are checked for degeneracy (and rejuvinated).
    // 1) In short this section does the following: compute some measure of "issues" with the weights.
    // 2) If "issues" is sufficiently small (indicated by a cutoff value), then we know that the weights are
    // degenerate or approaching degeneracy.
    // 3) Fix the degeneracy by resampling all parameter draws with probability proportional to weigths and
    //    reset all weights to 1 / n. This will be triggered a number of times during runtime to repair
    //    weights that are overly concentrated near zero, and keep going. Without this step, the algorithm 
    //    will not work.
    n_eff = arma::accu(exp(w) % exp(w));     // check weights for issues using sum of squared weights
    if(n_eff < tol) {   // if weights are under the minimum that we decide is tolerable, then do 2) and 3).
      // update each parameter draw with the resampled values
      arma::ivec sindex = arma::regspace<arma::ivec>((t-1)*n, t*n - 1); //0  // vector from 0, ..., t*n - 1 (zero up to current)
      arma::vec prev_probs = arma::vec(n, arma::fill::zeros);  //t*n  // get previous probabilities for the resampling step
      for(int w_ = (t-1)*n; w_ < t*n; w_++) { //w_=0  // for every previous probability, use the weight which was used earlier 
        prev_probs(w_ - (t-1)*n) = w_history(w_);
      }
      prev_probs = prev_probs / arma::accu(prev_probs);
      // Now, we complete 3) by resampling using weight history
      arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, n, true, prev_probs); // t*n
      for(int x_t = (t-1)*n; x_t < t*n; x_t++) { //t*n
        x(x_t) = x(samp_indices(x_t - (t-1)*n));
      }
      for(int n_ = 0; n_ < n; n_++) {   // weight = 1.0 / n for each weight. This is a reset of weights with equal value per
        w(n_) = 1.0 / n;              // weight. It's like starting the weighting over again.
      }
    }   // End of the weight rejuvination step.
    t += 1;  // That was for one parameter (t = 1), increment t by one and repeat for the rest of the parameters
  }
  return x_weighted;  // return vector of weighted estimates (will have to sum across draws for the expected value (posterior mean))
}

// go back to this recent version, I don't think it was updating with the 
// exponentiated one because it had a bug and failed to compile!
// still works!
















// // [[Rcpp::export]]
// arma::vec sis5(arma::mat & data, int n_dimensions, arma::vec dimension_start, arma::vec dimension_end, int n, double tol) {
//   // Draw from posterior distribution by sequential update with resampling
//   // using the sequential importance sampling with resampling angorithm 
//   // Reference: 
//   // Arguments: 
//   //  data - input matrix of observations
//   //  n_dimensions - number of first-order dimensions
//   //  dimension_start - start index (from 1) of each dimension. The dimension ends at dimension_start(n+1) - 1
//   //  dimension_end - the end index of each dimension.
//   //  n - number of resampling iterations
//   //  tol - tolerance for resampling. If sum of squared weight falls below tol,
//   //        resampling is performed and weights reset.
//   // Value:
//   //  vector of weighted estimates for each parameter over each of the n resampling iterations
//   double eps = 1e-16;                       // epsilon for log() when performing likelihood calculation
//   int d_nrow = data.n_rows;                 // number of rows in data
//   int d_ncol = data.n_cols;                 // number of cols in data
//   // n_dimensions                           // number of latent (ability) dimensions
//   int n_theta = d_nrow*(n_dimensions + 1);  // number of theta parameters (theta = latent ability)
//   int n_delta = d_ncol;                     // number of deltas, difficult parameters
//   int n_alpha = d_ncol;                     // number of thetas, ability parameters
//   int n_lambda = n_dimensions;              // number of lambda parameters (loadings for second-order theta on each of the first-order thetas)
//   int n_param = n_theta + n_delta + n_alpha + n_lambda;  // total number of parameters to estimate
//   // pre-allocation of memory for model estimation
//   arma::ivec seeds = arma::randi<arma::ivec>(n_param*n, arma::distr_param(+(1), +(999999999))); // random seeds for every draw of every parameter
//   arma::vec x = arma::vec(n_param*n, arma::fill::zeros);           // store all draws
//   arma::vec x_weighted = arma::vec(n_param*n, arma::fill::zeros);  // store all weighted draws
//   arma::vec w = arma::vec(n, arma::fill::ones);                    // store current weights
//   arma::vec w_history = arma::vec(n_param*n, arma::fill::ones);    // store all weights
//   arma::vec u = arma::vec(n, arma::fill::zeros);                   // store updates for current iteration
//   // pre-allocation of useful values
//   int t = 1;                 // t indexes parameter during SIS iterations
//   double theta_prob = 0.0;   // likelihood propability
//   double n_eff = 0.0;        // pre allocation of memory for tolerance test
//   // loop over each parameter
//   //   for each parameter, resample n times
//   //   update likelihood (weights), save history, check for degeneracy
//   //   if weights are near degenerate, resample parameters and reset weights
//   for(int p_ = 0; p_ < n_param; p_++) {   // p_ indexes current parameter to update
//     for(int n_ = 0; n_ < n; n_++) {         // n_ indexes the resampling iterations (sample each "p_" "n_" times)
//       if(t <= n_theta) {                  // if t <= n_theta, sample theta prior
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;

//       } else if(t <= n_theta + n_delta) {                // if t > n_theta and t <= n_theta + n_delta, then we sample difficuty prior
//         x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*2.0;
//       } else if(t <= n_theta + n_delta + n_alpha) {      // if t > n_theta + n_delta and t <= n_theta + n_delta + n_alpha, then we sample ability prior
//         //x((t-1)*n + n_) = truncated_normal_ab_sample(1.0, 2.0, 0.0, , seeds((t-1)*n + n_));
//         double a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0;
//         while (a_samp < 0.0) {
//           a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0; // Probably switch this back to the truncated normal dist rather than inefficiency of over-sampling
//         }
//         x((t-1)*n + n_) = a_samp;
//       } else {                                           // if t > n_theta + n_delta + n_alpha and t <= n_theta + n_delta + n_alpha + n_lambda, then we sample lambdas (loadings for latent dimensions)
//         x((t-1)*n + n_) = truncated_normal_ab_sample(0.0, 2.0, -5.0, 5.0, seeds((t-1)*n + n_));
//       }
//       double lk = 0.0;  // pre-allocate likelihood memory
//       // --------------------------
//       if(t <= n_theta) {
//         // first, estimate the second order theta
//         if(t <= d_nrow) {
//           // for(int j_ = 0; j_ < d_ncol; j_++) {  // compute the likelihood for the given step (theta step, second-order dimension)
//           //   int theta_current_index = t;
//           //   theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_current_index-1)*n + n_))));
//           //   lk += data(theta_current_index - 1, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_current_index - 1, j_))*std::log(1.0 - theta_prob + eps);
//           // }
//         } else {    // second, loop over each first order dimension of theta
//           // calculate which dimension we are in
//           int current_dimension = (t-1) / d_nrow;      // 0, 1, ... zero is the general dimension, 1, 2, 3, ... first-order
//           for(int j_ = dimension_start(current_dimension - 1) - 1; j_ < dimension_end(current_dimension - 1) - 1; j_++) {  // compute the likelihood for the given step (theta step, all first-order dimensions)
//             int theta_current_index = t; //current_dimension*d_nrow + t - (current_dimension - 1)*d_nrow;
//             int theta_row_index = (theta_current_index - 1) % d_nrow;
//             theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_current_index-1)*n + n_))));
//             lk += data(theta_row_index, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_row_index, j_))*std::log(1.0 - theta_prob + eps);
//           }
//         }
//       } else if(t <= n_theta + n_delta) {
//         for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (delta step)
//           int d_current_index = t - n_theta;
//           // loop through dimension indexes to discover which dimension the given question t belongs to
//           // store the dimension d_ for calling up the right theta to condition the alpha likelihood on it
//           int theta_dim = 0; // theta_dim is the dimension of theta (at the first-order, i.e. the order that directly affects item responses for an individual)
//           for(int d_ = 0; d_ < n_dimensions; d_++) { 
//             if(d_current_index <= dimension_end(d_)) {
//               theta_dim = d_ + 1;
//               break;
//             }
//           }
//           theta_prob = 1.0 / (1.0 + std::exp(-(x((theta_dim*d_nrow + i_)*n + n_) - x((t-1)*n + n_))));
//           lk += data(i_, d_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, d_current_index - 1))*std::log(1.0 - theta_prob + eps);
//         }
//       } else if(t <= n_theta + n_delta + n_alpha) {
//         for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (alpha step)
//           int a_current_index = t - n_theta - n_delta;
//           // loop through dimension indexes to discover which dimension the given question t belongs to
//           // store the dimension d_ for calling up the right theta to condition the alpha likelihood on it
//           int theta_dim = 0; // theta_dim is the dimension of theta (at the first-order, i.e. the order that directly affects item responses for an individual)
//           for(int d_ = 0; d_ < n_dimensions; d_++) {
//             if(a_current_index <= dimension_end(d_)) {
//               theta_dim = d_ + 1;
//               break;
//             }
//           }
//           theta_prob = 1.0 / (1.0 + std::exp(-(   (x((t-1)*n + n_)*x((theta_dim*d_nrow + i_)*n + n_)  )  - x((t-n_delta-1)*n + n_)     ))); //t-n_delta-1???
//           lk += data(i_, a_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, a_current_index - 1))*std::log(1.0 - theta_prob + eps);
//         }
//       } else {
//         // // loop through dimension indexes to discover which dimension the given question t belongs to
//         // // store the dimension d_ for calling up the right theta to condition the alpha likelihood on it
//         // int theta_dim = 0; // theta_dim is the dimension of theta (at the first-order, i.e. the order that directly affects item responses for an individual)
//         // int theta_id = t - n_delta - n_alpha; // which theta are we discussing?
//         // for(int d_ = 0; d_ < n_dimensions; d_++) {
//         //   if(t <= dimension_end(d_)) {
//         //     theta_dim = d_ + 1;
//         //     break;
//         //   }
//         // }
//         // for(int i_ = 0; i_ < d_nrow; i_++) { // likelihood part for relating the first order theta dimension to the second-order theta dimension
//         //   int theta_current_index = i_ + 1;

//         //   lk += log_normd_dx(x((theta_dim*d_nrow + i_)*n + n_), x((t-1)*n + n_)*x((theta_current_index-1)*n + n_), 1.0);
//         // }
//       }
//       u(n_) = lk;   // update using the likelihood at "prior distribution" draw (log scale)
//       w(n_) = w(n_) + u(n_);  // update weights on the log scale
//     }
//     // This section transforms the log weights by exponentiating them.
//     // We cannot just exponentiate the weights, so we take the 
//     // difference between weight and max weight,
//     // and then exponentiate the difference minus logadd.
//     // The weight adjustment method was taken from the SIR method of LaplacesDemon R package
//     // -----------------------------------------
//     double md = arma::max(w);
//     arma::vec lw = w - md;
//     arma::vec probs = arma::exp(lw - logadd(lw));
//     // -----------------------------------------
//     // End weight scale change
//     // Update weight history
//     for(int w_ = 0; w_ < n; w_++) {
//       w_history((t-1)*n + w_) = probs(w_);
//     }
//     double norm_const = arma::accu(probs);   // sum of probs should be 1 since probs are probabilities
//     // weight observations by their importance. This is the step by which the prior draws
//     // are weighted by the likelihood, yielding draws from the joint posterior distribution.
//     for(int n_ = 0; n_ < n; n_++) { 
//       x_weighted((t-1)*n + n_) = x((t-1)*n + n_) * probs(n_) / norm_const;
//     }
//     // This is the stage where weights are checked for degeneracy (and rejuvinated).
//     // 1) In short this section does the following: compute some measure of "issues" with the weights.
//     // 2) If "issues" is sufficiently small (indicated by a cutoff value), then we know that the weights are
//     // degenerate or approaching degeneracy.
//     // 3) Fix the degeneracy by resampling all parameter draws with probability proportional to weigths and
//     //    reset all weights to 1 / n. This will be triggered a number of times during runtime to repair
//     //    weights that are overly concentrated near zero, and keep going. Without this step, the algorithm 
//     //    will not work.
//     n_eff = arma::accu(exp(w) % exp(w));     // check weights for issues using sum of squared weights
//     if(n_eff < tol) {   // if weights are under the minimum that we decide is tolerable, then do 2) and 3).
//       // update each parameter draw with the resampled values
//       arma::ivec sindex = arma::regspace<arma::ivec>((t-1)*n, t*n - 1); //0  // vector from 0, ..., t*n - 1 (zero up to current)
//       arma::vec prev_probs = arma::vec(n, arma::fill::zeros);  //t*n  // get previous probabilities for the resampling step
//       for(int w_ = (t-1)*n; w_ < t*n; w_++) { //w_=0  // for every previous probability, use the weight which was used earlier 
//         prev_probs(w_ - (t-1)*n) = w_history(w_);
//       }
//       prev_probs = prev_probs / arma::accu(prev_probs);
//       // Now, we complete 3) by resampling using weight history
//       arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, n, true, prev_probs); // t*n
//       for(int x_t = (t-1)*n; x_t < t*n; x_t++) { //t*n
//         x(x_t) = x(samp_indices(x_t - (t-1)*n));
//       }
//       for(int n_ = 0; n_ < n; n_++) {   // weight = 1.0 / n for each weight. This is a reset of weights with equal value per
//         w(n_) = 1.0 / n;              // weight. It's like starting the weighting over again.
//       }
//     }   // End of the weight rejuvination step.
//     t += 1;  // That was for one parameter (t = 1), increment t by one and repeat for the rest of the parameters
//   }
//   return x_weighted;  // return vector of weighted estimates (will have to sum across draws for the expected value (posterior mean))
// }





























// [[Rcpp::export]]
arma::vec sis5(arma::mat & data, int n_dimensions, arma::vec dimension_start, arma::vec dimension_end, int n, double tol) {
  // Approximate draw from joint posterior distribution by sequential update with resampling
  // using the sequential importance sampling with resampling angorithm 
  // Reference: 
  // Arguments: 
  //  data - input matrix of observations
  //  n_dimensions - number of first-order dimensions
  //  dimension_start - start index (from 1) of each dimension. The dimension ends at dimension_start(n+1) - 1
  //  dimension_end - the end index of each dimension.
  //  n - number of resampling iterations
  //  tol - tolerance for resampling. If sum of squared weight falls below tol,
  //        resampling is performed and weights reset.
  // Value:
  //  vector of weighted estimates for each parameter over each of the n resampling iterations
  double eps = 1e-16;                       // epsilon for log() when performing likelihood calculation
  int d_nrow = data.n_rows;                 // number of rows in data
  int d_ncol = data.n_cols;                 // number of cols in data
  // n_dimensions                           // number of latent (ability) dimensions
  int n_theta = d_nrow*(n_dimensions + 1);  // number of theta parameters (theta = latent ability)
  int n_delta = d_ncol;                     // number of deltas, difficult parameters
  int n_alpha = d_ncol;                     // number of thetas, ability parameters
  int n_lambda = n_dimensions;              // number of lambda parameters (loadings for second-order theta on each of the first-order thetas)
  int n_param = n_lambda + n_theta + n_delta + n_alpha;  // total number of parameters to estimate
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
  double n_eff = 0.0;        // pre allocation of memory for tolerance test
  // loop over each parameter
  //   for each parameter, resample n times
  //   update likelihood (weights), save history, check for degeneracy
  //   if weights are near degenerate, resample parameters and reset weights
  for(int p_ = 0; p_ < n_param; p_++) {   // p_ indexes current parameter to update
    for(int n_ = 0; n_ < n; n_++) {         // n_ indexes the resampling iterations (sample each "p_" "n_" times)
      if(t <= n_lambda) {                   // if t <= n_lambda, sample lambda prior
        x((t-1)*n + n_) = truncated_normal_ab_sample(0.0, 2.0, -5.0, 5.0, seeds((t-1)*n + n_));
      } else if(t <= n_lambda + n_theta) {  // if t <= n_lambda + n_theta, sample theta prior
        if(t <= n_lambda + d_nrow) {        // if t is first dimension of theta (second-order theta, no prior on the prior's mean)
          x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0;
        } else {      // second, loop over each first order dimension of theta and set prior mean to second-order theta
          // calculate which dimension we are in
          int current_dimension = (t-n_lambda-1) / d_nrow;          // locate the current lambda parameter to load on proir 
          int thetag_index = ((t-n_lambda-1) % d_nrow) + n_lambda;  // locate the current second-order theta parameter to load on prior
          // compute prior draw for the given step (theta step, all first-order dimensions)
          // conditional on second-order theta
          x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*1.0 + (x((current_dimension)*n + n_) * x((thetag_index)*n + n_));
        }
      } else if(t <= n_lambda + n_theta + n_delta) {
        x((t-1)*n + n_) = normal_01_sample(seeds((t-1)*n + n_))*2.0;
      } else { //if(t <= n_lambda + n_theta + n_delta + n_alpha)
        //x((t-1)*n + n_) = truncated_normal_ab_sample(1.0, 2.0, 0.0, , seeds((t-1)*n + n_));
        double a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0;
        while (a_samp < 0.0) {
          a_samp = normal_01_sample(seeds((t-1)*n + n_))*2.0 + 1.0; // Probably switch this back to the truncated normal dist rather than inefficiency of over-sampling
        }
        x((t-1)*n + n_) = a_samp;
      }
      double lk = 0.0;  // pre-allocate likelihood memory
      // --------------------------
      if(t <= n_lambda) {
        for(int j_ = dimension_start(t-1) - 1; j_ < dimension_end(t-1) - 1; j_++) {
          for(int i_ = 0; i_ < d_nrow; i_++) {
            theta_prob = 1.0 / (1.0 + std::exp(-(x((t-1)*n + n_))));
            lk += data(i_, j_)*std::log(theta_prob + eps) + (1.0 - data(i_, j_))*std::log(1.0 - theta_prob + eps);
          }
        }
      } else if(t <= n_lambda + n_theta) {
        // first, estimate the second order theta
        if(t <= n_lambda + d_nrow) {
          for(int j_ = 0; j_ < d_ncol; j_++) {  // compute the likelihood for the given step (theta step, second-order dimension)
            int theta_current_index = t-n_lambda-1;
            theta_prob = 1.0 / (1.0 + std::exp(      -(    x((t-1)*n + n_)   )          ));
            lk += data(theta_current_index, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_current_index, j_))*std::log(1.0 - theta_prob + eps);
          }
        } else {    // second, loop over each first order dimension of theta
          // calculate which dimension we are in
          int current_dimension = (t-n_lambda-1) / d_nrow;
          int theta_row_index = (t-n_lambda-1) % d_nrow;
          for(int j_ = dimension_start(current_dimension - 1) - 1; j_ < dimension_end(current_dimension - 1) - 1; j_++) {  // compute the likelihood for the given step (theta step, all first-order dimensions)
            theta_prob = 1.0 / (1.0 + std::exp(      -(    x((t-1)*n + n_)    )             ));
            lk += data(theta_row_index, j_)*std::log(theta_prob + eps) + (1.0 - data(theta_row_index, j_))*std::log(1.0 - theta_prob + eps);
          }
        }
      } else if(t <= n_lambda + n_theta + n_delta) {
        for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (delta step)
          int d_current_index = t - n_lambda - n_theta;
          // loop through dimension indexes to discover which dimension the given question t belongs to
          // store the dimension d_ for calling up the right theta to condition the alpha likelihood on it
          int theta_dim = 0; // theta_dim is the dimension of theta (at the first-order, i.e. the order that directly affects item responses for an individual)
          for(int d_ = 0; d_ < n_dimensions; d_++) { 
            if(d_current_index <= dimension_end(d_)) {
              theta_dim = d_ + 1;
              break;
            }
          }
          theta_prob = 1.0 / (1.0 + std::exp(-(    x((n_lambda + (theta_dim*d_nrow) + i_)*n + n_)    - x((t-1)*n + n_)         ));
          lk += data(i_, d_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, d_current_index - 1))*std::log(1.0 - theta_prob + eps);
        }
      } else {  // if(t <= n_lambda + n_theta + n_delta + n_alpha)
        for(int i_ = 0; i_ < d_nrow; i_++) {  // compute the likelihood for the given step (alpha step)
          int a_current_index = t - n_lambda - n_theta - n_delta;
          // loop through dimension indexes to discover which dimension the given question t belongs to
          // store the dimension d_ for calling up the right theta to condition the alpha likelihood on it
          int theta_dim = 0; // theta_dim is the dimension of theta (at the first-order, i.e. the order that directly affects item responses for an individual)
          for(int d_ = 0; d_ < n_dimensions; d_++) {
            if(a_current_index <= dimension_end(d_)) {
              theta_dim = d_ + 1;
              break;
            }
          }
          theta_prob = 1.0 / (1.0 + std::exp(-(   (x((t-1)*n + n_) * x((n_lambda + (theta_dim*d_nrow) + i_)*n + n_)  )  - x((t-n_delta-1)*n + n_)     ))); //t-n_delta-1???
          lk += data(i_, a_current_index - 1)*std::log(theta_prob + eps) + (1.0 - data(i_, a_current_index - 1))*std::log(1.0 - theta_prob + eps);
        }
      }
      u(n_) = lk;   // update using the likelihood at "prior distribution" draw (log scale)
      w(n_) = w(n_) + u(n_);  // update weights on the log scale
    }
    // This section transforms the log weights by exponentiating them.
    // We cannot just exponentiate the weights, so we take the 
    // difference between weight and max weight,
    // and then exponentiate the difference minus logadd.
    // The weight adjustment method was taken from the SIR method of LaplacesDemon R package
    // -----------------------------------------
    double md = arma::max(w);
    arma::vec lw = w - md;
    arma::vec probs = arma::exp(lw - logadd(lw));
    // -----------------------------------------
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
    // This is the stage where weights are checked for degeneracy (and rejuvinated).
    // 1) In short this section does the following: compute some measure of "issues" with the weights.
    // 2) If "issues" is sufficiently small (indicated by a cutoff value), then we know that the weights are
    // degenerate or approaching degeneracy.
    // 3) Fix the degeneracy by resampling all parameter draws with probability proportional to weigths and
    //    reset all weights to 1 / n. This will be triggered a number of times during runtime to repair
    //    weights that are overly concentrated near zero, and keep going. Without this step, the algorithm 
    //    will not work.
    n_eff = arma::accu(exp(w) % exp(w));     // check weights for issues using sum of squared weights
    if(n_eff < tol) {   // if weights are under the minimum that we decide is tolerable, then do 2) and 3).
      // update each parameter draw with the resampled values
      arma::ivec sindex = arma::regspace<arma::ivec>((t-1)*n, t*n - 1); //0  // vector from 0, ..., t*n - 1 (zero up to current)
      arma::vec prev_probs = arma::vec(n, arma::fill::zeros);  //t*n  // get previous probabilities for the resampling step
      for(int w_ = (t-1)*n; w_ < t*n; w_++) { //w_=0  // for every previous probability, use the weight which was used earlier 
        prev_probs(w_ - (t-1)*n) = w_history(w_);
      }
      prev_probs = prev_probs / arma::accu(prev_probs);
      // Now, we complete 3) by resampling using weight history
      arma::ivec samp_indices = Rcpp::RcppArmadillo::sample(sindex, n, true, prev_probs); // t*n
      for(int x_t = (t-1)*n; x_t < t*n; x_t++) { //t*n
        x(x_t) = x(samp_indices(x_t - (t-1)*n));
      }
      for(int n_ = 0; n_ < n; n_++) {   // weight = 1.0 / n for each weight. This is a reset of weights with equal value per
        w(n_) = 1.0 / n;              // weight. It's like starting the weighting over again.
      }
    }   // End of the weight rejuvination step.
    t += 1;  // That was for one parameter (t = 1), increment t by one and repeat for the rest of the parameters
  }
  return x_weighted;  // return vector of weighted estimates (will have to sum across draws for the expected value (posterior mean))
}