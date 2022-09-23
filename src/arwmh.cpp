
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::export]]
double min(double a, double b) {
  return 0.5 * (a + b - std::abs(a - b));
}

// [[Rcpp::export]]
double log_likelihood(arma::vec & x, arma::mat & data, arma::uvec & alpha_pos,
                      arma::uvec & delta_pos, arma::uvec & theta_pos) {
  
  int n = data.n_rows;
  int j = data.n_cols;
  
  arma::mat theta_mat = arma::diagmat(x(theta_pos));
  arma::mat ones = arma::mat(n, j, arma::fill::ones);
  
  theta_mat = theta_mat * ones;
  
  arma::mat delta_mat = arma::diagmat(x(delta_pos));
  delta_mat = ones * delta_mat;
  
  arma::mat alpha_mat = arma::diagmat(x(alpha_pos));
  alpha_mat = ones * alpha_mat;
  
  theta_mat = 1 / (1 + arma::exp(-alpha_mat % (theta_mat - delta_mat)));
  
  arma::mat data_refl = data;
  
  data_refl = -(data_refl - 1);
  
  theta_mat = theta_mat % data + data_refl - (theta_mat % data_refl); // this may be wrong?? p^theta * (1-p)^(1-theta)???????????? Does this do that?
  // did I do all of this on the log scale??????????????????????????????//
  
  //arma::mat temp = arma::mat(n, j, arma::fill::randu);
  
  double llk = arma::accu(arma::log(theta_mat + 0.0001)) + 
    arma::accu(arma::log(arma::normpdf(x(theta_pos), 0.0, 1.0) + 0.0001)) + 
    arma::accu(arma::log(arma::normpdf(x(delta_pos), 0.0, 2.0) + 0.0001)) + 
    arma::accu(arma::log(arma::log_normpdf(x(alpha_pos), 0.0, 2.0) + 20.0)) + 10000.0; // + 0.0000000001
  
  // Rcpp::Rcout << x(alpha_pos) << std::endl;
  // Rcpp::Rcout << arma::log_normpdf(x(alpha_pos), 0.0, 2.0) << std::endl;
  // Rcpp::Rcout << arma::log(arma::log_normpdf(x(alpha_pos), 0.0, 2.0) + 20.0000000001) << std::endl;
  // 
  //Rcpp::Rcout << llk << std::endl;
  
  return llk;
}

// [[Rcpp::export]]
SEXP arwmh(arma::mat & x, int p, arma::vec x_current, int iter,
           const double & a, arma::mat & data, arma::uvec & alpha_pos,
           arma::uvec & delta_pos, arma::uvec & theta_pos) {
  
  // case switch the model to fit
  
  double l_def = 2.38 * 2.38 / 100;
  
  arma::vec prop_mu_t = arma::vec(p, arma::fill::zeros);
  arma::vec prop_sigma_t = arma::vec(p, arma::fill::ones);
  
  arma::vec prop_mu_t1 = arma::vec(p, arma::fill::zeros);
  arma::vec prop_sigma_t1 = arma::vec(p, arma::fill::ones);
  
  arma::vec l_scaling_t = arma::vec(p, arma::fill::value(l_def));
  l_scaling_t = arma::log(l_scaling_t);
  arma::vec l_scaling_t1 = arma::vec(p, arma::fill::value(l_def));
  
  //arma::vec gam_correct = arma::linspace<arma::vec>(1, iter, iter + 1);
  
  arma::vec gam_correct = arma::vec(p, arma::fill::ones);
  //gam_correct -= 1;
  //gam_correct = 1 / ((gam_correct/10) + 1);
  
  // random normal draws for proposal distribution
  arma::vec rn_prop_scale01 = arma::randn(iter, arma::distr_param(0, 1));
  arma::vec ru_prop_mh = arma::randu(iter, arma::distr_param(0, 1));
  
  arma::ivec p_update_index = arma::randi(iter, arma::distr_param(-0, +(p - 1)));
  arma::vec x_proposal = x_current;
  
  double R = 1;
  double r_accept = 0;
  double x_current_diff_prop_mu = 0;
  
  for(int it = 0; it < iter; it++) {
    
    //int p_current_select = p_update_index(it);
    int p_current_select = 0;
    
    gam_correct(p_current_select) = 1.0 / (it + 2.0); // / (gam_correct(p_current_select) + 1);
    Rcpp::Rcout << "first gam_correct(p_current_select): " << gam_correct(p_current_select) << std::endl;
    
    // // 1. Sample candidate value for component p_i
    Rcpp::Rcout << "here: " << std::endl;
    Rcpp::Rcout << "rn_prop_scale01(it): " << rn_prop_scale01(it) << std::endl;
    Rcpp::Rcout << "l_scaling_t(p_current_select): " << l_scaling_t(p_current_select) << std::endl;
    Rcpp::Rcout << "rn_prop_scale01(it): " << rn_prop_scale01(it) << std::endl;
    Rcpp::Rcout << "x_current(p_current_select): " << x_current(p_current_select) << std::endl;
    
    x_proposal(p_current_select) = rn_prop_scale01(it) * (std::exp(l_scaling_t(p_current_select)) * prop_sigma_t(p_current_select)) + x_current(p_current_select); 
    
    Rcpp::Rcout << "Accepted: x_proposal(p_current_select): " << x_proposal(p_current_select) << std::endl;
    Rcpp::Rcout << "Accepted: x_current(p_current_select): " << x_current(p_current_select) << std::endl;
    
    R = min(
        std::exp(log_likelihood(x_proposal, data, alpha_pos, delta_pos, theta_pos) - 
          log_likelihood(x_current, data, alpha_pos, delta_pos, theta_pos)), 
        1); // log likelihood should have
    // the data and functions, indexes already contained, so the only argument is the parameter vector

    r_accept = R > ru_prop_mh(it);

    // 2. Select the value for x(t+1) according to MH criterion
    if(r_accept == 1) {
      x_current(p_current_select) = x_proposal(p_current_select);
      //Rcpp::Rcout << "accepted" << std::endl;
      
      
    }
    // } else {
    //   //x_proposal(p_current_select) = x_current(p_current_select);
    //   
    //   // // 3. Adaptation step: update proposal distribution variance in two
    //   // // steps and update component scaling parameter
    //   // x_current_diff_prop_mu = x_current(p_current_select) - prop_mu_t(p_current_select);
    //   // Rcpp::Rcout << "Rejected: x_proposal(p_current_select): " << x_proposal(p_current_select) << std::endl;
    //   // Rcpp::Rcout << "Rejected: prop_mu(p_current_select): " << prop_mu_t(p_current_select) << std::endl;
    //   // Rcpp::Rcout << "Rejected: x_current_diff_prop_mu: " << x_current_diff_prop_mu << std::endl;
    //   
    // }
    
    // 3. Adaptation step: update proposal distribution variance in two
    // steps and update component scaling parameter
    x_current_diff_prop_mu = x_current(p_current_select) - prop_mu_t(p_current_select);
    Rcpp::Rcout << "Accepted: x_current(p_current_select): " << x_current(p_current_select) << std::endl;
    Rcpp::Rcout << "Accepted: prop_mu(p_current_select): " << prop_mu_t(p_current_select) << std::endl;
    Rcpp::Rcout << "Accepted: x_current_diff_prop_mu: " << x_current_diff_prop_mu << std::endl;

    x.col(it) = x_current;

    
    
    //Rcpp::Rcout << "prop_mu(p_current_select) before: " << prop_mu(p_current_select) << std::endl;
    //Rcpp::Rcout << "gam_correct(p_current_select): " << gam_correct(p_current_select) << std::endl;
    //Rcpp::Rcout << "x_current_diff_prop_mu: " << x_current_diff_prop_mu << std::endl;
    prop_mu_t1(p_current_select) = prop_mu_t(p_current_select) + (gam_correct(p_current_select) * x_current_diff_prop_mu);
    //Rcpp::Rcout << "prop_mu(p_current_select): " << prop_mu(p_current_select) << std::endl;
    
    Rcpp::Rcout << "prop_sigma_t(p_current_select): " << prop_sigma_t(p_current_select) << std::endl;
    Rcpp::Rcout << "gam_correct(p_current_select): " << gam_correct(p_current_select) << std::endl;
    Rcpp::Rcout << "x_current_diff_prop_m: " << x_current_diff_prop_mu << std::endl;
    prop_sigma_t1(p_current_select) = prop_sigma_t(p_current_select) + (gam_correct(p_current_select) * ((x_current_diff_prop_mu * x_current_diff_prop_mu)) - prop_sigma_t(p_current_select));
    Rcpp::Rcout << "prop_sigma_t1(p_current_select) after: " << prop_sigma_t1(p_current_select) << std::endl;
    
    l_scaling_t1(p_current_select) = l_scaling_t(p_current_select) + gam_correct(p_current_select) * (R - a);
    
    Rcpp::Rcout << "l_scaling_t1(p_current_select): " << l_scaling_t1(p_current_select) << std::endl;
    
    
    Rcpp::Rcout << "    " << std::endl;
    
    prop_mu_t(p_current_select) = prop_mu_t1(p_current_select);
    prop_sigma_t(p_current_select) = prop_sigma_t1(p_current_select);
    l_scaling_t(p_current_select) = l_scaling_t1(p_current_select);
    
    //x_proposal(p_current_select) = x_current(p_current_select);
  }
  
  arma::vec lambda_p = arma::vec(p, arma::fill::value(l_def));

  return R_NilValue;
}

/*** R
set.seed("7794153")
n = 300
j = 20
p = 2 * j + n
iterations = 10

alpha = exp(runif(j))
delta = rnorm(j, 0, 2)
theta = rnorm(n)

true = c(alpha, delta, theta)

data = matrix(0, nrow = n, ncol = j)

for(i in 1:n) {
  for(jj in 1:j) {
    data[i, jj] = (1 / (1 + exp(-alpha[jj]*(theta[i] - delta[jj])))) > runif(1)
  }
}

x = matrix(data = 0, p, iterations)

start = c(c(0.7, alpha[-1]), delta, theta)
#start = c(alpha + exp(runif(j, 0, 0.1)), delta + rnorm(j, 0, 0.2), theta + rnorm(n, 0, 0.1))

#c(alpha, delta, theta)[548]

arwmh(x, p = p, x_current = start, iter = iterations, a = 0.44, data = data,
      alpha_pos = 0:(j - 1), delta_pos = j:(2 * j - 1), theta_pos = (2 * j):(p - 1))

# 0:(j - 1) + 1
#x2 = x[, 80000:100000]

#apply(x2, 1, mean)
#post.mean = apply(x2, 1, mean)
# post.mean[is.infinite(post.mean)] = NA
#post.compare = cbind.data.frame(true, post.mean)
# post.compare = post.compare[complete.cases(post.compare),]
#cor(post.compare)
#cor(cbind.data.frame(true, start))

# log_likelihood(x = start, data = data, alpha_pos = 0:(j - 1),
#       delta_pos = j:(2 * j - 1), theta_pos = (2 * j):(p - 1))
# log_likelihood(x = c(alpha, delta, theta), data = data, alpha_pos = 0:(j - 1), 
#                delta_pos = j:(2 * j - 1), theta_pos = (2 * j):(p - 1))

# -1.731417
# start_make_worse = c(alpha, delta, theta)
# start_make_worse[548+1] = -1.2
# start_make_less_worse = c(alpha, delta, theta)
# start_make_less_worse[548+1] = -1.74
# true = c(alpha, delta, theta)
# 
# log_likelihood(x = true, data = data, alpha_pos = 0:(j - 1),
#                delta_pos = j:(2 * j - 1), theta_pos = (2 * j):(p - 1))
# log_likelihood(x = start_make_worse, data = data, alpha_pos = 0:(j - 1), 
#                delta_pos = j:(2 * j - 1), theta_pos = (2 * j):(p - 1))
# log_likelihood(x = start_make_less_worse, data = data, alpha_pos = 0:(j - 1), 
#                delta_pos = j:(2 * j - 1), theta_pos = (2 * j):(p - 1))

*/
