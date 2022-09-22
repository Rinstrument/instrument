#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::export]]
double min(double a, double b) {
  return 0.5 * (a + b - std::abs(a - b));
}

// [[Rcpp::export]]
arma::uvec get(std::list<arma::uvec> list, int l){
  std::list<arma::uvec>::iterator it = list.begin();
  for(int i = 0; i < l; i++){
    it++;
  }
  return *it;
}

// [[Rcpp::export]]
double llk_2pl_logit(arma::vec & x, arma::mat & data, std::list<arma::uvec> parameter_indexes) {
  // three sets of parameters:
  // parameter_indexes 1: alpha
  //                   2: beta
  //                   3: theta
  arma::uvec alpha_pos = get(parameter_indexes, 0);
  arma::uvec delta_pos = get(parameter_indexes, 1);
  arma::uvec theta_pos = get(parameter_indexes, 2);
  
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
  
  double llk = arma::accu(arma::log(theta_mat + 0.000000000000001)) + 
    arma::accu(arma::log(arma::normpdf(x(theta_pos), 0.0, 1.0) + 0.000000000000001)) + 
    arma::accu(arma::log(arma::normpdf(x(delta_pos), 0.0, 2.0) + 0.000000000000001)) + 
    arma::accu(arma::log(arma::log_normpdf(x(alpha_pos), 0.0, 2.0) + 200.0) + 0.000000000000001);
  
  return llk;
}

// [[Rcpp::export]]
SEXP arwmh(arma::mat & x, int p, arma::vec x_current, int iter, const double & a, arma::mat & data,
    int log_likelihood_selector, std::list<arma::uvec> parameter_indexes, arma::vec & accept) {

  double l_def = 2.38 * 2.38 / 1;

  arma::vec prop_mu_t = arma::vec(p, arma::fill::zeros);
  arma::vec prop_sigma_t = arma::vec(p, arma::fill::ones);

  arma::vec prop_mu_t1 = arma::vec(p, arma::fill::zeros);
  arma::vec prop_sigma_t1 = arma::vec(p, arma::fill::ones);

  arma::vec l_scaling_t = arma::vec(p, arma::fill::value(l_def));
  l_scaling_t = arma::log(l_scaling_t);
  arma::vec l_scaling_t1 = arma::vec(p, arma::fill::value(l_def));
  arma::vec gam_correct = arma::vec(p, arma::fill::ones);
  arma::vec gam_correct_iter = arma::vec(p, arma::fill::zeros);
  
  // random normal draws for proposal distribution
  arma::vec rn_prop_scale01 = arma::randn(iter, arma::distr_param(0, 1));
  arma::vec ru_prop_mh = arma::randu(iter, arma::distr_param(0, 1));

  arma::ivec p_update_index = arma::randi(iter, arma::distr_param(-0, +(p - 1)));
  arma::vec x_proposal = x_current;

  double R = 1;
  double r_accept = 0;
  double x_current_diff_prop_mu = 0;
  
  double (* log_likelihood)(arma::vec &, arma::mat &, std::list<arma::uvec>);
  
  if(log_likelihood_selector == 0) { // 1-P logit
    log_likelihood = llk_2pl_logit;
  } else if(log_likelihood_selector == 1) { // 2-P logit
    log_likelihood = llk_2pl_logit;
  }

  for(int it = 0; it < iter; it++) {

    int p_current_select = p_update_index(it);
    //gam_correct(p_current_select) = 1.0 / (it + 2.0);
    gam_correct(p_current_select) = 1.0 / (gam_correct_iter(p_current_select) + 2.0);
    gam_correct_iter(p_current_select) += 1.0;
    
    // // 1. Sample candidate value for component p_i
    x_proposal(p_current_select) = rn_prop_scale01(it) * (std::exp(l_scaling_t(p_current_select)) * prop_sigma_t(p_current_select)) + x_current(p_current_select);
    R = min(std::exp(log_likelihood(x_proposal, data, parameter_indexes) - log_likelihood(x_current, data, parameter_indexes)), 1);
    r_accept = R > ru_prop_mh(it);
    
    // Rcpp::Rcout << "    " << std::endl;
    // 
    // Rcpp::Rcout << "(p_current_select)" << p_current_select << std::endl;
    // Rcpp::Rcout << "x_current(p_current_select)" << x_current(p_current_select) << std::endl;
    // Rcpp::Rcout << "x_proposal(p_current_select)" << x_proposal(p_current_select) << std::endl;
    // Rcpp::Rcout << "log_likelihood(x_proposal, data, parameter_indexes)" << log_likelihood(x_proposal, data, parameter_indexes) << std::endl;
    // Rcpp::Rcout << "log_likelihood(x_current, data, parameter_indexes)" << log_likelihood(x_current, data, parameter_indexes) << std::endl;
    // 
    // Rcpp::Rcout << "    " << std::endl;
    
    // 2. Select the value for x(t+1) according to MH criterion
    if(r_accept == 1) {
      x_current(p_current_select) = x_proposal(p_current_select);
      accept(it) = 1.0;
      Rcpp::Rcout << "accepted" << std::endl;
      Rcpp::Rcout << "total: " << arma::accu(accept) << std::endl;
    }

    // 3. Adaptation step: update proposal distribution variance in two
    // steps and update component scaling parameter
    x_current_diff_prop_mu = x_current(p_current_select) - prop_mu_t(p_current_select);
    
    x.col(it) = x_current;
    
    prop_mu_t1(p_current_select) = prop_mu_t(p_current_select) + (gam_correct(p_current_select) * x_current_diff_prop_mu);
    prop_sigma_t1(p_current_select) = prop_sigma_t(p_current_select) + gam_correct(p_current_select) * ((x_current_diff_prop_mu * x_current_diff_prop_mu) - prop_sigma_t(p_current_select));
    l_scaling_t1(p_current_select) = l_scaling_t(p_current_select) + gam_correct(p_current_select) * (R - a);
    
    prop_mu_t(p_current_select) = prop_mu_t1(p_current_select);
    prop_sigma_t(p_current_select) = prop_sigma_t1(p_current_select);
    l_scaling_t(p_current_select) = l_scaling_t1(p_current_select);
    
    // Rcpp::Rcout << "    " << std::endl;
    // 
    // Rcpp::Rcout << "prop_mu_t(p_current_select)" << prop_mu_t(p_current_select) << std::endl;
    // Rcpp::Rcout << "prop_sigma_t(p_current_select)" << prop_sigma_t(p_current_select) << std::endl;
    // Rcpp::Rcout << "l_scaling_t(p_current_select)" << l_scaling_t(p_current_select) << std::endl;
    // 
    // Rcpp::Rcout << "    " << std::endl;
    
  }

  arma::vec lambda_p = arma::vec(p, arma::fill::value(l_def));

  return R_NilValue;
}

/*** R
set.seed("7794153")
n = 300
j = 20
p = 2 * j + n
iterations = 200000

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
accept = rep(0, iterations)

#start = c(c(0.7, alpha[-1]), delta, theta)
#start = c(alpha, delta, theta)
#llk_2pl_logit(x = start, data = data, list(0:(j - 1), j:(2 * j - 1), (2 * j):(p - 1)))
start = c(exp(runif(j)), rnorm(j, 0, 2), rnorm(n))
14.4264
true2 = true
true2[16] = 14
llk_2pl_logit(x = true2, data = data, list(0:(j - 1), j:(2 * j - 1), (2 * j):(p - 1)))

arwmh(x, p, start, iterations, a = 0.44, data, log_likelihood_selector = 0,
      parameter_indexes = list(0:(j - 1), j:(2 * j - 1), (2 * j):(p - 1)),
      accept)

*/
