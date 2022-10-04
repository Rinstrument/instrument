
#data[i, jj] = (1 / (1 + exp(-(alpha[jj]*theta[i] - delta[jj])))) > runif(1)

# mle2(llk_2pl_logit_np, start = list(x = start_vals), data = list(data = dat, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta = ))
# llk_2pl_logit_np(c(alpha, delta), data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta_start)
# llk_2pl_logit_np(c(alpha, delta), data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta)
# 
# mle2(llk_2pl_logit_np, start = list(x = start_vals),
#      data = list(data = dat, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta = theta_start))
# 
# mle2(llk_2pl_logit_np, start = list(x = start_vals),
#      data = list(data = dat, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta = theta_start))
# # apply(data, 2, mean)
# # apply(data, 1, mean)
# # coef(mirt(data, model = 1), simplify = TRUE)$items
# # fscores(mirt(data, model = 1))
# # cor(fscores(mirt(data, model = 1)), theta)
# # theta_start = apply(data, 1, mean)
# # f <- function(x, a) { 1 / (1 + exp(-((a[1] * x) - d[2]))) }
# # nlm(f, x = data[,1], a = c(3,5))
# # lm(data[,1] ~ theta_start)
# # glm(data[,1] ~ theta_start, family = binomial(link = "logit"))
# # dt = cbind(data[,4], theta)
# # colnames(dt) = c("y", "x")
# # m = nls(y ~ 1 / (1 + exp(-((a * x) - d))),
# #                  data = as.data.frame(dt),
# #                  start = list(a = 1, d = 0),
# #                  algorithm = "plinear")
# # summary(m)
# # get_start.lmHOIRT = function(data) {
# #   n = nrow(data)
# #   j = ncol(data)
# #   theta_start = rowMeans(data)
# #   alpha_start = rep(1, j)
# #   delta_start = rep(0, j)
# #   for(jj in 1:j) {
# #     jj = 19
# #     dt = cbind(data[, jj], rowMeans(data[, -jj]))
# #     #lm(dt[,1] ~ dt[,2])
# #     colnames(dt) = c("y", "x")
# #     dt[, 1] = ifelse(dt[, 1] < 0.5, dt[, 1] + runif(n, 0.01, 0.1), dt[, 1] - runif(n, 0.01, 0.1))
# #     #dt[,1] = abs(runif(dt[,1], -0.05, 0.05))
# #     #dt[,1] = ifelse(dt[,1] < 0.5, dt[,1] + runif(n, 0.01, 0.1), dt[,1] - runif(n, 0.01, 0.1))
# #     #lm(dt[,1], dt[,2])
# #     m = nls(y ~ 1 / (1 + exp(-((a * x) - d))), data = as.data.frame(dt),
# #             start = list(a = 2, d = -2), algorithm = "plinear")
# #     mcf = coef(m)
# #     mcf
# #     alpha_start[jj] = mcf["a"]
# #     delta_start[jj] = mcf["b"]
# #   }
# #   return(list(theta_start, alpha_start, beta_start))
# # }
#
# get_start.lmHOIRT(data)
# #start = c(rep(1, j), rep(0, j), theta_start)
# start = c(rep(1, j), rep(0, j), rep(0, n))

# 
#indices = indices[(burn + 1):iterations, ]
# get_draws = function(x, param, indices) {
#   indices = indices + 1
#   x[param, indices == param]
# }
#llk_2pl_logit(true, data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1), (2 * j):(p - 1)))
# optim(par = c(rep(1, j), rep(0, j)), fn = llk_2pl_logit_np, method = "L-BFGS-B",
#       data = data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta = rowMeans(data))
# llk_2pl_logit_np(true, data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), rowMeans(data))
# llk_2pl_logit_np(true, data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1)), theta)


// [[Rcpp::export]]
double llk_2pl_logit_np(arma::vec & x, arma::mat & data, std::list<arma::uvec> parameter_indexes,
                        arma::vec & theta) {
  // three sets of parameters:
    // parameter_indexes 0: alpha
  //                   1: beta
  //                   2: theta
  arma::uvec alpha_pos = get(parameter_indexes, 0);
  arma::uvec delta_pos = get(parameter_indexes, 1);
  //arma::uvec theta_pos = get(parameter_indexes, 2);
  double eps = 1e-16;
  int n = data.n_rows;
  int j = data.n_cols;
  
  arma::mat theta_mat = arma::diagmat(theta);
  arma::mat ones = arma::mat(n, j, arma::fill::ones);
  theta_mat = theta_mat * ones;
  arma::mat delta_mat = arma::diagmat(x(delta_pos));
  delta_mat = ones * delta_mat;
  arma::mat alpha_mat = arma::diagmat(x(alpha_pos));
  alpha_mat = ones * alpha_mat;
  
  theta_mat = 1.0 / (1.0 + arma::exp(-((alpha_mat % theta_mat) - delta_mat)));
  arma::mat data_refl = data;
  data_refl = -(data_refl - 1.0);
  theta_mat = arma::log((theta_mat % data) + data_refl - (theta_mat % data_refl) + eps);
  double llk = -arma::accu(theta_mat);
  // + arma::accu(log_normd(x(theta_pos), 0.0, 1.0)) + 
    // arma::accu(log_normd(x(delta_pos), 0.0, 2.0)) + arma::accu(arma::log(lognorm_dens(x(alpha_pos), 0.0, 1.0) + eps));
  return llk;
}

// [[Rcpp::export]]
double llk_2pl_logit(arma::vec & x, arma::mat & data, std::list<arma::uvec> parameter_indexes) {
  // three sets of parameters:
    // parameter_indexes 0: alpha
  //                   1: beta
  //                   2: theta
  arma::uvec alpha_pos = get(parameter_indexes, 0);
  arma::uvec delta_pos = get(parameter_indexes, 1);
  arma::uvec theta_pos = get(parameter_indexes, 2);
  double eps = 1e-16;
  int n = data.n_rows;
  int j = data.n_cols;
  
  arma::mat theta_mat = arma::diagmat(x(theta_pos));
  arma::mat ones = arma::mat(n, j, arma::fill::ones);
  theta_mat = theta_mat * ones;
  arma::mat delta_mat = arma::diagmat(x(delta_pos));
  delta_mat = ones * delta_mat;
  arma::mat alpha_mat = arma::diagmat(x(alpha_pos));
  alpha_mat = ones * alpha_mat;
  
  theta_mat = 1.0 / (1.0 + arma::exp(-((alpha_mat % theta_mat) - delta_mat)));
  arma::mat data_refl = data;
  data_refl = -(data_refl - 1.0);
  theta_mat = arma::log((theta_mat % data) + data_refl - (theta_mat % data_refl) + eps);
  double llk = arma::accu(theta_mat) + arma::accu(log_normd(x(theta_pos), 0.0, 1.0)) + 
    arma::accu(log_normd(x(delta_pos), 0.0, 2.0)) + arma::accu(arma::log(lognorm_dens(x(alpha_pos), 0.0, 1.0) + eps));
  return llk;
}



# benchmarking
# library(rbenchmark)
# 
benchmark("first_version" = {
  for(ii in 1:10) {
    llk_2pl_logit(true, data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1), (2 * j):(p - 1)))
  }
},
"second_version" = {
  for(ii in 1:10) {
    llk_2pl_logit_v2(true, data, parameter_indexes = list(0:(j - 1), j:(2 * j - 1), (2 * j):(p - 1)))
  }
},
replications = 1000,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))



