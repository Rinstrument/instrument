# INIRT higher-order IRT tests
devtools::install(dependencies = FALSE)
library(rstan)
stanc(file = "./inst/stan/inirt_hoirt2.stan", verbose = TRUE)
n = 1000
d = 4
j = 10*d
ncat = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
k = 0
uk = 0
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n)))
b_alpha = 1
alpha_dstart = c(1, 11, 21, 31)
alpha_dend = c(10, 20, 30, 40)
for(dd in 1:d) {
  alpha[dd, alpha_dstart[dd]:alpha_dend[dd]] = sort(runif(j/d, -1.5, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 1
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta_g = rnorm(n, 0, 2)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, sqrt(1.5)) #0.1??????
}
beta = NULL
predictors = NULL
start_index = 1
beta_dstart = NULL
beta_dend = NULL
# for(dd in 1:1) {
#   beta_dstart[dd] = start_index
#   beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
#   start_index = start_index + length(predictors[[dd]])
# }
# beta_mat = matrix(0, nrow = k, ncol = 1)
# index = 1
# for(dd in 1:1) {
#   for(i in beta_dstart[dd]:beta_dend[dd]) {
#     beta_mat[i, dd] = beta[index]
#     index = index + 1
#   }
# }
# x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
lambda_ind = rep(1:d, each = j/d)
lambda = c(0.3, 0.4, 0.5, -0.2)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) { #                                                                      + x[i, ] %*% beta_mat
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(lambda[lambda_ind[jj]]*(theta_g[i]) + theta[i, lambda_ind[jj]])) - (delta[jj, 1:ncategi[jj]] + b_delta*d_design[i,])))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
  }
}
apply(data, 2, table)
# remove_gaps = function(x) {
#   ord = order(x); vec = sort(x)
#   old = unique(vec); replace = 1:length(unique(vec))
#   names(replace) = old; names(vec) = vec
#   new = replace[names(vec)]; names(new) = NULL
#   return(new[ord])
# }
# data = apply(data, 2, remove_gaps)
# apply(data, 2, table)
# data = cbind(data, x)
colnames(data) = c(paste0("x", 1:j)) #, paste0("z", 1:k)
# for(dd in 1:1) {
#   predictors[[dd]] = predictors[[dd]] + j
# }
# dims = 3
# h2_dims = 1
# h2_dim_id = list(1:20, 21:40, 41:60)
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta, lambda = lambda,
    theta_g = theta_g)
fit_data = list(data = data)
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()
# fit = inirt::inirt(data = fit_data$data, model = fit_data$model, predictors = fit_data$predictors, dims = fit_data$dims, 
#     h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, 
#     init = fit_data$init)

library(devtools)
library(Rcpp)
compileAttributes()
load_all()
data = fit_data$data
model = "thetag ~ theta1 + theta2 + theta3 + theta4
         theta1 = c(1:10)
         theta2 = c(11:20)
         theta3 = c(21:30)
         theta4 = c(31:40)"
itype = "2pl"
method = "vb"
iter = 10000
tol_rel_obj = 1e-4
exploratory = FALSE
weights = NULL

data = fit_data$data
colnames(data)
fit = theta2::theta2(
  data = data,
  model = "thetag ~ theta1 + theta2 + theta3 + theta4
           theta1 = c(1:10)
           theta2 = c(11:20)
           theta3 = c(21:30)
           theta4 = c(31:40)",
  itype = "2pl",
  method = "vb",
  iter = 10000,
  tol_rel_obj = 1e-4)






# fit = inirt::inirt(data = fit_data$data, model = fit_data$model, predictors = fit_data$predictors, dims = fit_data$dims, 
#     h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, iter = fit_data$iter, 
#     init = fit_data$init)


# summary(fit, pars = "alpha")$summary[,"mean"]
aest = matrix(summary(fit, pars = "alpha")$summary[,"mean"], nrow = dims, byrow = TRUE)
aest
cor(aest[1,1:20], alpha[1,1:20])
cor(aest[2,21:40], alpha[2,21:40])
cor(aest[3,41:60], alpha[3,41:60])

dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,"mean"], nrow = 60, byrow = TRUE)
cor(dest[,1], delta[,2])
cor(dest[,2], delta[,3])
cor(dest[,3], delta[,4])

tgest = summary(fit, pars = c("theta_g"))$summary[,"mean"]
cor(tgest, theta_g)
plot(tgest, theta_g)
test = matrix(summary(fit, pars = c("theta_resid"))$summary[,"mean"], nrow = n, byrow = TRUE)
cor(tgest + test[,1], theta[,1])
cor(test[,1], theta[,1])
cor(test[,2], theta[,2])
cor(test[,3], theta[,3])

cor(tgest + test[,2], theta[,2])
cor(tgest + test[,3], theta[,3])

summary(fit, pars = c("lambda"))$summary[,"mean"]
