# INIRT with univariate theta dimension
# test of the new alpha regression setup
devtools::install(dependencies = FALSE)
library(rstan)
stanc(file = "./inst/stan/inirt_unidim.stan", verbose = TRUE)
# library(devtools)
# Rcpp::compileAttributes()
# load_all()
# Test 1: simplest possible settings
n = 1000
ncat = 2
j = 35
d = 1
k = 0
uk = 0
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = rnorm(n, 0.1)))
b_alpha = c(0.4, 0)
for(dd in 1:d) {
  alpha[dd, ] = (sort(runif(j, -1.5, 1.5)))
}
alpha[1, ] = alpha[1, ] - mean(alpha[1, ])
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 0.7
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
eta = runif(j, 0.02, 0.4)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, sqrt(1.5))
}
beta = NULL
predictors = NULL
start_index = 1
beta_dstart = NULL
beta_dend = NULL
# Lz = 20
# z = matrix(rep(diag(20), each = 20), nrow = Lz*20) # for n = 800
# z = cbind(z, z[gtools::permute(1:nrow(z)), ])
# zeta_sd = 0.5
# zeta = rnorm(Lz, 0, sd = zeta_sd)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                                                   z[i, ] %*% zeta
    prb = eta[jj] + ((1 - eta[jj]) * (1 / (1 + exp(-(sum(exp(alpha[, jj] + as.vector(b_alpha%*%a_design[i,]))*(theta[i, ])) - (delta[jj, ] + as.vector(b_delta%*%d_design[i,])))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
apply(data, 2, table)
data = cbind(data, a_design)
data = as.data.frame(data)
colnames(data) = c(paste0("x", 1:j), paste0("ap", 1:2))
# data = cbind(data, z_fac = rep(paste0("z", 1:(Lz)), each = 20))
# ranef_id = c(rep(1, ncol(z)/2)) #, rep(2, ncol(z)/2)
dims = 1
item_id = 1:j
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
fit_data = list(data = data, item_id = item_id, model = NULL, predictors = predictors,
    n_pranef_cor = 0,
    dims = dims, h2_dims = 0, h2_dim_id = NULL, structural_design = list(alpha = a_design, delta = d_design), 
    method = "vb", weights = NULL, tol_rel_obj = 0.0002, iter = 5e3, init = "random")
# data = data; item_id = item_id; model = NULL; predictors = predictors; predictors_ranef = NULL; ranef_id = NULL; 
# predictors_ranef_corr = NULL; n_pranef_cor = NULL;
# dims = dims; h2_dims = 0; h2_dim_id = NULL; structural_design = list(alpha = a_design, delta = d_design); 
# method = "vb"; weights = NULL; tol_rel_obj = 0.0002; iter = 5e3; init = "random";
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

library(devtools)
library(Rcpp)
compileAttributes()
load_all()
data = fit_data$data
model = "theta = c(1:35)"
itype = "3pl"
method = "vb"
iter = 15000
tol_rel_obj = 1e-4
exploratory = FALSE
method = "vb"
weights = NULL

data = fit_data$data
colnames(data)
fit = theta2::theta2(
  data = data,
  model = "theta = c(1:35)",
  itype = "3pl",
  method = "vb", 
  iter = 10000, 
  tol_rel_obj = 1e-4)

fit = theta2::theta2(
  data = data,
  model = "theta = c(1:35)",
  itype = "3pl",
  method = "hmc", iter = 30, warmup = 15,
  chains = 1
  )


cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 35, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])
plot(dest[,1], sim_data$delta[,2])


rstan::summary(fit, pars = "alpha_l")$summary
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha_l")$summary[,1])
plot(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha_l")$summary[,1])

rstan::summary(fit, pars = "alpha_r_l")$summary
exp(rstan::summary(fit, pars = "alpha_r_l")$summary)


rstan::summary(fit, pars = "delta_r_l")$summary
fit@model_pars

cor(sim_data$zeta, (rstan::summary(fit, pars = "aeta_l")$summary[,1]))
plot(sim_data$zeta, (rstan::summary(fit, pars = "aeta_l")$summary[,1]))

rstan::summary(fit, pars = "aeta_l_sd")$summary





















# INIRT with univariate theta dimension
# test of the new alpha regression setup (correlated random effect)
devtools::install(dependencies = FALSE)
library(rstan)
stanc(file = "./inst/stan/inirt_unidim.stan", verbose = TRUE)
# ?stanc
# library(devtools)
# Rcpp::compileAttributes()
# load_all()
# Test 1: simplest possible settings
n = 800
ncat = 3
j = 35
d = 1
k = 0
uk = 0
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = rnorm(n, 0.1), x3 = rnorm(n, 0.1), 
                                x4 = rnorm(n, 0.1), x5 = sample(c(1,0), n, replace = TRUE), x6 = sample(c(1,0), n, replace = TRUE)))
b_alpha = c(0.8, 1.4, 0.7, 0.2, 0.4, -0.8)
exp(c(0.2, -1.4, 0.7, 0.2, 0.4, -0.8))
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, -1, 1.5))
  alpha[dd, ] = alpha[dd, ] - mean(alpha[dd, ])
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 1.3
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 2.5))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = NULL
predictors = NULL
start_index = 1
beta_dstart = NULL
beta_dend = NULL
Lz = 40
sigma = matrix(c(4, 3, 3, 3.5), nrow = 2)
int_slope = mvtnorm::rmvnorm(Lz, mean = c(0, 0), sigma = sigma) / 5
cov(int_slope)
z_c_int = matrix(rep(diag(40), each = 20), nrow = Lz*20)
zc_is = z_c_int %*% int_slope[,1]
z_c_slope = matrix(rep(diag(40), each = 20), nrow = Lz*20)
z_c_slope[z_c_slope == 1] = a_design[,2] #runif(20*20)
zc_is = zc_is + z_c_slope %*% int_slope[,2]
z_c = cbind(z_c_int, z_c_slope)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                               + zc_is[i,]                     z[i, ] %*% zeta
    prb = (1 / (1 + exp(-(sum(exp(alpha[, jj] + as.vector(b_alpha%*%a_design[i,]))*(theta[i, ])) - (delta[jj, ] + as.vector(b_delta%*%d_design[i,]))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
apply(data, 2, table)
data = cbind(data, a_design)
data = as.data.frame(data)
colnames(data) = c(paste0("x", 1:j), paste0("ap", 1:6)) #, paste0("z", 1:(2*Lz))
data = cbind(data, z_fac = rep(paste0("z", 1:(Lz)), each = 20))
# ranef_id = c(rep(1, ncol(z)/2)) #, rep(2, ncol(z)/2)
dims = 1
# item_id = 1:j           
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta, int_slope = int_slope)
fit_data = list(data = data, model = NULL, predictors = predictors,
    n_pranef_cor = 0,
    dims = dims, h2_dims = 0, h2_dim_id = NULL, structural_design = list(alpha = a_design, delta = d_design), 
    method = "vb", weights = NULL, tol_rel_obj = 0.0002, iter = 5e3, init = "random")
# data = data; item_id = item_id; model = NULL; predictors = predictors; predictors_ranef = NULL; ranef_id = NULL; 
# predictors_ranef_corr = NULL; n_pranef_cor = NULL;
# dims = dims; h2_dims = 0; h2_dim_id = NULL; structural_design = list(alpha = a_design, delta = d_design); 
# method = "vb"; weights = NULL; tol_rel_obj = 0.0002; iter = 5e3; init = "random";
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

# data = fit_data$data
# fit = theta2::theta2(
#   data = data,
#   model = "theta = c(1:25)
#            theta ~ 0
#            alpha ~ 1
#            delta ~ 1",
#   method = "hmc", iter = 300, warmup = 150,
#   chains = 1
#   )


data = fit_data$data
colnames(data)

library(devtools)
load_all()
weights = NULL
exploratory = FALSE
data = data
model = "theta = c(1:35)"
#   pre_start = FALSE,
itype = "1pl"
method = "vb"; iter = 15000; tol_rel_obj = 5e-4



# fit = theta2::theta2(
#   data = data,
#   model = "theta = c(1:35)",
# #   pre_start = FALSE,
#   itype = "2pl",
#   method = "vb", iter = 1500, tol_rel_obj = 5e-4
#   )

fit = theta2::theta2(
  data = data,
  model = "theta = c(1:35)
           alpha ~ 1 + ap2 + ap3 + ap4 + ap5 + ap6",
#   pre_start = FALSE,
  itype = "2pl",
  method = "vb", iter = 15000, tol_rel_obj = 5e-4
  )
 
fit = theta2::theta2(
  data = data,
  model = "theta = c(1:35)
           alpha ~ 1 + ap2 + (1 + ap2|z_fac)",
#   pre_start = FALSE,
  method = "vb", iter = 15000, tol_rel_obj = 5e-4
  )

model = "theta = c(1:35)
        alpha ~ 1 + ap2 + (1 + ap2|z_fac)"
method = "vb"; iter = 10000; tol_rel_obj = 2e-4
exploratory = FALSE
weights = NULL
library(devtools)
load_all()

model = "theta = c(1:25)
         theta ~ 0
         alpha ~ 1 + ap2 + (1|z_fac)
         delta ~ 1"
method = "vb"; iter = 5000; tol_rel_obj = 2e-4
library(devtools)
load_all()
weights = NULL
exploratory = FALSE


cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 35, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

rstan::summary(fit, pars = "alpha_r_l")$summary
exp(rstan::summary(fit, pars = "alpha_r_l")$summary)
exp(c(0.8, 1.4, 0.7, 0.2, 0.4, -0.8))

fit@model_pars

cor(sim_data$zeta, (rstan::summary(fit, pars = "aeta_l")$summary[,1]))
plot(sim_data$zeta, (rstan::summary(fit, pars = "aeta_l")$summary[,1]))

rstan::summary(fit, pars = "aeta_l_sd")$summary

fit@model_pars

omega = matrix(rstan::summary(fit, pars = "Omega_a")$summary[,1], nrow = 2)
tau = diag(rstan::summary(fit, pars = "tau_a")$summary[,1])
vcov = tau %*% omega %*% tau
vcov


aeta_c = matrix(rstan::summary(fit, pars = "aeta_c")$summary[,1], byrow = TRUE, ncol = 2)

cor((aeta_c[,1]), sim_data$int_slope[,1])
cor((aeta_c[,2]), sim_data$int_slope[,2])

plot((aeta_c[,1]), sim_data$int_slope[,1])
plot((aeta_c[,2]), sim_data$int_slope[,2])


fit@model_pars

rstan::summary(fit, pars = "alpha")$summary
exp(rstan::summary(fit, pars = "alpha")$summary[,1])

cor(sim_data$alpha[1,], exp(rstan::summary(fit, pars = "alpha")$summary[,1]))
plot(sim_data$alpha[1,], exp(rstan::summary(fit, pars = "alpha")$summary[,1]))

# What does a simple single fixed effect give with this
# new alpha parameterization?
# two groups + mean alpha

