# INIRT with univariate theta dimension
devtools::install(dependencies = FALSE)
# library(devtools)
# Rcpp::compileAttributes()
# load_all()
# Test 1: simplest possible settings
n = 100
ncat = 3
j = 25
d = 1
k = 0
uk = 0
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = rnorm(n)))
b_alpha = c(0.8, 1.2)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 1.3
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
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
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                                                  z[i, ] %*% zeta
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + as.vector(b_alpha%*%a_design[i,]))*(theta[i, ])) - (delta[jj, ] + as.vector(b_delta%*%d_design[i,]))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, a_design)
colnames(data) = c(paste0("x", 1:j), paste0("ap", 1:2))
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
# fit = inirt::inirt(data = fit_data$data, item_id = fit_data$item_id, model = fit_data$model, predictors = fit_data$predictors, 
#     dims = fit_data$dims, h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, init = fit_data$init)


# data = fit_data$data
# fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, method = "vb", iter = 5000, tol_rel_obj = 2e-4)

# new interface
# library(stringr)
# data = fit_data$data
# model = "theta = c(1:25)
#          theta ~ 0
#          alpha ~ 1
#          delta ~ 1"
# # mod = parse_model(model ,data)
# method = "vb"
# iter = 5000
# tol_rel_obj = 2e-4
# exploratory = FALSE
# weights = NULL
# source("R/parse_model.R")
# source("R/parse_regression_eq.R")
# source("R/parse_theta_eq.R")
# library(stringr)
# fit = inirt::inirt(
#   data = data,
#   model = "theta = c(1:25)
#            theta ~ 0
#            alpha ~ 1
#            delta ~ 1",
#   method = "vb", iter = 5000, tol_rel_obj = 2e-4
#   )
# method = "hmc"; iter = 300; warmup = 150; chains = 1
# exploratory = FALSE
# weights = NULL
# source("R/parse_model.R")
# source("R/parse_regression_eq.R")
# source("R/parse_theta_eq.R")
# library(stringr)

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
fit = theta2::theta2(
  data = data,
  model = "theta = c(1:25)
           theta ~ 0
           alpha ~ 1 + ap2
           delta ~ 1",
  method = "vb", iter = 5000, tol_rel_obj = 2e-4
  )



fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta_l")$summary

fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
plot(sim_data$alpha[1,], (rstan::summary(fit, pars = "alpha")$summary[,1]))

exp(rstan::summary(fit, pars = "alpha_r_l")$summary[,1])

# looked correct on this run


rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
rstan::summary(fit, pars = "alpha_l")$summary[,1]

cor(
cbind(
  as.vector(sim_data$b_alpha + sim_data$alpha),
  as.vector(exp(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1]))
)
)

plot(sim_data$alpha[1,], exp(rstan::summary(fit, pars = "alpha")$summary[,1] + 
  rstan::summary(fit, pars = "alpha_r_l")$summary[,1]))

exp(rstan::summary(fit, pars = "alpha_r_l")$summary[,1])

cbind(
  sim_data$alpha[1,],
  rstan::summary(fit, pars = "alpha")$summary[,1] + rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
)

rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])
plot(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]

library(bayesplot)
color_scheme_set("blue")
mcmc_trace(fit, pars = c("alpha_l[2]"))



# --------------------------
devtools::install(dependencies = FALSE)

# Test 2: beta regression, fixed effects only
n = 800
ncat = 4
j = 25
d = 1
k = 2
uk = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n)))
b_alpha = 0.8
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 1.3
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(1.0, 0.4)
predictors = list(c(1, 2))
start_index = 1
beta_dstart = numeric(d)
beta_dend = numeric(d)
for(dd in 1:d) {
  beta_dstart[dd] = start_index
  beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
  start_index = start_index + length(predictors[[dd]])
}
beta_mat = matrix(0, nrow = k, ncol = d)
index = 1
for(dd in 1:d) {
  for(i in beta_dstart[dd]:beta_dend[dd]) {
    beta_mat[i, dd] = beta[index]
    index = index + 1
  }
}
x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                                                  z[i, ] %*% zeta
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat)) - (delta[jj, ] + b_delta*d_design[i,])))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x)
colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2))
dims = 1
item_id = 1:j
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
fit_data = list(data = data, item_id = item_id, predictors = predictors)

rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

data = fit_data$data
fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
  method = "vb", iter = 5000, tol_rel_obj = 2e-4)
  method = "vb"
data = fit_data$data
iter = 5000
tol_rel_obj = 2e-4
exploratory = FALSE
weights = NULL
source("R/parse_model.R")
source("R/parse_regression_eq.R")
source("R/parse_theta_eq.R")
library(stringr)
fit = inirt::inirt(
  data = data,
  model = "theta = c(1:25)
           theta ~ pred1 + pred2
           alpha ~ 1
           delta ~ 1",
  method = "vb", iter = 5000, tol_rel_obj = 2e-4
  )

# run
data = fit_data$data
fit = inirt::inirt(
  data = data,
  model = "theta = c(1:25)
           theta ~ pred1 + pred2
           alpha ~ 1
           delta ~ 1",
  method = "hmc", iter = 20, chains = 1
  )


fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta")$summary
rstan::summary(fit, pars = "beta_l")$summary
# this one worked but fit the wrong model! - it should have two predictors
# but does not




fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]



# --------------------------
devtools::install(dependencies = FALSE)

# Test 3: beta regression, fixed effects only, alpha fixed effects, delta fixed effects
n = 800
ncat = 4
j = 25
d = 1
k = 2
uk = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_alpha = c(0.8, -0.6)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_delta = c(-1.3, 0.4)
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(1.0, 0.4)
predictors = list(c(1, 2))
start_index = 1
beta_dstart = numeric(d)
beta_dend = numeric(d)
for(dd in 1:d) {
  beta_dstart[dd] = start_index
  beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
  start_index = start_index + length(predictors[[dd]])
}
beta_mat = matrix(0, nrow = k, ncol = d)
index = 1
for(dd in 1:d) {
  for(i in beta_dstart[dd]:beta_dend[dd]) {
    beta_mat[i, dd] = beta[index]
    index = index + 1
  }
}
x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                                                  z[i, ] %*% zeta
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha %*% a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat)) - (delta[jj, ] + as.vector(b_delta %*% d_design[i,]))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x)
colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2))
dims = 1
item_id = 1:j
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
fit_data = list(data = data, item_id = item_id, predictors = predictors, alpha_data = a_design, 
  delta_data = d_design)

rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

data = fit_data$data
alpha_data = fit_data$alpha_data
delta_data = fit_data$delta_data
fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
  structural_design = list(alpha = alpha_data, delta = delta_data), method = "vb", 
  iter = 7000, tol_rel_obj = 2e-4)

fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta")$summary
rstan::summary(fit, pars = "beta_l")$summary





fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]





# --------------------------
devtools::install(dependencies = FALSE)

# Test 4: beta regression, fixed effects + two random effects (each with 20 dummies), alpha fixed effects, delta fixed effects
n = 200
ncat = 4
j = 25
d = 1
k = 2
uk = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_alpha = c(0.8, -0.6)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_delta = c(-1.3, 0.4)
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(1.0, 0.4)
predictors = list(c(1, 2))
start_index = 1
beta_dstart = numeric(d)
beta_dend = numeric(d)
for(dd in 1:d) {
  beta_dstart[dd] = start_index
  beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
  start_index = start_index + length(predictors[[dd]])
}
beta_mat = matrix(0, nrow = k, ncol = d)
index = 1
for(dd in 1:d) {
  for(i in beta_dstart[dd]:beta_dend[dd]) {
    beta_mat[i, dd] = beta[index]
    index = index + 1
  }
}
x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
Lz = 10
z = matrix(rep(diag(10), each = 20), nrow = Lz*20) # for n = 800
z = cbind(z, z[gtools::permute(1:nrow(z)), ])
zeta_sd = 2
zeta = rnorm(Lz*2, 0, sd = zeta_sd)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha %*% a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat + z[i, ] %*% zeta)) - (delta[jj, ] + as.vector(b_delta %*% d_design[i,]))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x, z)
colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2), paste0("z", 1:(2*Lz)))
ranef_id = c(rep(1, ncol(z)/2), rep(2, ncol(z)/2))
dims = 1
item_id = 1:j
predictors_ranef = list(1:ncol(z) + j + k)

sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, 
  beta = beta, theta = theta, zeta = zeta)
fit_data = list(data = data, item_id = item_id, predictors = predictors, 
  predictors_ranef = predictors_ranef, alpha_data = a_design, delta_data = d_design)

rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

data = fit_data$data
alpha_data = fit_data$alpha_data
delta_data = fit_data$delta_data
fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
  predictors_ranef = list(28:107), ranef_id = c(rep(1, 40), rep(2, 40)),
  structural_design = list(alpha = alpha_data, delta = delta_data), method = "vb", 
  iter = 7000, tol_rel_obj = 2e-4)

fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta")$summary
rstan::summary(fit, pars = "beta_l")$summary


cor(sim_data$zeta[1:40], rstan::summary(fit, pars = "zeta")$summary[1:40,1])
cor(sim_data$zeta[41:80], rstan::summary(fit, pars = "zeta")$summary[41:80,1])
plot(sim_data$zeta[1:40], rstan::summary(fit, pars = "zeta")$summary[1:40,1])
plot(sim_data$zeta[1:40], rstan::summary(fit, pars = "zeta")$summary[1:40,1])
rstan::summary(fit, pars = "zeta_l_sd")$summary


fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]



# --------------------------
devtools::install(dependencies = FALSE)

# Test 4: beta regression, fixed effects, alpha fixed effects + one random effects (with 20 levels), delta fixed effects
n = 800
ncat = 4
j = 25
d = 1
k = 2
uk = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_alpha = c(0.8, -0.6)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_delta = c(-1.3, 0.4)
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(1.0, 0.4)
predictors = list(c(1, 2))
start_index = 1
beta_dstart = numeric(d)
beta_dend = numeric(d)
for(dd in 1:d) {
  beta_dstart[dd] = start_index
  beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
  start_index = start_index + length(predictors[[dd]])
}
beta_mat = matrix(0, nrow = k, ncol = d)
index = 1
for(dd in 1:d) {
  for(i in beta_dstart[dd]:beta_dend[dd]) {
    beta_mat[i, dd] = beta[index]
    index = index + 1
  }
}
x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
Lz = 40
z = matrix(rep(diag(40), each = 20), nrow = Lz*20)
# z = cbind(z, z[gtools::permute(1:nrow(z)), ])
zeta_sd = 2
# zeta = rnorm(Lz*2, 0, sd = zeta_sd)
zeta = rnorm(Lz, 0, sd = zeta_sd)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha %*% a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat)) - (delta[jj, ] + as.vector(b_delta %*% d_design[i,]) + z[i, ] %*% zeta)))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x, z)
# colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2), paste0("z", 1:(2*Lz)))
colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2), paste0("z", 1:(Lz)))
# ranef_id = c(rep(1, ncol(z)/2), rep(2, ncol(z)/2))
ranef_id = c(rep(1, ncol(z)))
dims = 1
item_id = 1:j
predictors_ranef = list(1:ncol(z) + j + k)

sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, 
  beta = beta, theta = theta, zeta = zeta)
fit_data = list(data = data, item_id = item_id, predictors = predictors, 
  predictors_ranef = predictors_ranef, alpha_data = a_design, delta_data = d_design)

rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

data = fit_data$data
alpha_data = fit_data$alpha_data
delta_data = fit_data$delta_data

# data = data; item_id = 1:25; dims = 1; predictors = list(c(26, 27)); 
# structural_design = list(alpha = alpha_data, delta = delta_data);
# structural_design_ranef = list(a_predictors_ranef = data[, 28:67], a_ranef_id = c(rep(1, 40)));
# method = "vb"; iter = 7000; tol_rel_obj = 2e-4

fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
  structural_design = list(alpha = alpha_data, delta = delta_data),
  structural_design_ranef = list(d_predictors_ranef = data[, 28:67], d_ranef_id = c(rep(1, 40))),
  method = "vb", iter = 7000, tol_rel_obj = 2e-4)
# fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
#   predictors_ranef = list(28:107), ranef_id = c(rep(1, 40), rep(2, 40)),
#   structural_design = list(alpha = alpha_data, delta = delta_data), method = "vb", 
#   iter = 7000, tol_rel_obj = 2e-4)

fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta")$summary
rstan::summary(fit, pars = "beta_l")$summary

rstan::summary(fit, pars = "deta_l")$summary[,1]
rstan::summary(fit, pars = "deta_l_sd")$summary[,1]



cor(sim_data$zeta[1:40], rstan::summary(fit, pars = "deta_l")$summary[,1])
plot(sim_data$zeta[1:40], rstan::summary(fit, pars = "deta_l")$summary[,1])

rstan::summary(fit, pars = "deta_l_sd")$summary[,1]

cor(sim_data$zeta[1:40], rstan::summary(fit, pars = "zeta")$summary[1:40,1])
cor(sim_data$zeta[41:80], rstan::summary(fit, pars = "zeta")$summary[41:80,1])
plot(sim_data$zeta[1:40], rstan::summary(fit, pars = "zeta")$summary[1:40,1])
plot(sim_data$zeta[1:40], rstan::summary(fit, pars = "zeta")$summary[1:40,1])
rstan::summary(fit, pars = "zeta_l_sd")$summary


fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]






# --------------------------
devtools::install(dependencies = FALSE)

# Test 5: beta regression, fixed effects + two correlated random effects (intercept, slope), alpha fixed effects, delta fixed effects
n = 800
ncat = 4
j = 25
d = 1
k = 2
uk = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_alpha = c(0.8, -0.6)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_delta = c(-1.3, 0.4)
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(1.0, 0.4)
predictors = list(c(1, 2))
start_index = 1
beta_dstart = numeric(d)
beta_dend = numeric(d)
for(dd in 1:d) {
  beta_dstart[dd] = start_index
  beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
  start_index = start_index + length(predictors[[dd]])
}
beta_mat = matrix(0, nrow = k, ncol = d)
index = 1
for(dd in 1:d) {
  for(i in beta_dstart[dd]:beta_dend[dd]) {
    beta_mat[i, dd] = beta[index]
    index = index + 1
  }
}
x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
Lz = 40
sigma = matrix(c(4, 3, 3, 3.5), nrow = 2)
int_slope = mvtnorm::rmvnorm(Lz, mean = c(0, 0), sigma = sigma)
cov(int_slope)
z_c_int = matrix(rep(diag(40), each = 20), nrow = Lz*20)
zc_is = z_c_int %*% int_slope[,1]
z_c_slope = matrix(rep(diag(40), each = 20), nrow = Lz*20)
z_c_slope[z_c_slope == 1] = runif(40*20)
zc_is = zc_is + z_c_slope %*% int_slope[,2]
z_c = cbind(z_c_int, z_c_slope)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha %*% a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat + zc_is[i,])) - (delta[jj, ] + as.vector(b_delta %*% d_design[i,]))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x, z_c)
colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2), paste0("z", 1:(2*Lz)))
dims = 1
item_id = 1:j
predictors_ranef = list(1:ncol(z_c) + j + k)

sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, 
  beta = beta, theta = theta, int_slope = int_slope)
fit_data = list(data = data, item_id = item_id, predictors = predictors, 
  predictors_ranef = predictors_ranef, alpha_data = a_design, delta_data = d_design)

rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

data = fit_data$data
alpha_data = fit_data$alpha_data
delta_data = fit_data$delta_data
fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
  predictors_ranef_corr = list(28:107), n_pranef_cor = 2,
  structural_design = list(alpha = alpha_data, delta = delta_data), method = "vb", 
  iter = 7000, tol_rel_obj = 2e-4)

fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta")$summary
rstan::summary(fit, pars = "beta_l")$summary


betas = matrix(rstan::summary(fit, pars = "zeta_c")$summary[,1], ncol = 2, byrow = TRUE)
sim_data$int_slope
cor(sim_data$int_slope[,1], betas[,1])
cor(sim_data$int_slope[,2], betas[,2])

omega = matrix(rstan::summary(fit, pars = "Omega")$summary[,1], nrow = 2)
tau = diag(rstan::summary(fit, pars = "tau")$summary[,1])
vcov = tau %*% omega %*% tau
vcov

fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]







# --------------------------
devtools::install(dependencies = FALSE)

# Test 6: beta regression, fixed effects + two correlated random effects (intercept, slope), alpha fixed effects + same random effects, delta fixed effects + same random effects
n = 800
ncat = 4
j = 25
d = 1
k = 2
uk = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_alpha = c(0.8, -0.6)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n), x2 = runif(n,-1,1)))
b_delta = c(-1.3, 0.4)
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(1.0, 0.4)
predictors = list(c(1, 2))
start_index = 1
beta_dstart = numeric(d)
beta_dend = numeric(d)
for(dd in 1:d) {
  beta_dstart[dd] = start_index
  beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
  start_index = start_index + length(predictors[[dd]])
}
beta_mat = matrix(0, nrow = k, ncol = d)
index = 1
for(dd in 1:d) {
  for(i in beta_dstart[dd]:beta_dend[dd]) {
    beta_mat[i, dd] = beta[index]
    index = index + 1
  }
}
x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
Lz = 40
sigma = matrix(c(4, 3, 3, 3.5), nrow = 2)
int_slope = mvtnorm::rmvnorm(Lz, mean = c(0, 0), sigma = sigma)
cov(int_slope)
z_c_int = matrix(rep(diag(40), each = 20), nrow = Lz*20)
zc_is = z_c_int %*% int_slope[,1]
z_c_slope = matrix(rep(diag(40), each = 20), nrow = Lz*20)
z_c_slope[z_c_slope == 1] = runif(40*20)
zc_is = zc_is + z_c_slope %*% int_slope[,2]
z_c = cbind(z_c_int, z_c_slope)

sigma_a = matrix(c(4, 3, 3, 3.5), nrow = 2)
int_slope_a = mvtnorm::rmvnorm(Lz, mean = c(0, 0), sigma = sigma_a)
zc_is_a = z_c_int %*% int_slope_a[,1] + z_c_slope %*% int_slope_a[,2]

data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {                                                                                   #  + zc_is[i,]
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha %*% a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat)) - (delta[jj, ] + as.vector(b_delta %*% d_design[i,] + zc_is_a[i,]))))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x, z_c)
colnames(data) = c(paste0("x", 1:j), paste0("pred", 1:2), paste0("z", 1:(2*Lz)))
dims = 1
item_id = 1:j
predictors_ranef = list(1:ncol(z_c) + j + k)

sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, 
  beta = beta, theta = theta, int_slope = int_slope, int_slope_a = int_slope_a)
fit_data = list(data = data, item_id = item_id, predictors = predictors, 
  predictors_ranef = predictors_ranef, alpha_data = a_design, delta_data = d_design)

rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()

data = fit_data$data
alpha_data = fit_data$alpha_data
delta_data = fit_data$delta_data
d_predictors_ranef_corr = data[, 28:107]

# data = data; item_id = 1:25; dims = 1; predictors = list(c(26, 27)) 
# predictors_ranef_corr = list(28:107); n_pranef_cor = 2;
# structural_design = list(alpha = alpha_data, delta = delta_data); 
# structural_design_ranef = list(a_predictors_ranef_corr = a_predictors_ranef_corr, a_n_pranef_cor = 2)
# method = "vb"; iter = 7000; tol_rel_obj = 2e-4

fit = inirt::inirt(data = data, item_id = 1:25, dims = 1, predictors = list(c(26, 27)), 
 #predictors_ranef_corr = list(28:107), n_pranef_cor = 2,
  structural_design = list(alpha = alpha_data, delta = delta_data), 
  structural_design_ranef = list(d_predictors_ranef_corr = d_predictors_ranef_corr, d_n_pranef_cor = 2),
  method = "vb", iter = 7000, tol_rel_obj = 2e-4, vb_algorithm = "meanfield")

fit@model_pars
rstan::summary(fit, pars = "delta_r_l")$summary
rstan::summary(fit, pars = "delta_l")$summary
rstan::summary(fit, pars = "alpha_r_l")$summary
rstan::summary(fit, pars = "alpha_l")$summary

rstan::summary(fit, pars = "beta")$summary
rstan::summary(fit, pars = "beta_l")$summary


betas = matrix(rstan::summary(fit, pars = "zeta_c")$summary[,1], ncol = 2, byrow = TRUE)
sim_data$int_slope
cor(sim_data$int_slope[,1], betas[,1])
cor(sim_data$int_slope[,2], betas[,2])

aetas = matrix(rstan::summary(fit, pars = "aeta_c")$summary[,1], ncol = 2, byrow = TRUE)
sim_data$int_slope
cor(sim_data$int_slope_a[,1], aetas[,1])
cor(sim_data$int_slope_a[,2], aetas[,2])

detas = matrix(rstan::summary(fit, pars = "deta_c")$summary[,1], ncol = 2, byrow = TRUE)
sim_data$int_slope
cor(sim_data$int_slope_a[,1], detas[,1])
cor(sim_data$int_slope_a[,2], detas[,2])

omega_a = matrix(rstan::summary(fit, pars = "Omega_a")$summary[,1], nrow = 2)
tau_a = diag(rstan::summary(fit, pars = "tau_a")$summary[,1])
vcov_a = tau_a %*% omega_a %*% tau_a
vcov_a

omega_d = matrix(rstan::summary(fit, pars = "Omega_d")$summary[,1], nrow = 2)
tau_d = diag(rstan::summary(fit, pars = "tau_d")$summary[,1])
vcov_d = tau_d %*% omega_d %*% tau_d
vcov_d

omega = matrix(rstan::summary(fit, pars = "Omega")$summary[,1], nrow = 2)
tau = diag(rstan::summary(fit, pars = "tau")$summary[,1])
vcov = tau %*% omega %*% tau
vcov

fit@model_pars
cor(sim_data$alpha[1,], rstan::summary(fit, pars = "alpha")$summary[,1])
rstan::summary(fit, pars = "alpha_r_l")$summary[,1]
mean(sim_data$b_alpha)
mean(rstan::summary(fit, pars = "alpha_r_l")$summary[,1] + rstan::summary(fit, pars = "alpha")$summary[,1])


plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta)


dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 25, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])

mean(dest)
mean(sim_data$delta[,-1])
rstan::summary(fit, pars = c("delta_r_l"))$summary[,1]
dest[,1] + sim_data$b_delta
sim_data$delta[,2]

