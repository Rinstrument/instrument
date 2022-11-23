# INIRT with univariate theta dimension
devtools::install(dependencies = FALSE)
n = 800
j = 25
d = 1
k = 6
uk = 6
ncategi = c(rep(4, 25))
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
beta = c(-0.3, 0.2, 1, -0.6, -0.4, 0.6)
predictors = list(c(1:3,4:6))
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
lz_cor = TRUE
lz_both = TRUE
if(lz_both) {
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
  Lz = 40
  z = matrix(rep(diag(40), each = 20), nrow = Lz*20) # for n = 800
  zeta_sd = 0.7
  zeta = rnorm(Lz, 0, sd = zeta_sd)
} else {
  if(lz_cor) {
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
  } else {
    Lz = 40
    z = matrix(rep(diag(40), each = 20), nrow = Lz*20) # for n = 800
    zeta_sd = 0.7
    zeta = rnorm(Lz, 0, sd = zeta_sd)
  }
}
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                                                  z[i, ] %*% zeta
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat + zc_is[i,] + z[i, ] %*% zeta)) - (delta[jj, ] + b_delta*d_design[i,])))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
data = cbind(data, x, z_c, z)
colnames(data) = c(paste0("x", 1:j), paste0("k", 1:k), paste0("z", 1:(2*Lz)), paste0("z", 1:(Lz)))
for(dd in 1:d) {
  predictors[[dd]] = predictors[[dd]] + j
}
predictors_ranef = list(1:ncol(z) + ncol(z_c) + j + k)
predictors_ranef_corr = list(1:ncol(z_c) + j + k)
Lzeta_cor = ncol(z_c)
cor_z_item_ind = rep(1:Lz, 2)
cor_z_item_elem_ind = rep(1:2, each = Lz)
ranef_id = rep(1, ncol(z)) # ranef_id = 0 #rep(1, ncol(z))
dims = 1
item_id = 1:j
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta, 
     int_slope = int_slope, zeta = zeta) # zeta = zeta, zeta_sd = zeta_sd,
# fit_data = list(data = data, item_id = item_id, model = NULL, predictors = predictors, predictors_ranef = predictors_ranef, ranef_id = ranef_id, 
#     dims = dims, h2_dims = 0, h2_dim_id = NULL, structural_design = list(alpha = a_design, delta = d_design), 
#     method = "vb", weights = NULL, tol_rel_obj = 0.0002, iter = 5e3, init = "random")
fit_data = list(data = data, item_id = item_id, model = NULL, predictors = predictors, predictors_ranef = predictors_ranef, ranef_id = ranef_id, 
    predictors_ranef_corr = predictors_ranef_corr, n_pranef_cor = 2,
    dims = dims, h2_dims = 0, h2_dim_id = NULL, structural_design = list(alpha = a_design, delta = d_design), 
    method = "vb", weights = NULL, tol_rel_obj = 0.0002, iter = 5e3, init = "random")
data = data; item_id = item_id; model = NULL; predictors = predictors; predictors_ranef = predictors_ranef; ranef_id = ranef_id; 
predictors_ranef_corr = predictors_ranef_corr; n_pranef_cor = 2;
dims = dims; h2_dims = 0; h2_dim_id = NULL; structural_design = list(alpha = a_design, delta = d_design); 
method = "vb"; weights = NULL; tol_rel_obj = 0.0002; iter = 5e3; init = "random";
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()
# fit = inirt::inirt(data = fit_data$data, item_id = fit_data$item_id, model = fit_data$model, predictors = fit_data$predictors, 
#     predictors_ranef = fit_data$predictors_ranef, ranef_id = fit_data$ranef_id, dims = fit_data$dims, 
#     h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, 
#     init = fit_data$init)

fit = inirt::inirt(data = fit_data$data, item_id = fit_data$item_id, model = fit_data$model, predictors = fit_data$predictors, 
    predictors_ranef = fit_data$predictors_ranef, ranef_id = fit_data$ranef_id, predictors_ranef_corr = fit_data$predictors_ranef_corr, n_pranef_cor = fit_data$n_pranef_cor, dims = fit_data$dims, 
    h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
    method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, 
    init = fit_data$init)

# fit = inirt::inirt(data, dims = d, method = "hmc", chains = 1, iter = 300, init = 0)
Omega = matrix(rstan::summary(fit, pars = "Omega")$summary[, "mean"], 2, 2, byrow = TRUE)
diag_tau = diag(rstan::summary(fit, pars = "tau")$summary[,"mean"])
diag_tau %*% Omega %*% diag_tau
cov(sim_data$int_slope)

betas = matrix(rstan::summary(fit, pars = "zeta_c")$summary[,1], ncol = 2, byrow = TRUE)
sim_data$int_slope
cor(sim_data$int_slope[,1], betas[,1])
cor(sim_data$int_slope[,2], betas[,2])

cor(sim_data$zeta, rstan::summary(fit, pars = "zeta")$summary[,1])

rstan::summary(fit, pars = "alpha")$summary[,1]
fit@model_pars
sim_data$alpha
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

rstan::summary(fit, pars = "beta_l")$summary[,1]
sim_data$beta


rstan::summary(fit, pars = "zeta_l_sd")$summary[,]
rstan::summary(fit, pars = "zeta_l")$summary[,1]
cor(sim_data$zeta, rstan::summary(fit, pars = "zeta_l")$summary[,1])
