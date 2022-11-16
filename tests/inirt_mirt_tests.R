# INIRT higher-order IRT tests
devtools::install()
n = 800
d = 3
j = 40
k = 6
uk = 6
ncategi = c(rep(4, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n)))
b_alpha = 0.8
alpha_dominant = list(1:15, 16:27, 28:40)
for(dd in 1:d) {
  alpha[dd, alpha_dominant[[dd]]] = sort(runif(length(alpha_dominant[[dd]]), 1.7, 3.0))
  alpha[dd, setdiff(unlist(alpha_dominant), alpha_dominant[[dd]])] = sort(runif(length(unlist(alpha_dominant[-dd])), 0.2, 1.0))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = -1.1
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = c(-0.3, 0.2, 1, -0.6, -0.4, 0.6)
predictors = list(1:2, 3:4, 5:6)
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
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(theta[i, ] + x[i, ] %*% beta_mat)) - (delta[jj, 1:ncategi[jj]] + b_delta*d_design[i,])))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
  }
}
apply(data, 2, table)
remove_gaps = function(x) {
  ord = order(x); vec = sort(x)
  old = unique(vec); replace = 1:length(unique(vec))
  names(replace) = old; names(vec) = vec
  new = replace[names(vec)]; names(new) = NULL
  return(new[ord])
}
data = apply(data, 2, remove_gaps)
apply(data, 2, table)
data = cbind(data, x)
colnames(data) = c(paste0("x", 1:j), paste0("z", 1:k))
for(dd in 1:d) {
  predictors[[dd]] = predictors[[dd]] + j
}
dims = 3
h2_dims = 0
h2_dim_id = NULL
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
fit_data = list(data = data, model = NULL, predictors = predictors, dims = dims, h2_dims = h2_dims, h2_dim_id = h2_dim_id, 
    structural_design = list(alpha = a_design, delta = d_design), method = "vb", weights = NULL, 
    tol_rel_obj = 0.0002, iter = 5e3, init = "random")
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()
fit = inirt::inirt(data = fit_data$data, model = fit_data$model, predictors = fit_data$predictors, dims = fit_data$dims, 
    h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
    method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, 
    init = fit_data$init)

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
