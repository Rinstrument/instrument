library(rstan)
rm(list = ls())
devtools::install()
n = 800
d = 3
j = 20*d
ncategi = c(rep(4, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
alpha_dstart = c(1, 21, 41)
alpha_dend = c(20, 40, 60)
for(dd in 1:d) {
  alpha[dd, alpha_dstart[dd]:alpha_dend[dd]] = sort(runif(j/d, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta_g = rnorm(n, 0, 3)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 0.1)
}
lambda_ind = rep(1:d, each = j/d)
lambda = c(0.8, 1.3, 1.6)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(sum(alpha[, jj]*(lambda[lambda_ind[jj]]*theta_g[i] + theta[i, lambda_ind[jj]])) - delta[jj, 1:ncategi[jj]]))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
  }
}
colnames(data) = c(paste0("x", 1:j))

dims = 3
h2_dims = 1
method = "vb"
h2_dim_id = list(1:20, 21:40, 41:60)
fit = inirt::inirt(data, dims = dims, h2_dims = h2_dims, method = method, weights = NULL, 
    h2_dim_id = h2_dim_id, tol_rel_obj = 0.0002, iter = 1e4, init = "random")
# fit = inirt::inirt(data, model = NULL, predictors = predictors, dims = 1, method = "hmc", chains = 1, iter = 100, init = "random")

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
