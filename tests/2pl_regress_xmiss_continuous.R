library(rstan)
rm(list = ls())
devtools::install()
n = 800
j = 25
ncategi = c(rep(4, 25))
ncateg_max = max(ncategi)
k = 6
uk = 6
d = 1
# p = 2*j + n + k
# alpha = exp(runif(j, 0, 1.5)) - 0.2
alpha = matrix(0, d, j)
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
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
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(sum(alpha[, jj]*(theta[i, ] + x[i, ] %*% beta_mat)) - delta[jj, 1:ncategi[jj]]))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
  }
}
# randomly place some missing values
x_to_remove = sample(1:prod(dim(x)), size = 25)
x_removed_true = x[x_to_remove]
x[x_to_remove] = NA
data = cbind(data, x)
# coef(mirt(data, 1))
colnames(data) = c(paste0("x", 1:j), paste0("k", 1:uk))
# fit = inirt(data, method = "hmc", iter = 2000, chains = 1)
# rm(n, j, alpha, delta, theta, i, jj)
for(dd in 1:d) {
  predictors[[dd]] = predictors[[dd]] + j
}
# fit = inirt::inirt(data, predictors = predictors, dims = d, method = "vb", tol_rel_obj = 0.0001, iter = 1e4, init = 0)
fit = inirt::inirt(data, model = NULL, predictors = predictors, dims = 1, method = "vb", weights = NULL, tol_rel_obj = 0.0005, iter = 1e4, init = "random")
fit = inirt::inirt(data, model = NULL, predictors = predictors, dims = 1, method = "hmc", chains = 1, iter = 100, init = "random")

summary(fit, pars = "alpha")$summary[,"mean"]
cor(summary(fit, pars = "alpha")$summary[,"mean"], alpha[1,])
dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,"mean"], nrow = 25, byrow = TRUE)
cor(dest[,1], delta[,2])

cor(summary(fit, pars = c("theta"))$summary[,"mean"], theta[,1])
summary(fit, pars = c("x_l"))
x_removed_true
cor(summary(fit, pars = c("x_l"))$summary[,"mean"], x_removed_true)
ords = unique(as.vector((standata$reg_miss)))[-1]
cor(summary(fit, pars = c("x_l"))$summary[,"mean"], x_removed_true[ords])

summary(fit, pars = c("beta"))
beta

x2 = x
x2[is.na(x)] = 1:20
x2[!is.na(x)] = 0

x_miss_id = 1:20
reg_miss = is.na(x) * 1
reg_miss = as.vector(t(reg_miss))
reg_miss[reg_miss == 1] = x_miss_id
reg_miss = matrix(reg_miss, nrow = nrow(x), byrow = TRUE)
ord = reg_miss[reg_miss != 0]
