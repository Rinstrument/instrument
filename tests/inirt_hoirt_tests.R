# INIRT higher-order IRT tests
devtools::install(dependencies = FALSE)
library(rstan)
stanc(file = "./inst/stan/theta2_soirt.stan", verbose = TRUE)
n = 1000
d = 4
j = 20*d
ncat = 2
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
k = 0
uk = 0
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n)))
b_alpha = 1
# alpha_dstart = c(1, 6, 11, 16)
# alpha_dend = c(5, 10, 15, 20)
alpha_dstart = c(1, 21, 41, 61)
alpha_dend = c(20, 40, 60, 80)
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
theta_g = rnorm(n, 0, 5) + 100 # seq(-2, 2, length.out = n) #rnorm(n, 0, 2) # sqrt(1.5)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 2.5) #0.1??????
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
lambda = c(0.6, 0.5, 0.4, -0.3)
f_eq = matrix(nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    f_eq[i, jj] = (lambda[lambda_ind[jj]]*(theta_g[i])) + theta[i, lambda_ind[jj]] # do I need some jitter on this???
  }
}
# summary(lm(f_eq[,4] ~ theta_g))
# tg = (theta_g - mean(theta_g))/sd(theta_g)
# summary(lm(f_eq[,4] ~ tg))
# rescale f_eq columns
for(jj in 1:j) {
  f_eq[, jj] = (f_eq[, jj] - mean(f_eq[, jj])) / sd(f_eq[, jj])
}
tg = (theta_g - mean(theta_g))/sd(theta_g)
summary(lm(f_eq[,1] ~ tg))
cor(f_eq[,1], tg)
simmed_cors = c(cor(f_eq[,1], tg), cor(f_eq[,21], tg), cor(f_eq[,41], tg), cor(f_eq[,61], tg))
beta = 1.2
x = matrix(data = runif(n,-1,1), nrow = n, ncol = 1)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) { #                                                                      + x[i, ] %*% beta_mat
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*( f_eq[i, jj] + beta*x[i,])) - (delta[jj, 1:ncategi[jj]] + b_delta*d_design[i,])))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
  }
}
apply(data, 2, table)
data = cbind(data, x)
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
colnames(data) = c(paste0("x", 1:j), "p1") #, paste0("z", 1:k)
# for(dd in 1:1) {
#   predictors[[dd]] = predictors[[dd]] + j
# }
# dims = 3
# h2_dims = 1
# h2_dim_id = list(1:20, 21:40, 41:60)
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta, lambda = lambda,
    theta_g = theta_g, simmed_cors = simmed_cors, beta = beta)
fit_data = list(data = data)
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()
# fit = inirt::inirt(data = fit_data$data, model = fit_data$model, predictors = fit_data$predictors, dims = fit_data$dims, 
#     h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, 
#     init = fit_data$init)

library(devtools)
library(Rcpp)
# delfiles <- dir(pattern = "*_pattern.csv")
# file.remove(file.path(mydir, delfiles))
compileAttributes()
load_all()
data = fit_data$data
model = "thetag = theta1 + theta2 + theta3 + theta4
         theta1 = c(1:20)
         theta2 = c(21:40)
         theta3 = c(41:60)
         theta4 = c(61:80)"
         # thetag ~ p1"
# itype = "2pl"
# method = "vb"
# iter = 10000
# tol_rel_obj = 1e-4
# exploratory = FALSE
# weights = NULL

#---
# itype = "2pl"
# method = "hmc"
# iter = 500
# warmup = 300
# chains = 1
# fweights = NULL
# cores = 1

# data = fit_data$data
# colnames(data)
# model = "thetag = theta1 + theta2 + theta3 + theta4
#          theta1 = c(1:20)
#          theta2 = c(21:40)
#          theta3 = c(41:60)
#          theta4 = c(61:80)"
fit = theta2::theta2(data = data, model = model, itype = "2pl", method = "hmc", 
  iter = 50, warmup = 30, chains = 1, cores = 1)

library(rstan)

sim_data$theta
sim_data$theta_g
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta_g)
plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta_g)
theta_resid = matrix(rstan::summary(fit, pars = c("theta_resid"))$summary[,1], ncol = 4, byrow = TRUE)
tg = rstan::summary(fit, pars = c("theta"))$summary[,1]
lam = rstan::summary(fit, pars = c("lambda_identify"))$summary[,1]

rstan::summary(fit, pars = c("sig_sq_thetag_reg"))

cor(tg, sim_data$theta_g)

plot(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])

plot(theta_resid[,1], sim_data$theta[,1])
plot(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
plot(theta_resid[,2], sim_data$theta[,2])


plot(tg*lam[4] + theta_resid[,4], sim_data$theta_g*sim_data$lambda[4] + sim_data$theta[,4])


hist(tg*lam[1] + theta_resid[,1])

# plot(theta_resid[,1], sim_data$theta[,1])
cor(tg, sim_data$theta_g)
plot(tg, sim_data$theta_g)
plot(theta_resid[,1], sim_data$theta[,1])
cor(theta_resid[,1], sim_data$theta[,1])
cor(theta_resid[,2], sim_data$theta[,2])
cor(theta_resid[,3], sim_data$theta[,3])
cor(theta_resid[,4], sim_data$theta[,4])

df = data.frame(x1 = tg*lam[1] + theta_resid[,1], 
                x2 = sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
ggplot(df) + 
  aes(x = x1, y = x2) + 
  geom_point() + 
  xlim(-1.5, 1.5)

cor(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
cor(tg*lam[2] + theta_resid[,2], sim_data$theta_g*sim_data$lambda[2] + sim_data$theta[,2])
cor(tg*lam[3] + theta_resid[,3], sim_data$theta_g*sim_data$lambda[3] + sim_data$theta[,3])
cor(tg*lam[4] + theta_resid[,4], sim_data$theta_g*sim_data$lambda[4] + sim_data$theta[,4])

sim_data$lambda
rstan::summary(fit, pars = c("lambda"))$summary[,1]

rstan::summary(fit, pars = c("lambda_identify"))$summary[,1]

library(rstan)
traceplot(fit, pars = c("lambda_identify"))

traceplot(fit, pars = c("theta_resid[1,1]", "theta_resid[2,1]"))

rstan::summary(fit, pars = c("sig_sq_thetag_reg"))

# fit = inirt::inirt(data = fit_data$data, model = fit_data$model, predictors = fit_data$predictors, dims = fit_data$dims, 
#     h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, iter = fit_data$iter, 
#     init = fit_data$init)


# summary(fit, pars = "alpha")$summary[,"mean"]
aest = matrix(rstan::summary(fit, pars = "alpha")$summary[,1], nrow = 4, byrow = TRUE)
aest
plot(exp(aest[aest != 0.0]), exp(sim_data$alpha[sim_data$alpha != 0.0]))
cor(aest[1,1:5], sim_data$alpha[1,1:5])
plot(exp(aest[1,1:10]), exp(sim_data$alpha[1,1:10]))
cor(aest[2,21:40], alpha[2,21:40])
cor(aest[3,41:60], alpha[3,41:60])

dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 20, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])
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
