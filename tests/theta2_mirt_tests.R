# INIRT higher-order IRT tests
devtools::install(dependencies = FALSE)
library(rstan)
stanc(file = "./inst/stan/theta2_mirt.stan", verbose = TRUE)
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
alpha_dominant = list(1:20, 21:40, 41:60, 61:80)
for(dd in 1:d) {
  alpha[dd, alpha_dominant[[dd]]] = sort(runif(length(alpha_dominant[[dd]]), 1.7, 3.0))
  alpha[dd, setdiff(unlist(alpha_dominant), alpha_dominant[[dd]])] = sort(runif(length(unlist(alpha_dominant[-dd])), 0.2, 1.0))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 1
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
  for(jj in 1:j) { #                                                                      + x[i, ] %*% beta_mat
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(theta[i, ])) - (delta[jj, 1:ncategi[jj]] + b_delta*d_design[i,])))))
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
# data = cbind(data, x)
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
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
fit_data = list(data = data)
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()
# library(devtools)
# library(Rcpp)
# compileAttributes()
# load_all()
data = fit_data$data
model = "theta1 = c(1:80)
         theta2 = c(1:80)
         theta3 = c(1:80)
         theta4 = c(1:80)"

#---
# itype = "2pl"
# method = "hmc"
# iter = 500
# warmup = 300
# chains = 1
# fweights = NULL
# cores = 1
# exploratory = TRUE

fit = theta2::theta2(data = data, model = model, itype = "2pl", exploratory = TRUE, 
  method = "hmc", iter = 100, warmup = 50, chains = 1, cores = 1)
class(fit)
summary(fit)
print(fit)
traceplot(fit, param = "")









#array(c(0.0), dim = 1)
library(rstan)

sim_data$theta
sim_data$theta_g
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta_g)
plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta_g)
theta_resid = matrix(rstan::summary(fit, pars = c("theta_resid"))$summary[,1], ncol = 4, byrow = TRUE)
tg = rstan::summary(fit, pars = c("theta"))$summary[,1]
lam = rstan::summary(fit, pars = c("lambda_identify"))$summary[,1]

rstan::summary(fit, pars = c("sig_sq_thetag_reg"))

plot(tg, sim_data$theta_g)

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

plot(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
cor(tg*lam[2] + theta_resid[,2], sim_data$theta_g*sim_data$lambda[2] + sim_data$theta[,2])
cor(tg*lam[3] + theta_resid[,3], sim_data$theta_g*sim_data$lambda[3] + sim_data$theta[,3])
cor(tg*lam[4] + theta_resid[,4], sim_data$theta_g*sim_data$lambda[4] + sim_data$theta[,4])

sim_data$lambda
rstan::summary(fit, pars = c("lambda"))$summary[,1]


rstan::summary(fit, pars = c("sig_thetag_reg"))$summary[,1]

rstan::summary(fit, pars = c("lambda_identify"))$summary[,1]

rstan::summary(fit, pars = c("beta_l"))$summary[,1]

library(rstan)
traceplot(fit, pars = c("lambda_identify"))
traceplot(fit, pars = c("sig_thetag_reg"))
traceplot(fit, pars = c("theta_resid[1,1]", "theta_resid[2,1]"))

rstan::summary(fit, pars = c("sig_thetag_reg"))

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
