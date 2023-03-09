# this file contains incorrect likelihood calc in simulation!!!

# # library(ggplot2)
# # library(mirt)
# library(devtools)
# library(rstan)
# rm(list = ls())
# # pkgbuild::compile_dll(force = TRUE)
# # load_all()
# devtools::install()
# # set.seed(1234565322)
# n = 800
# j = 25
# ncateg = 4
# k = 6
# uk = 6
# d = 1
# # p = 2*j + n + k
# # alpha = exp(runif(j, 0, 1.5)) - 0.2
# alpha = matrix(0, d, j)
# for(dd in 1:d) {
#   alpha[dd, ] = sort(runif(j, 0.2, 1.5))
# }
# alpha[2, ] = sort(alpha[2, ], decreasing = TRUE)
# alpha[2, 1] = 0
# #delta = rnorm(j, 0, 1)
# delta = matrix(nrow = j, ncol = ncateg - 1)
# for(jj in 1:j) {
#   delta[jj, ] = sort(rnorm(ncateg - 1, 0, 1))
# }
# delta = cbind(0, delta)
# # theta = rnorm(n, 0.0, sqrt(1.5))
# # theta = rnorm(n, 0, 1)
# theta = matrix(0, nrow = n, ncol = d)
# for(dd in 1:d) {
#   theta[, dd] = rnorm(n, 0, 1)
# }
# beta = c(-0.3, 0.2, 1, -0.6, -0.4, 0.6)
# # predictors = list(c(1:3), c(1,2,4))
# predictors = list(c(1:3,4:6))
# start_index = 1
# beta_dstart = numeric(d)
# beta_dend = numeric(d)
# for(dd in 1:d) {
#   beta_dstart[dd] = start_index
#   beta_dend[dd] = start_index + length(predictors[[dd]]) - 1
#   start_index = start_index + length(predictors[[dd]])
# }
# beta_mat = matrix(0, nrow = k, ncol = d)
# index = 1
# for(dd in 1:d) {
#   for(i in beta_dstart[dd]:beta_dend[dd]) {
#     beta_mat[i, dd] = beta[index]
#     index = index + 1
#   }
# }
# x = matrix(data = runif(uk*n,-1,1), nrow = n, ncol = uk)
# # xrp = x[,c(1,2,3,1,2,4)]
# # true = c(alpha, delta, theta)
# data = matrix(0, nrow = n, ncol = j)
# for(i in 1:n) {
#   for(jj in 1:j) {
#     prb = (1 / (1 + exp(-(sum(alpha[, jj]*(theta[i, ] + x[i,] %*% beta_mat)) - delta[jj, ]))))
#     prb = prb / sum(prb)
#     data[i, jj] = sample(1:ncateg, 1, prob = prb)
#   }
# }
# data = cbind(data, x)
# # coef(mirt(data, 1))
# colnames(data) = c(paste0("x", 1:j), paste0("k", 1:uk))
# # fit = inirt(data, method = "hmc", iter = 2000, chains = 1)
# # rm(n, j, alpha, delta, theta, i, jj)
# for(dd in 1:d) {
#   predictors[[dd]] = predictors[[dd]] + j
# }
# fit = inirt::inirt(data, predictors = predictors, dims = d, method = "vb", tol_rel_obj = 0.0001, iter = 1e4, init = 0)
# fit = inirt::inirt(data, predictors = predictors, dims = d, method = "hmc", chains = 1, iter = 300, init = 0)


# fit = inirt::inirt(data, predictors = predictors, dims = 1, method = "vb", tol_rel_obj = 0.0001, iter = 1e4)
# fit = inirt::inirt(data, predictors = predictors, dims = 1, method = "hmc", chains = 1, iter = 300)
# # smry = rstan::summary(fit)
# # hist(smry[[1]][1:100,"mean"])
# # hist(smry[[1]][101:120,"mean"])
# # hist(smry[[1]][121:140,"mean"])
# # cor(smry[[1]][1:500,"mean"], theta)
# ttas = summary(fit, pars = c("theta"))[[1]][,1]
# ttas1 = ttas[c(T, F)]
# ttas2 = ttas[c(F, T)]
# plot(ttas, theta[,1])
# plot(ttas2, theta[,2])
# cor(summary(fit, pars = c("theta"))[[1]][,1], theta)
# tal = matrix(summary(fit, pars = c("alpha"))[[1]][,1], nrow = 2, byrow = TRUE)
# plot(tal[1,], alpha[1,])
# plot(tal[2,], alpha[2,])
# cor(summary(fit, pars = c("delta"))[[1]][,1], delta)
# btas = summary(fit, pars = c("beta"))[[1]][,1]
# btas = matrix(btas, nrow = k, byrow = TRUE)
