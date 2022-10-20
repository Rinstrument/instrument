# --------------------------------------------------------------
# With no delta parameters

# library(ggplot2)
# # library(mirt)
library(devtools)
# library(truncnorm)
rm(list = ls())
load_all()
# set.seed(1234565322)
n = 1000
j = 20
data = matrix(0, nrow = n, ncol = j)
theta = seq(-1.5, 1.5, length.out = n) #rnorm(n, 0 , 1)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(theta[i]))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:j)
# cor(rowMeans(data), theta)
x = sis_theta_model(data, n = 10000)
x = matrix(as.vector(x), nrow = n, byrow = TRUE)
post = apply(x, 1, sum)
# post
# lm = sapply(as.data.frame(t(data)), function(x) {glm(data.frame(x)[,1] ~ 1, family = binomial())$coef}) 
# lm
# cor(post, lm)
# View(cbind(post, theta))
# rowMeans(data)
# theta
cor(post, theta)
cor(post, rowMeans(data))
# cor(lm, rowMeans(data))
# cor(lm, theta)

# hist(p_means)
# # cor(p_means, c(alpha, delta, theta))
# mean(p_means)
# glm(data[,1] ~ 1, family = binomial())$coef






# --------------------------------------------------------------
# With delta parameters
library(devtools)
library(mirt)
rm(list = ls())
load_all()
# set.seed(1234565322)
n = 200
j = 20
p = n + j
data = matrix(0, nrow = n, ncol = j)
theta = seq(-1.5, 1.5, length.out = n) #rnorm(n, 0 , 1)
delta = rnorm(j, 0, 1)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(theta[i] - delta[jj]))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:j)
mirt.coefs = coef(mirt(data, itemtype = "Rasch"), simplify = TRUE)$items
mirt.delta = mirt.coefs[, "d"]
cor(mirt.delta, delta)
# cor(rowMeans(data), theta)
x = sis_theta_model2(data, n = 5000, tol = 2.0e-30)
x = matrix(as.vector(x), nrow = p, byrow = TRUE)
post = apply(x, 1, sum)
length(post)
post[201:220]
cor(post[201:220], delta)
View(cbind(post[201:220], delta))
cor(post[1:200], theta)
View(cbind(post[1:200], theta))
cor(post[201:220], mirt.delta)












# --------------------------------------------------------------
# With delta & alpha parameters
library(devtools)
library(mirt)
library(truncnorm)
rm(list = ls())
load_all()
# set.seed(1234565322)
n = 200
j = 20
p = n + 2*j
data = matrix(0, nrow = n, ncol = j)
theta = seq(-1.5, 1.5, length.out = n) #rnorm(n, 0 , 1)
alpha = rtruncnorm(j, a = 0, b = 10, mean = 1.5)
delta = rnorm(j, 0, 1)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-((alpha[jj]*theta[i]) - delta[jj]))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:j)
apply(data, 2, mean)
mirt.coefs = coef(mirt(data, itemtype = "2PL"), simplify = TRUE)$items
mirt.delta = mirt.coefs[, "d"]
mirt.alpha = mirt.coefs[, "a1"]
cor(mirt.delta, delta)
cor(mirt.alpha, alpha)
# cor(rowMeans(data), theta)
x = sis4(data, n = 2000, tol = 20)
x = matrix(as.vector(x), nrow = p, byrow = TRUE)
post = apply(x, 1, sum)
length(post)
# post[201:220]
cor(post[201:220], delta)
# plot(post[201:220], delta)
# abline(a = 0, b = 1)
# View(cbind(post[201:220], delta))
cor(post[1:200], theta)
# plot(post[1:200], theta)
# abline(a = 0, b = 1)
# View(cbind(post[1:200], theta))
# cor(post[201:220], mirt.delta)
cor(post[221:240], alpha)
# View(cbind(post[221:240], alpha))
# cor(mirt.alpha, alpha)
# plot(post[221:240], alpha)
# abline(a = 0, b = 1)












# --------------------------------------------------------------
# With delta, alpha, second-order theta, lambda (loading) parameters
library(devtools)
# library(mirt)
library(truncnorm)
rm(list = ls())
load_all()
# set.seed(1234565322)
n = 1000
d = 5
j = 20
p = (d + 1)*n + 2*j*d + d
data = matrix(0, nrow = n, ncol = d*j)
lambda = c(0.7, -0.8, -0.9, 1.0, 0.6) #/ sqrt(1.5)
thetag = seq(-4, 4, length.out = n) #rnorm(n, 0, 3) #seq(-1.5, 1.5, length.out = n) #rnorm(n, 0 , 1)
theta = matrix(0, nrow = n, ncol = d)
for(i in 1:d) {
  theta[,i] = rnorm(n, lambda[i]*thetag, sd(thetag)/10)
  # theta[,i] = theta[,i] / sqrt(1.5)
}
alpha = rtruncnorm(d*j, a = 0, b = 7, mean = 1.5)
delta = rnorm(d*j, 0, 1)
jj_d = rep(1:d, each = j)
for(i in 1:n) {
  for(jj in 1:(d*j)) {
    prb = (1 / (1 + exp(-((alpha[jj]*theta[i, jj_d[jj]]) - delta[jj]))))
    # prb = (1 / (1 + exp(-((theta[i, jj_d[jj]])))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:(d*j))
apply(data, 2, mean)
d_index = rep(1:j, d)
d_start = which(d_index == 1)
d_end = which(d_index == j)
x = sis5(data, n_dimensions = d, dimension_start = d_start, dimension_end = d_end, n = 2000, tol = 400)
x = matrix(as.vector(x), nrow = p, byrow = TRUE)
post = apply(x, 1, sum)
length(post)
# is general (second-order theta correct)?
cor(post[1:n], thetag)
plot(post[1:n], thetag)
abline(a = 0, b = 1)
# st_tg = post[1:400]/(sd(post[1:400]))
# st_est = thetag/sd(thetag)
# summary(lm(st_tg ~ st_est))
# plot(post[1:400], thetag)
# abline(a = 0, b = 1)
post[c((n+1):(n+d))]
cor(post[(n+d+1):(2*n+d)], theta[,1])
# plot(post[(n+d+1):(2*n+d)], theta[,1])
# abline(a = 0, b = 1)
cor(post[(n+d+1):(2*n+d)], theta[,1])
cor(post[(2*n+d+1):(3*n+d)], theta[,2])
cor(post[(3*n+d+1):(4*n+d)], theta[,3])
cor(post[(4*n+d+1):(5*n+d)], theta[,4])
cor(post[(5*n+d+1):(6*n+d)], theta[,5])
# cor((post[n+1]*post[1:n]) + post[(n+d+1):(2*n+d)], theta[,1])
# cor((post[402]*post[1:400]) + post[804:1203], theta[,2])
# cor((post[403]*post[1:400]) + post[1204:1603], theta[,3])
# plot((post[403]*post[1:400]) + post[1204:1603], theta[,3])
# abline(a = 0, b = 1)
cor(post[((d + 1)*n + d + 1):((d + 1)*n + d + j*d)], alpha)
plot(post[((d + 1)*n + d + 1):((d + 1)*n + d + j*d)], alpha)
abline(a = 0, b = 1)
# plot(post[((d + 1)*n + d):((d + 1)*n + d + j*d - 1)], alpha)

# cor(post[((d + 1)*n + d + j*d):((d + 1)*n + d + 2*j*d - 1)], delta)
cor(post[((d + 1)*n + d + d*j + 1):((d + 1)*n + d + 2*j*d)], delta)
plot(post[((d + 1)*n + d + d*j + 1):((d + 1)*n + d + 2*j*d)], delta)
abline(a = 0, b = 1)




















# Reversing order (from the second-order solution).
# We have the second-order IRT model
# Now let's make sis5 work for any special case,
# starting with the simple univariate IRT model
library(devtools)
library(mirt)
library(truncnorm)
rm(list = ls())
load_all()
# set.seed(1234565322)
n = 200
j = 20
p = n + 2*j
data = matrix(0, nrow = n, ncol = j)
theta = seq(-3, 3, length.out = n) #rnorm(n, 0 , 1)
alpha = rtruncnorm(j, a = 0, b = 10, mean = 1.5)
delta = rnorm(j, 0, 1)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-((alpha[jj]*theta[i]) - delta[jj]))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:j)
apply(data, 2, mean)
mirt.coefs = coef(mirt(data, itemtype = "2PL"), simplify = TRUE)$items
mirt.delta = mirt.coefs[, "d"]
mirt.alpha = mirt.coefs[, "a1"]
cor(mirt.delta, delta)
cor(mirt.alpha, alpha)
# cor(rowMeans(data), theta)
x = sis5(data, n_dimensions = 1, dimension_start = 1, dimension_end = ncol(data), 
         n_second_order = 0, n = 2000, tol = 0.2*2000)
x = matrix(as.vector(x), nrow = p, byrow = TRUE)
post = apply(x, 1, sum)
length(post)
cor(post[1:n], theta)
cor(post[(n+1):(n+j)], alpha)
cor(post[(n+1+j):(n+2*j)], delta)












# Multidimensional IRT model
library(devtools)
library(mirt)
library(truncnorm)
rm(list = ls())
load_all()
n = 1000
d = 2
j = 20
p = d*n + j*d*d + j*d
data = matrix(0, nrow = n, ncol = d*j)
theta = matrix(rnorm(d*n, 0, 1), nrow = n, byrow = FALSE)
alpha = rtruncnorm(d*j*d, a = 0, b = 10, mean = 1.5)
for(dd in 1:d) {
  for(jj in 1:(d*j)) {
    if(jj < dd) {
      alpha[(dd-1)*d*j + jj] = 0
    }
  }
}
alpha_as_mat = matrix(alpha, ncol = j*d, byrow = TRUE)
delta = rnorm(d*j, 0, 1)
for(i in 1:n) {
  for(jj in 1:(d*j)) {
    prb = (1 / (1 + exp(-(sum(alpha_as_mat[,jj]*theta[i,]) )))) #  # - delta[jj]
    # prb = (1 / (1 + exp(-((theta[i, jj_d[jj]])))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:(d*j))
apply(data, 2, mean)
mirt.mod = mirt(data, 2, method = "MHRM")
mod.pars = coef(mirt.mod, simplify = TRUE)$items
cor(alpha_as_mat[1,], -mod.pars[,1])
cor(alpha_as_mat[2,], -mod.pars[,2])
# cor(delta, -mod.pars[,3])

# d_index = rep(1:j, d)
# d_start = which(d_index == 1)
# d_end = which(d_index == j)
# x = sis5(data, n_dimensions = 5, dimension_start = 0, dimension_end = 0, # take care of dimensions??
        #  n_second_order = 0, n = 20, tol = 2)
x = sis5(data, n_dimensions = d, dimension_start = 0, dimension_end = 0, # take care of dimensions??
         n_second_order = 0, n = 5000, tol = 0.2*5000)
x = matrix(as.vector(x), nrow = p, byrow = TRUE)
post = apply(x, 1, sum)
length(post)
cor(post[1:n], theta[,1])
cor(post[(n+1):(2*n)], theta[,2])
cor(post[(2*n+1):(3*n)], theta[,3])
cor(post[(3*n+1):(4*n)], theta[,4])
cor(post[(4*n+1):(5*n)], theta[,5])
# post[d*n + j*d + 1]
post_alpha_as_mat = matrix(post[(d*n+1):(d*n+j*d*d)], ncol = j*d, byrow = TRUE)
post_alpha_as_mat = matrix(post[(d*n+1):(d*n+j*d*d)], nrow = d, byrow = TRUE)

# post_alpha_as_mat = matrix(post[(d*n+1):(d*n+j*d*d)], ncol = j*d, byrow = FALSE)

# alpha_as_mat[1,]
cor(post_alpha_as_mat[1,], alpha_as_mat[1,])
# plot(post_alpha_as_mat[1,], alpha_as_mat[1,])
cor(post_alpha_as_mat[2,-d*j], alpha_as_mat[2,-d*j])
cor(post_alpha_as_mat[3,], alpha_as_mat[3,])
cor(post_alpha_as_mat[4,], alpha_as_mat[4,])
cor(post_alpha_as_mat[5,], alpha_as_mat[5,])

cor(post[(d*n+j*d*d+1):(d*n+j*d*d+j*d)], delta)


cor(alpha_as_mat[1,], mod.pars[,1])
cor(alpha_as_mat[2,], mod.pars[,2])
cor(post_alpha_as_mat[1,], -mod.pars[,1])
cor(post_alpha_as_mat[2,], mod.pars[,2])
