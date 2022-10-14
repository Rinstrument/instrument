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
x = sis3(data, n = 2000, tol = 20)
x = matrix(as.vector(x), nrow = p, byrow = TRUE)
post = apply(x, 1, sum)
length(post)
# post[201:220]
cor(post[201:220], delta)
plot(post[201:220], delta)
abline(a = 0, b = 1)
# View(cbind(post[201:220], delta))
cor(post[1:200], theta)
plot(post[1:200], theta)
abline(a = 0, b = 1)
# View(cbind(post[1:200], theta))
# cor(post[201:220], mirt.delta)
cor(post[221:240], alpha)
# View(cbind(post[221:240], alpha))
# cor(mirt.alpha, alpha)
plot(post[221:240], alpha)
abline(a = 0, b = 1)
