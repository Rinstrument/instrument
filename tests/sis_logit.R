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
x = sis_theta_model(data, n = 2000)
x = matrix(as.vector(x), nrow = n, byrow = TRUE)
post = apply(x, 1, sum)
# post
lm = sapply(as.data.frame(t(data)), function(x) {glm(data.frame(x)[,1] ~ 1, family = binomial())$coef}) 
# lm
# cor(post, lm)
View(cbind(post, theta))
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






