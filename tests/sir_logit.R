library(ggplot2)
# library(mirt)
library(devtools)
library(truncnorm)
rm(list = ls())
load_all()
set.seed(1234565322)
n = 400
j = 1
p = 1
beta = 0.25
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  # for(jj in 1:j) {
    # prb = (1 / (1 + exp(-(alpha[jj]*theta[i] - delta[jj]))))
    prb = (1 / (1 + exp(-(beta))))
    data[i, 1] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  # }
}
colnames(data) = paste0("x", 1:j)

m_iterations = 5e5
n_iterations = 100 #m_iterations / 10
p_means = sir_logit(data, m_iterations = m_iterations, n_iterations = n_iterations)
hist(p_means)
# cor(p_means, c(alpha, delta, theta))
mean(p_means)
glm(data[,1] ~ 1, family = binomial())$coef
