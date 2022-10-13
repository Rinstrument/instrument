library(ggplot2)
# library(mirt)
library(devtools)
library(truncnorm)
rm(list = ls())
load_all()
set.seed(1234565322)
n = 400
j = 20
p = 2*j + n
model_indices = list(0:(j-1), (j):(2*j-1), (2*j):(p-1))
alpha = rtruncnorm(j, a = 0, mean = 1, sd = 0.5) #exp(runif(j, 0, 1.5)) - 0.2
delta = rnorm(j, 1, 0.5) #rnorm(j, 0, 1)
#theta = rnorm(n, 0.0, sqrt(1.5))
theta = rnorm(n, 0, sqrt(1.5))
true = c(alpha, delta, theta)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(alpha[jj]*theta[i] - delta[jj]))))
    # prb = (1 / (1 + exp(-(alpha[jj]*(theta[i] - delta[jj])))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:j)

m_iterations = 1e6
n_iterations = 200 #m_iterations / 10
p_means = sir_irt3(data, m_iterations = m_iterations, n_iterations = n_iterations)
# cor(p_means, c(alpha, delta, theta))
hist(p_means[,1])
hist(exp(p_means[,1]/100 - logadd(p_means[,1]/100)))
hist(exp(p_means[,1]))

any(is.na(p_means))
post = apply(p_means, 1, mean)
hist(p_means[100,])
cor(post[1:20], alpha)
cor(post[21:40], delta)
cor(post[41:440], theta)
