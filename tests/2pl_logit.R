library(ggplot2)
library(mirt)
rm(list = ls())
load_all()
set.seed(1234565322)
n = 500
j = 20
p_reg = 1
p = 2 * j + n + p_reg
model_indices = list(0:(j-1), (j):(2*j-1), (2*j):(p-p_reg-1), (p-p_reg-1) + p_reg)
alpha = exp(runif(j, 0, 1.5)) - 0.2
delta = rnorm(j, 0, 1)
beta = 0.5
z = sample(0:1, n, replace = TRUE)
theta = rnorm(n, z*beta, sqrt(1.5)) #seq(-3, 3, length.out = n)
summary(lm(theta ~ z))
true = c(alpha, delta, theta, beta)
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {
    prb = (1 / (1 + exp(-(alpha[jj]*theta[i] - delta[jj]))))
    data[i, jj] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
  }
}
colnames(data) = paste0("x", 1:j)
data = cbind(data, "z1" = z)
start = c(rep(1, j), rep(0, j), rep(0, n), rep(0, 1))
iterations = 4e5
burn = 3e5
greedy_iterations = 2e5
accept = rep(0, iterations)
iter_save = (iterations - burn)
iter_save = 1.2 * iter_save
x = matrix(data = 0, p, iter_save)
indices = rep(0, p)
validation_indexes = 0:(j-1)
validation_lower = rep(0, j)
validation_upper = rep(Inf, j)
amc(x = x, x_start = start, iter = iterations, burn = burn, greedy_iterations = greedy_iterations, 
    a = 0.234, a_greedy = 1.0, data = data, lp_select = 0, accept = accept, validation_indexes = validation_indexes, 
    validation_lower = validation_lower, validation_upper = validation_upper, gam_correct_iter_post_burn = indices, p_reg = 1)
get_draws = function(x, param, indices) {
  y = x[param, 1:indices[param]]
  return(y)
}
param_id = 541
drws = get_draws(x, param_id, indices)
pdat = data.frame(it = 1:length(drws), y = drws)
ggplot(pdat) + aes(it, y) + geom_line()
# mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# cor(
#   c(alpha, delta, theta),
#   mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# )

post_mean = mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
cor(
  alpha,
  post_mean[model_indices[[1]] + 1]
)
cor(
  delta,
  post_mean[model_indices[[2]] + 1]
)
cor(
  theta,
  post_mean[model_indices[[3]] + 1]
)
View(cbind(
  theta,
  post_mean[model_indices[[3]] + 1]
))
post_mean[model_indices[[4]] + 1]

lp_2pl_logit_reg(start, data, 1)
lp_2pl_logit_reg(post_mean, data, 1)
lp_2pl_logit_reg(true, data, 1)




# c(beta0, sd)
# mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# cor(
#   c(beta0, beta, sd),
#   mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# )
# post = mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# lp_lm(post, data)
# lp_lm(c(post[-length(post)], 1), data)
