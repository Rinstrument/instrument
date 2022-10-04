library(ggplot2)
library(mirt)
set.seed(1234565322)
n = 500
j = 20
p_reg = 1
p = 2 * j + n + p_reg
model_indices = list(0:(j-1), (j):(2*j-1), (2*j):(p-1), p_reg)
alpha = exp(runif(j, 0, 1.5)) - 0.2
delta = rnorm(j, 0, 1)
beta = 0.4
z = sample(0:1, n, replace = TRUE)
theta = rnorm(n, z*beta, sqrt(1.5)) #seq(-3, 3, length.out = n)
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
iterations = 1000000
burn = 800000
greedy_iterations = 500000
accept = rep(0, iterations)
iter_save = (iterations - burn)
iter_save = 1.2 * iter_save
x = matrix(data = 0, p, iter_save)
logPostPtr = postFunc("lp2plr")
indices <- amc(x = x, x_start = start, iter = iterations, burn = burn, greedy_iterations = greedy_iterations, 
               a = 0.234, data = data, logPostPtr = logPostPtr, accept = accept,
               validation_indexes = 0:(j-1), validation_lower = rep(0, j), p_reg = 1)
get_draws = function(x, param, indices) {
  y = x[param, 1:indices[param, ]]
  return(y)
}
param_id = 5
drws = get_draws(x, param_id, indices)
pdat = data.frame(it = 1:length(drws), y = drws)
ggplot(pdat) + aes(it, y) + geom_line()
mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
cor(
  c(alpha, delta, theta),
  mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
)

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

c(beta0, sd)
mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
cor(
  c(beta0, beta, sd),
  mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
)
post = mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
lp_lm(post, data)
lp_lm(c(post[-length(post)], 1), data)
