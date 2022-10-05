library(ggplot2)
set.seed(1234565322)
n = 500
data = matrix(rep(0, n), n, 1)
j = 1
beta = 2
sd = 0.7
for(i in 1:n) {
  data[i] = rnorm(1, beta, sd)
}
summary(lm(data ~ 1))
miss = sample(1:n, size = floor(n / 20), replace = FALSE)
true_miss = data[miss]
data[miss] = NA
n_miss = length(miss)
start = c(rep(0, 1), rep(1, 1), rep(0, n_miss))
p = length(start)
iterations = 6e5
burn = 5e5
greedy_iterations = 4e5
accept = rep(0, iterations)
iter_save = (iterations - burn)
iter_save = 1.2 * iter_save
x = matrix(data = 0, p, iter_save)
indices = rep(0, p)
validation_indexes = c(1)
validation_lower = c(0.001)
validation_upper = c(Inf)
lp_lm(start, data, n_miss)
amc(x = x, x_start = start, iter = iterations, burn = burn, greedy_iterations = greedy_iterations, 
    a = 0.234, data = data, lp_select = 0, accept = accept, validation_indexes = validation_indexes, 
    validation_lower = validation_lower, validation_upper = validation_upper, 
    gam_correct_iter_post_burn = indices, p_reg = n_miss)
get_draws = function(x, param, indices) {
  y = x[param, 1:indices[param]]
  return(y)
}
param_id = 5
drws = get_draws(x, param_id, indices)
pdat = data.frame(it = 1:length(drws), y = drws)
ggplot(pdat) + aes(it, y) + geom_line()
post_mean = mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
post_mean
#cor(post_mean[-(1:2)], true_miss)
