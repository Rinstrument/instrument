# test amc algorithm with 2pl 2-level HOIRT model
library(devtools)
library(ggplot2)
rm(list = ls())
load_all()
set.seed(1234565322)
n = 500
n1l = 3
dims1l = c(12, 14, 15)
j = sum(dims1l)
p = 2*j + (n1l + 1)*n + n1l
#model_indices = list(0:(j-1), (j):(2*j-1), (2*j):(p-p_reg-1), (p-p_reg-1) + p_reg)
alpha = exp(runif(j, 0, 1.5)) - 0.2
delta = rnorm(j, 0, 1)
lambda = c(0.5, 0.6, 0.8)
#beta = 0.5
#z = sample(0:1, n, replace = TRUE)
thetag = rnorm(n, 0, sqrt(1.5))
thetas = matrix(data = 0, nrow = n, ncol = n1l + 1)
for(i in 1:n1l) {
  thetas[,i] = rnorm(n, thetag * lambda[i], sqrt(1.5))
}
thetas[,n1l + 1] = thetag
#theta = rnorm(n, z*beta, sqrt(1.5)) #seq(-3, 3, length.out = n)
#summary(lm(theta ~ z))
true = c(alpha, delta, as.vector(thetas), lambda)
data = matrix(0, nrow = n, ncol = j + 1 + n1l)
for(i in 1:n) {
  d_start = 1
  for(d in 1:n1l) {
    d_end = d_start + dims1l[d] - 1
    for(jj in d_start:d_end) { #(2 + n1l):(j + 1 + n1l)
      prb = (1 / (1 + exp(-(alpha[jj]*thetas[i,d] - delta[jj]))))
      data[i, jj  + 1 + n1l] = sample(c(1, 0), 1, prob = c(prb, 1 - prb))
    }
    d_start = d_start + dims1l[d]
  }
}
data[1,1] = n1l
data[1, 2:(n1l+1)] = dims1l
#colnames(data) = paste0("x", 1:j)
#data = cbind(data, "z1" = z)
start = c(rep(1, j), rep(0, j), rep(0, (n1l + 1)*n), rep(0, n1l))
# thetas_start_data = data[, -(1:(n1l+1))]
# thetas_start = matrix(data = 0, nrow = n, ncol = n1l + 1)
# thetas_start[,1] = apply(thetas_start_data[, 1:12], 1, \(x) {mean(x)})
# thetas_start[,1] = (thetas_start[,1] - mean(thetas_start[,1])) * 3
# thetas_start[,2] = apply(thetas_start_data[, 13:(13+14-1)], 1, \(x) {mean(x)})
# thetas_start[,2] = (thetas_start[,2] - mean(thetas_start[,2])) * 3
# thetas_start[,3] = apply(thetas_start_data[, (13+14):(13+14+15-1)], 1, \(x) {mean(x)})
# thetas_start[,3] = (thetas_start[,3] - mean(thetas_start[,3])) * 3
# thetas_start[,4] = rowMeans(thetas_start[,1:3])
# thetas_start[,4] = (thetas_start[,4] - mean(thetas_start[,4])) * 3
# cor(thetas_start[,1], thetas[,1])
# cor(thetas_start[,2], thetas[,2])
# cor(thetas_start[,3], thetas[,3])
# cor(thetas_start[,4], thetas[,4])
# start = c(rep(1, j), rep(0, j), as.vector(thetas_start), rep(0.5, 3))
iterations = 3.8e6
burn = 3.5e6
greedy_iterations = 1.5e6
accept = rep(0, iterations)
iter_save = (iterations - burn)
iter_save = 1.2 * iter_save
x = matrix(data = 0, p, iter_save)
indices = rep(0, p)
model_indices = list(0:(j-1), j:(2*j - 1), (2*j):(2*j + 1*n - 1), (2*j + 1*n):(2*j + 2*n - 1), (2*j + 2*n):(2*j + 3*n - 1), (2*j + 3*n):(2*j + 4*n - 1), (2*j + (n1l + 1)*n):(2*j + (n1l + 1)*n + n1l - 1))
validation_indexes = c(0:(j-1),j:(2*j - 1), (2*j + (n1l + 1)*n):(2*j + (n1l + 1)*n + n1l - 1))
validation_lower = c(rep(0, j), rep(-20, j), rep(-10, n1l))
validation_upper = c(rep(10, j), rep(20, j), rep(10, n1l))
lp_2pl_ho2l(start, data, 0)
lp_2pl_ho2l(true, data, 0)
#ml solution
# opt = optim(par = start, fn = llk_2pl_logit, method = "L-BFGS-B",
#     data = data, theta = theta_start)
# par = opt$par
amc(x = x, x_start = start, iter = iterations, burn = burn, greedy_iterations = greedy_iterations, 
    a = 0.234, a_greedy = 1.0, data = data, lp_select = 0, accept = accept, validation_indexes = validation_indexes, 
    validation_lower = validation_lower, validation_upper = validation_upper, 
    gam_correct_iter_post_burn = indices, p_reg = 1)
get_draws = function(x, param, indices) {
  y = x[param, 1:indices[param]]
  return(y)
}
param_id = 4
drws = get_draws(x, param_id, indices)
pdat = data.frame(it = 1:length(drws), y = drws)
ggplot(pdat) + aes(it, y) + geom_line()
# mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# cor(
#   c(alpha, delta, theta),
#   mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
# )
#head(x[48,])
post_mean = mapply(\(row, ind_max) {mean(x[row, 1:ind_max])}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])
post_tail = mapply(\(row, ind_max) {tail(x[row, 1:ind_max], 1)}, row = 1:nrow(x), ind_max = indices[1:nrow(x)])

cor(
  alpha,
  post_mean[1:41]
)
View(cbind(alpha, post_mean[1:41]))
cor(
  delta,
  post_mean[42:(42 + 41 - 1)]
)
cor(thetas[,1], post_mean[model_indices[[3]] + 1])
cor(thetas[,2], post_mean[model_indices[[4]] + 1])
cor(thetas[,3], post_mean[model_indices[[5]] + 1])
cor(thetas[,4], post_mean[model_indices[[6]] + 1])

cor(post_mean,
    true)
# cor(
#   theta,
#   post_mean[83:()]
# )
# post_mean[model_indices[[4]] + 1]

hist(post_mean[model_indices[[3]] + 1])
post_mean[model_indices[[7]] + 1]

lp_2pl_ho2l(start, data, 1)
lp_2pl_ho2l(true, data, 1)
lp_2pl_ho2l(post_mean, data, 1)
lp_2pl_ho2l(post_tail, data, 1)
lp_2pl_ho2l(c(post_mean[-(2083:2085)], c(0.5, 0.5, 0.5)), data, 1)
