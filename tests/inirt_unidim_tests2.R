# INIRT with univariate theta dimension
devtools::install(dependencies = FALSE)

# Test 1: simplest possible settings
n = 800
ncat = 4
j = 25
d = 1
k = 0
uk = 0
ncategi = c(rep(ncat, j))
ncateg_max = max(ncategi)
alpha = matrix(0, d, j)
a_design = as.matrix(data.frame(x1 = rep(1, n)))
b_alpha = 0.8
for(dd in 1:d) {
  alpha[dd, ] = sort(runif(j, 0.2, 1.5))
}
delta = matrix(nrow = j, ncol = ncateg_max - 1)
d_design = as.matrix(data.frame(x1 = rep(1, n)))
b_delta = 1.3
for(jj in 1:j) {
  delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
}
delta = cbind(0, delta)
theta = matrix(0, nrow = n, ncol = d)
for(dd in 1:d) {
  theta[, dd] = rnorm(n, 0, 1)
}
beta = NULL
predictors = NULL
start_index = 1
beta_dstart = NULL
beta_dend = NULL
data = matrix(0, nrow = n, ncol = j)
for(i in 1:n) {
  for(jj in 1:j) {#                                                                                  z[i, ] %*% zeta
    prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(theta[i, ])) - (delta[jj, ] + b_delta*d_design[i,])))))
    prb[1] = 1.0
    prb = c(prb, 0)
    prb = prb[-length(prb)] - prb[2:length(prb)]
    data[i, jj] = sample(1:ncategi[[jj]], 1, prob = prb)
  }
}
colnames(data) = c(paste0("x", 1:j))
dims = 1
item_id = 1:j
sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
fit_data = list(data = data, item_id = item_id, model = NULL, predictors = predictors,
    n_pranef_cor = 2,
    dims = dims, h2_dims = 0, h2_dim_id = NULL, structural_design = list(alpha = a_design, delta = d_design), 
    method = "vb", weights = NULL, tol_rel_obj = 0.0002, iter = 5e3, init = "random")
# data = data; item_id = item_id; model = NULL; predictors = predictors; predictors_ranef = NULL; ranef_id = NULL; 
# predictors_ranef_corr = NULL; n_pranef_cor = NULL;
# dims = dims; h2_dims = 0; h2_dim_id = NULL; structural_design = list(alpha = a_design, delta = d_design); 
# method = "vb"; weights = NULL; tol_rel_obj = 0.0002; iter = 5e3; init = "random";
rm(list = setdiff(ls(), c("fit_data", "sim_data")))
ls()
fit = inirt::inirt(data = fit_data$data, item_id = fit_data$item_id, model = fit_data$model, predictors = fit_data$predictors, 
    dims = fit_data$dims, h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
    method = fit_data$method, weights = fit_data$weights, tol_rel_obj = fit_data$tol_rel_obj, iter = fit_data$iter, init = fit_data$init)
