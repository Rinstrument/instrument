# simulation function for mirt data
# type = c("irt", "mirt", "soirt", "birt")
# "irt", 
# "mirt", 
# "soirt", 
# "birt"

sim_mirt_pars = function(type) {

	n = 1000
	d = 4
	j = 20*d
	
	ncat = 2
	ncategi = c(rep(ncat, j))
	ncateg_max = max(ncategi)
	
	k = 0
	uk = 0
	
	alpha = matrix(0, d, j)
	a_design = as.matrix(data.frame(x1 = rep(1, n)))
	
	b_alpha = 1
	alpha_dominant = list(1:20, 21:40, 41:60, 61:80)
	
	for(dd in 1:d) {
		alpha[dd, alpha_dominant[[dd]]] = sort(runif(length(alpha_dominant[[dd]]), 1.7, 3.0))
		alpha[dd, setdiff(unlist(alpha_dominant), alpha_dominant[[dd]])] = sort(runif(length(unlist(alpha_dominant[-dd])), 0.2, 1.0))
	}
	
	delta = matrix(nrow = j, ncol = ncateg_max - 1)
	d_design = as.matrix(data.frame(x1 = rep(1, n)))
	b_delta = 1
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

	return(list(n = n, d = d, j = j, ncat = ncat, ncategi = ncategi, 
		ncateg_max = ncateg_max, k = k, uk = uk, alpha = alpha, a_design = a_design, 
		b_alpha = b_alpha, alpha_dominant = alpha_dominant, delta = delta, 
		d_design = d_design, b_delta = b_delta, theta = theta, beta = beta, 
		predictors = predictors, start_index = start_index, beta_dstart = beta_dstart, 
		beta_dend = beta_dend))

}

sim_mirt_data = function(type, pars) {

	n = pars$n
	d = pars$d
	j = pars$j 
	ncat = pars$ncat
	ncategi = pars$ncategi
	ncateg_max = pars$ncateg_max
	k = pars$k
	uk = pars$uk
	alpha = pars$alpha
	a_design = pars$a_design
	b_alpha = pars$b_alpha
	alpha_dominant = pars$alpha_dominant
	delta = pars$delta
	d_design = pars$d_design
	b_delta = pars$b_delta
	theta = pars$theta
	beta = pars$beta
	predictors = pars$predictors
	start_index = pars$start_index
	beta_dstart = pars$beta_dstart
	beta_dend = pars$beta_dend

	# set.seed(seed_replication_level)

	data = matrix(0, nrow = n, ncol = j)
	for(i in 1:n) {
		for(jj in 1:j) { #                                                                      + x[i, ] %*% beta_mat
			prb = (1 / (1 + exp(-(sum((alpha[, jj] + b_alpha*a_design[i,])*(theta[i, ])) - (delta[jj, 1:ncategi[jj]] + b_delta*d_design[i,])))))
			prb[1] = 1.0
			prb = c(prb, 0)
			prb = prb[-length(prb)] - prb[2:length(prb)]
			data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
		}
	}
	
	apply(data, 2, table)
	
	remove_gaps = function(x) {
		ord = order(x); vec = sort(x)
		old = unique(vec); replace = 1:length(unique(vec))
		names(replace) = old; names(vec) = vec
		new = replace[names(vec)]; names(new) = NULL
		return(new[ord])
	}
	
	data = apply(data, 2, remove_gaps)
	# data = cbind(data, x)
	# remove_gaps = function(x) {
	#   ord = order(x); vec = sort(x)
	#   old = unique(vec); replace = 1:length(unique(vec))
	#   names(replace) = old; names(vec) = vec
	#   new = replace[names(vec)]; names(new) = NULL
	#   return(new[ord])
	# }
	# data = apply(data, 2, remove_gaps)
	# apply(data, 2, table)
	# data = cbind(data, x)
	
	colnames(data) = c(paste0("x", 1:j)) #, paste0("z", 1:k)
	# for(dd in 1:1) {
	#   predictors[[dd]] = predictors[[dd]] + j
	# }
	# dims = 3
	# h2_dims = 1
	# h2_dim_id = list(1:20, 21:40, 41:60)
	sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, beta = beta, theta = theta)
	fit_data = list(data = data)

	# structure the alpha's in long format for comparison with estiamtes
	alpha = data.table::data.table(true = as.vector(sim_data$alpha))
	a_names = c()
	for(jj in 1:j) {
		for(i in 1:d) {
			a_names = c(a_names, paste0("alpha[", i, ",", jj, "]"))
		}
	}
	
	alpha[, parameter := a_names]

	# structure the theta's in long format for comparison with estiamtes
	theta = data.table::data.table(true = as.vector(sim_data$theta))
	t_names = c()
	for(i in 1:n) {
		for(dd in 1:d) {
			t_names = c(t_names, paste0("theta[", i, ",", dd, "]"))
		}
	}
	
	theta[, parameter := t_names]

	# structure the delta's in long format for comparison with estiamtes
	delta = data.table::data.table(true = as.vector(sim_data$delta))
	d_names = c()
	for(ncat in ncateg_max:1) {
		for(jj in 1:j) {
			d_names = c(d_names, paste0("delta[", jj, ",", ncat, "]"))
		}
	}
	
	delta[, parameter := d_names]

	# bind true parameter values together
	sim_data$true = data.table::rbindlist(list(alpha, delta, theta))

	return(list(sim_data = sim_data, fit_data = fit_data))

}
