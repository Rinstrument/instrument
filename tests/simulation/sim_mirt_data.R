
# Simulate MIRT parameters given model type
# type - type of mirt model to simulate
# type = c("irt", "mirt", "soirt", "birt")
# "irt"
# "mirt" 
# "soirt" 
# "birt"
sim_mirt_pars = function(type) {

	# n = number of observations (sample size)
	n = 1000

	# d = number of dimensions (1st order)
	d = 4

	# j is number of total questions (20 q's per dimension)
	j = 20*d
	
	# number of response categories
	ncat = 5

	# number of categories per item
	ncategi = c(rep(ncat, j))

	# maximum number of categories
	ncateg_max = max(ncategi)
	
	# number of fixed effect regression pars on theta
	k = 0
	uk = 0

	# two common fixed effects shared across theta, delta, and alpha
	ff = as.matrix(data.frame(x1 = sample(c(0, 1), n, TRUE),
		x2.1 = sample(c(0, 1), n, TRUE), x2.2 = sample(c(0, 1), n, TRUE)))
	
	# alpha - discrimination parameter (slope from factor analysis)
	alpha = matrix(0, d, j)

	# predictor design matrix of the alpha's
	a_design = ff
	
	# effects for alpha regression pars (defined in ff)
	b_alpha = c(0.5, 1.0, 1.5)

	# define the dominant alphas for each dimension
	alpha_dominant = list(1:20, 21:40, 41:60, 61:80)
	
  # Dominant alphas follow Unif(1.7, 3.0)
	# non-dominant alphas follow Unif(0.2, 1.0)
	for(dd in 1:d) {
		alpha[dd, alpha_dominant[[dd]]] = sort(runif(length(alpha_dominant[[dd]]), 2, 5))
		alpha[dd, setdiff(unlist(alpha_dominant), alpha_dominant[[dd]])] = sort(runif(length(unlist(alpha_dominant[-dd])), 0.2, 1.5))
	}
	
	# delta - difficulty parameters
	delta = matrix(nrow = j, ncol = ncateg_max - 1)

	# fixed effect design for the deltas
	d_design = ff

	# effects on difficult due to the fixed effects of ff
	b_delta = c(-1.5, 1, -0.5)

	# deltas are ordered and follow N(0, 1)
	for(jj in 1:j) {
		delta[jj, 1:(ncategi[jj]-1)] = sort(rnorm(ncategi[jj] - 1, 0, 1))
	}

	# first level of deltas are zero for estimability
	delta = cbind(0, delta)
	
	# theta parameters - ability at individual level
	theta = matrix(0, nrow = n, ncol = d)
	
	# ability follows N(0, 1)
	for(dd in 1:d) {
		theta[, dd] = rnorm(n, 0, 1)
	}
	
	# fixed effects on the theta's
	beta = c(1.5, 1, -1)

	# 
	# start_index = 1
	# beta_dstart = NULL
	# beta_dend = NULL

	return(list(n = n, d = d, j = j, ncat = ncat, ncategi = ncategi, 
		ncateg_max = ncateg_max, k = k, uk = uk, alpha = alpha, a_design = a_design, 
		b_alpha = b_alpha, alpha_dominant = alpha_dominant, delta = delta, 
		d_design = d_design, b_delta = b_delta, theta = theta, beta = beta,
		ff = ff))

}

# Given the set of parameters we generated, sample a new data set
# pars = sim_mirt_pars()
sim_mirt_data = function(type, pars) {

	# pull values from sim_mirt_pars() output
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
	ff = pars$ff
	# start_index = pars$start_index
	# beta_dstart = pars$beta_dstart
	# beta_dend = pars$beta_dend

	# set.seed(seed_replication_level)

	# sample data set
	data = matrix(0, nrow = n, ncol = j)

	# for individual i and question j, calculate the probability for responding at
	# each level of the response
	for(i in 1:n) {
		for(jj in 1:j) {
			# alpha 
			ar = as.vector(b_alpha %*% a_design[i,])
			# delta 
			dr = as.vector(b_delta %*% d_design[i,])
			# theta
			tr = as.vector(beta %*% ff[i, ])
			# regression equation
			nu = sum((t(alpha[, jj] + ar) %*% (theta[i, ] + tr))/5) - (delta[jj, 1:ncategi[jj]] + dr)
			# inverse logit
			prb = (1 / (1 + exp(-(nu))))
			# first level is 1
			prb[1] = 1.0
			# last level is 0
			prb = c(prb, 0)
			# difference between the two levels
			prb = prb[-length(prb)] - prb[2:length(prb)]
			# sample given probability vector
			data[i, jj] = sample(1:ncategi[jj], 1, prob = prb)
		}
	}
	
	# preview the data
	apply(data, 2, table)
	
	# remove gaps in response options, e.g., no 1,2,4,5 -> 1,2,3,4
	remove_gaps = function(x) {
		ord = order(x); vec = sort(x)
		old = unique(vec); replace = 1:length(unique(vec))
		names(replace) = old; names(vec) = vec
		new = replace[names(vec)]; names(new) = NULL
		return(new[ord])
	}
	
	# apply the remove gaps function
	data = apply(data, 2, remove_gaps)

	# preview the data
	apply(data, 2, table)

	# bind sampled data and predictor variables
	data = cbind(data, ff)
	
	# assign variable names to data
	colnames(data) = c(paste0('x', 1:j), paste0('z', 1:ncol(ff))) #, paste0("z", 1:k)
	
	# store simulated data
	sim_data = list(alpha = alpha, b_alpha = b_alpha, delta = delta, b_delta = b_delta, 
		beta = beta, theta = theta)

	# store data for fitting the model
	fit_data = list(data = data)

	# structure the alpha's in long format for comparison with estiamtes
	alpha = data.table::data.table(true = as.vector(sim_data$alpha))

	# update alpha parameter names for later use
	a_names = c()
	for(jj in 1:j) {
		for(i in 1:d) {
			a_names = c(a_names, paste0("alpha[", i, ",", jj, "]"))
		}
	}
	
	# assign new names to long parameter data
	alpha[, parameter := a_names]

	# structure the theta's in long format for comparison with estiamtes
	theta = data.table::data.table(true = as.vector(sim_data$theta))

	# update alpha parameter names for later use
	t_names = c()
	for(i in 1:n) {
		for(dd in 1:d) {
			t_names = c(t_names, paste0("theta[", i, ",", dd, "]"))
		}
	}
	
	# assign new names to long parameter data
	theta[, parameter := t_names]

	# structure the delta's in long format for comparison with estiamtes
	delta = data.table::data.table(true = as.vector(sim_data$delta))

	# update alpha parameter names for later use
	d_names = c()
	for(ncat in ncateg_max:1) {
		for(jj in 1:j) {
			d_names = c(d_names, paste0("delta[", jj, ",", ncat, "]"))
		}
	}
	
	# assign new names to long parameter data
	delta[, parameter := d_names]

	# bind true parameter values together
	sim_data$true = data.table::rbindlist(list(alpha, delta, theta))

	# list of return values
	return(list(sim_data = sim_data, fit_data = fit_data))

}
