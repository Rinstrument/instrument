# retrieve all columns of the same name from the output from each iteration
# of the simulation
retrieve_by_name = function(res, name) {

	n = length(res)

	out = vector(mode = "list", length = n)

	for(i in 1:n) {
    out[[i]] = res[[i]][, ..name]
	}

	out = data.table::data.table(do.call(cbind, out))

	return(out)

}

# Evaluate simulation output using the standard metrics like mse, bias, etc...
evaluate_model = function(res) {

	n_sim = length(res) # length of simulation

	# Parameter names
	pnames = res[[1]][, 'parameter']

	# true values are equal across the res[[i]]'s since they were sampled once
	truth = res[[1]][, 'true']

	# posterior estimate table for all parameters
	estimates = retrieve_by_name(res, 'mean')

	# distinct names for all mean parameters
	data.table::setnames(estimates, paste0('mean', 1:n_sim))

	# CI lower bounds for each iteration
	ci_lower = retrieve_by_name(res, 'quantile_0.025')
	data.table::setnames(ci_lower, paste0('quantile_0.025_', 1:n_sim))

	# CI upper bound estimates for each iteration
	ci_upper = retrieve_by_name(res, 'quantile_0.975')
	data.table::setnames(ci_upper, paste0('quantile_0.975_', 1:n_sim))

	# merge CI data
	ci = cbind(ci_lower, ci_upper, truth)

	# Merge the truth
	estimates = cbind(estimates, truth)

	# Average result over simulation iterations
	estimates[ , `:=`(avg_sim = rowMeans(.SD)), .SDcols = grep('mean', colnames(estimates))]

	# Bias
  estimates[ , `:=`(bias = avg_sim - true)]

	# Bias Monte Carlo SE
	estimates[, `:=`(bias_mcse = sqrt((1 / (n_sim*(n_sim-1)))*rowSums((.SD - avg_sim)^2))), .SDcols = grep('mean', colnames(estimates))]

	# MSE
	estimates[, `:=`(mse = (1/n_sim)*rowSums((.SD - true)^2)), .SDcols = grep('mean', colnames(estimates))]

	# MSE Monte Carlo SE
	estimates[, `:=`(mse_mcse = sqrt(rowSums((((.SD - true)^2) - mse)^2) / (n_sim*(n_sim-1)))), .SDcols = grep('mean', colnames(estimates))]

	# Empirical SE
	estimates[, `:=`(emp_se = sqrt((1/(n_sim-1))*rowSums((.SD - avg_sim)^2))), .SDcols = grep('mean', colnames(estimates))]

	# Empirical SE Monte Carlo SE
	estimates[, `:=`(emp_se_mcse = emp_se/sqrt(2*(n_sim-1)))]

	# Confidence interval coverage
	for(i in 1:n_sim) {
		lower = paste0('quantile_0.025_', i)
		upper = paste0('quantile_0.975_', i)
		name = paste0('cover_', i)
		ci[, (name) := 1*(..lower <= true & true <= ..upper)]
	}

	# Estimated coverage
	ci[, `:=`(avg_cov = rowMeans(.SD)), .SDcols = grep('cover', colnames(ci))]

	# Estimated coverage Monte Carlo SE
	ci[, `:=`(avg_cov_mcse = sqrt((avg_cov*(1-avg_cov))/n_sim))]

	# return results to the user
	out = cbind(pnames, estimates, ci)

	return(out)

}
