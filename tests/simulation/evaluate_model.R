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

	# true values are equal across the res[[i]]'s since they were sampled once
	truth = res[[1]][, 'true']

	# posterior estimate table for all parameters
	estimates = retrieve_by_name(res, 'mean')

	# distinct names for all mean parameters
	data.table::setnames(estimates, paste0('mean', 1:n_sim))

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

	estimates[, `:=`(avg_sim = rowMeans(.SD)), .SDcols = grep('mean', colnames(estimates))]

	estimates[, lapply(.SD, \(x, y) {x - y}, y = avg_sim), .SDcols = grep('mean', colnames(estimates))]

	estimates[, `:=`(new = .SD - 2), .SDcols = grep('mean', colnames(estimates))]

	estimates[ , `:=`(avg_sim = rowMeans(.SD)), .SDcols = grep('mean', colnames(estimates))][
						 , `:=`(bias = avg_sim - true)][
						 , `:=`(sq_dev_mean = rowSums((.SD))), .SDcols = grep('mean', colnames(estimates))
						 ][
						 , `:=`(sq_dev_mean = rowSums((.SD - avg_sim)^2)), .SDcols = grep('mean', colnames(estimates))][
						 , `:=`(bias_mcse = sqrt((1 / (n_sim*(n_sim-1)))*rowSums(.SD))), .SDcols = grep('mean', colnames(estimates))]

	# mean of model estimates
	est_mean_nres = est_sum_nres / n_res

	bias = est_mean_nres - truth

	bias_mcse = sqrt((1 / (n_res*(n_res - 1))) * sum((est_mean_nres)^2))

	mse = mean((est_mean_nres - truth)^2)

	rmse = sqrt(mse)



}