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

	# true values are equal across the res[[i]]'s since they were sampled once
	truth = res[[1]][, 'true']

	# posterior estimate table for all parameters
	estimates = retrieve_by_name(res, 'mean')

	

	# mean of model estimates
	est_mean_nres = est_sum_nres / n_res

	bias = est_mean_nres - truth

	bias_mcse = sqrt((1 / (n_res*(n_res - 1))) * sum((est_mean_nres)^2))

	mse = mean((est_mean_nres - truth)^2)

	rmse = sqrt(mse)



}