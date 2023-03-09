
evaluate_model = function(res) {

	# number of parallel apply operations
	n_res = length(res)

	# first element of summary
	est_sum_nres = res[[1]][, 'mean']

	truth = res[[1]][, 'true']

	# 2 to n_res elements fo summary
	for(i in 2:n_res) {
		est_sum_nres = est_sum_nres + res[[i]][, "mean"] 
	}

	# mean of model estimates
	est_mean_nres = est_sum_nres / n_res

	bias = est_mean_nres - truth

	bias_mcse = sqrt((1 / (n_res*(n_res - 1))) * sum((est_mean_nres)^2))

	mse = mean((est_mean_nres - truth)^2)

	rmse = sqrt(mse)



}