
evaluate_model = function(res, n_cores) {

	# number of parallel apply operations
	n_res = length(res)

	# first element of summary
	est_sum_nres = res[[1]][, 'mean']

	# 2 to n_res elements fo summary
	for(i in 2:n_res) {
		est_sum_nres = est_sum_nres + res[[i]][, "mean"] 
	}

	# mean of model estimates
	est_mean_nres = est_sum_nres / n_res

	

}