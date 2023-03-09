
get_bias = function(sim_data, fit_smy) {

}

merge_estimate_truth = function(fit_smy, mirt_data) {

	fit_data = mirt_data$fit_data
	sim_data = mirt_data$sim_data

	true = sim_data$true

	# merge true values into fit_smy table
	fit_smy = fit_smy[true, on = 'parameter']

	# # Bias, 
	# bias = get_bias(sim_data, fit_smy)

	# # empirical SE, 
	# emp_se = get_emp_se(sim_data, fit_smy)

	# # mean-squared error, 
	# mse = get_mse(sim_data, fit_smy)
	
	# # coverage
  # coverage = get_coverage(sim_data, fit_smy)

	return(fit_smy)

}