
make_gt_tables = function(mod_assessment) {

  # Subset column names relevant to summary
  cnames = c('parameter', 'mse', 'mse_mcse', 'emp_se', 'emp_se_mcse', 'bias',
    'bias_mcse', 'avg_cov', 'avg_cov_mcse')

  # subset cols
  data = mod_assessment[, ..cnames]

  # Remove rows with Inf values
  data = data[!is.nan(rowSums(data[, -1])), ]

  # begin with alpha parameters
  data[grep('alpha', parameter), ] %>% 
    gt::gt() %>% 
    gt::gtsave(filename = 'tests/simulation/results/alpha_recovery.docx')

  # delta parameters
  data[grep('delta', parameter), ] %>% 
    gt::gt() %>% 
    gt::gtsave(filename = 'tests/simulation/results/delta_recovery.docx')

  # theta parameters
  data[grep('theta', parameter), ] %>% 
    gt::gt() %>% 
    gt::gtsave(filename = 'tests/simulation/results/theta_recovery.docx')

}