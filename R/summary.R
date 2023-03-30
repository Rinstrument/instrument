#' Summarize a theta2 object
#' 
#' Produce summary data table for theta2 output.
#' 
#' @param object theta2Obj model object
#' @param pars "default" will give the default parameters of the model.
#' Other options will give only the specified vector of parameter names.
#' For example, pars = c("theta", "alpha").
#' @param probs probability vector for quantiles of posterior estiamtes
#' @param ... not used
#' 
#' @return a `data.table` of draws
#' 
#' @importFrom rstan summary
#' @export
#' 
summary.theta2Obj = function(object, pars = 'default', probs = c(0.025, 0.50, 0.975), 
  ...) {
  
  stanfit = object$stanfit

  # column names
  all_names = stanfit@sim$fnames_oi

  all_par_names = stanfit@model_pars

  # remove unused parameter names
  if(pars == 'default') {
    all_par_names = setdiff(all_par_names, c('alpha_l', 'eta3pl_l', 'delta_l', 
      'eta3pl', 'db', 'ab', 'xb', 'nu', 'c', 'lp__'))
  } else {
    all_par_names = pars
  }

  # subset parameter names to keep
  all_names = all_names[grep(paste0(all_par_names, collapse = '|'), all_names)]

  # remove _l parameters (unformatted)
  all_names = all_names[!grepl('_l', all_names)]
  
  # extract draws and coerce to data.table (fast)
  draws = data.table::as.data.table(rstan::extract(stanfit, pars = all_par_names, 
    permute = FALSE))

  # parameters pivoted longer
  draws = data.table::dcast(draws, iterations + chains ~ parameters, 
    value.var = 'value')

  par_names = colnames(draws)[-c(1:2)]

  # find where the parameter names are to re-order them
  match_order = match(all_names, par_names)

  # alpha parameters need to be transformed with exp(a), since they were estimated
  # with exp(a) * (theta - d)
  alpha_par_names = par_names[grep('alpha', par_names)]

  draws = draws[
      , (alpha_par_names) := lapply(.SD, \(x) { exp(x) }), .SDcols = alpha_par_names
    ]

  # posterior summary function
  summary_theta2Obj = function(x, probs) { 
      c(mean(x), sd(x), quantile(x, probs = probs)) 
    }

  # apply posterior summary to all parameters
  draws = draws[
      , lapply(.SD[, -c(1, 2)], summary_theta2Obj, probs = probs)
    ]

  # summary names
  draws = draws[
      , summary := c('mean', 'sd', paste0('quantile_', probs))
    ]

  # transpose data to long format (column are summary statistics)
  draws = data.table::transpose(draws, keep.names = 'parameter', 
    make.names = 'summary')

  # reorder result for convenience
  draws = draws[
      , match := order(match_order)
    ]
  
  data.table::setorder(draws, match)

  draws[
      , -match
    ]

  return(draws)

}

# write a function that takes a pattern and gives the parameter summaries for filtering
# using my data.table summary function... getPar("theta", dim = 1), getPar("theta")
# just translates into the right filter function