#' Summarize a theta2 object
#' 
#' Produce summary matrix for theta2 output.
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
summary.theta2Obj = function(object, pars = "default", probs = c(0.025, 0.50, 0.975), 
  ...) {
  
  stanfit = object$stanfit

  str(stanfit)
  all_names = stanfit@sim$fnames_oi

  
  all_names = all_names[grep(paste0(all_par_names, collapse = '|'), all_names)]

  all_par_names = names(stanfit@par_dims)

  if(pars == "default") {
    all_par_names = setdiff(all_par_names, c("alpha_l", "eta3pl_l", "delta_l", 
      "eta3pl", "db", "ab", "xb", "nu", "c", "lp__"))
  } else {
    all_par_names = pars
  }

  # w = data.table(a = letters[1:10], b = 1:10)
  # match(w$a, letters[10:1])

  # w[a == order(letters[1:10], rev = TRUE), ]
  # setorder()
  
  draws = as.data.table(rstan::extract(stanfit, pars = all_par_names, permute = FALSE))

  ordered_par_names = unique(draws[['parameters']])

  draws = data.table::as.data.table(draws)

  draws = data.table::dcast(draws, iterations + chains ~ parameters, 
    value.var = "value")

  par_names = colnames(draws)[-c(1:2)]

  match_order = match(ordered_par_names, par_names)

  # draws = draws[
  #   , match_order := ..match_order
  # ]

  alpha_par_names = par_names[grep('alpha', par_names)]

  draws = draws[
      , (alpha_par_names) := lapply(.SD, \(x) { exp(x) }), .SDcols = alpha_par_names
    ]

  summary_theta2Obj = function(x, probs) { 
      c(mean(x), sd(x), quantile(x, probs = probs)) 
    }

  draws = draws[, lapply(.SD[, -c(1, 2)], summary_theta2Obj, probs = probs)]

  draws = draws[, summary := c("mean", "sd", paste0("quantile_", probs))]

  draws = data.table::transpose(draws, keep.names = "parameter", make.names = "summary")

  draws = draws[
      , match := match_order
    ]

  setorder(draws, match)
  # draws = draws[order(parameter)

  return(draws)

}

# write a function that takes a pattern and gives the parameter summaries for filtering
# using my data.table summary function... getPar("theta", dim = 1), getPar("theta")
# just translates into the right filter function