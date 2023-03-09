
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
summary.theta2Obj = function(object, pars = "default", probs = c(0.025, 0.50, 0.975), ...) {

  stanfit = object$stanfit

  all_par_names = names(stanfit@par_dims)

  if(pars == "default") {
    all_par_names = setdiff(all_par_names, c("alpha_l", "eta3pl_l", "delta_trans", 
      "eta3pl", "db", "ab", "xb", "nu", "c", "lp__"))
  } else {
    all_par_names = pars
  }
  
  draws = data.table::as.data.table(rstan::extract(stanfit, pars = all_par_names, permute = FALSE))

  draws = data.table::dcast(draws, iterations + chains ~ parameters, 
    value.var = "value")

  summary_theta2Obj = function(x, probs) {c(mean(x), quantile(x, probs = probs))}

  draws = draws[, lapply(.SD[, -c(1, 2)], summary_theta2Obj, probs = probs)]

  return(draws)

}