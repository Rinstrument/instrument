
#' Summarize a theta2 object
#' 
#' @importFrom rstan summary
#' @export
#' 
summary.theta2Obj = function(object, pars = "default", probs = c(0.025, 0.50, 0.975), ...) {
  probs = c(0.025, 0.50, 0.975)
  object = fit
  stanfit = object$stanfit

  all_par_names = names(stanfit@par_dims)

  all_par_names = setdiff(all_par_names, c("alpha_l", "eta3pl_l", "delta_trans", 
    "eta3pl", "db", "ab", "xb", "nu", "c", "lp__"))
  
  draws = data.table::as.data.table(rstan::extract(stanfit, pars = all_par_names, permute = FALSE))

  draws = data.table::dcast(draws, iterations + chains ~ parameters, 
    value.var = "value")

  summary_theta2Obj = function(x, probs) {c(mean(x), quantile(x, probs = probs))}

  draws = draws[, lapply(.SD[, -c(1, 2)], summary_theta2Obj, probs = probs)]

  return(draws)

}