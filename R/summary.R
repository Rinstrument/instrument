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
summary.instrumentObj = function(object, pars = 'default', probs = c(0.025, 0.50, 0.975), 
  ...) {
  
  stanfit = object$stanfit

  # # column names
  # all_names = stanfit@sim$fnames_oi

  # all_par_names = stanfit@model_pars

  # better column names
  all_par_names = 
    names(object$stanfit@par_dims)[
      sapply(object$stanfit@par_dims,
             \(x) {
               !any(x == 0) 
             })
      ]

  # remove unneeded parameters
  all_par_names = all_par_names[!grepl('delta_l|alpha_l|beta_l|zeta_l\\[|db|ab|xb|nu|c|eta3pl|_elong|lp__', all_par_names)]

  name_sorter = function(x) {

    # x = draws$parameter

    # x = c('a[2,2]', 'b[1,2]','a[2,1]', 'b[1,1]','a[1,2]', 'b[2,1]')

    y = str_split(x, '\\[', simplify = TRUE)

    name = y[,1]

    index = y[,2]

    ppm = matrix(as.numeric(str_remove_all(str_split(index, ',', simplify = TRUE), pattern = '\\[|\\]')), ncol = 2)

    sort_order = cbind.data.frame(
      raw = x,
      name = name,
      ind1 = ppm[,1],
      ind2 = ppm[,2]
    )

    x = sort_order[with(sort_order, order(name, ind1, ind2)), ]$raw

    return(x)

  }

  # # remove unused parameter names
  # if(pars == 'default') {
  #   all_par_names = setdiff(all_par_names, c('alpha_l', 'eta3pl_l', 'delta_l', 
  #     'eta3pl', 'db', 'ab', 'xb', 'nu', 'c', 'zeta_l', 'lp__')) # this was zeta_l[ ?? Why?

  #   # # More parameters to remove here
  #   # grepl('zeta_', c('zeta_1'))

  #   # grepl(c('zeta_l_sd_elong_'), all_par_names)

  #   # all_par_names = all_par_names[str_detect(pattern = c('zeta_l_sd_elong'), string = all_par_names, negate = TRUE)]



  # } else {
  #   all_par_names = pars
  # }

  # # subset parameter names to keep
  # all_names = all_names[grep(paste0(all_par_names, collapse = '|'), all_names)]

  # # remove _l parameters (unformatted)
  # all_names = all_names[!grepl('_l', all_names)]
  
  # extract draws and coerce to data.table (fast)
  draws = data.table::as.data.table(rstan::extract(stanfit, pars = all_par_names, 
    permute = FALSE))

  # parameters pivoted longer
  draws = data.table::dcast(draws, iterations + chains ~ parameters, 
    value.var = 'value')

  par_names = colnames(draws)[-c(1:2)]

  # find where the parameter names are to re-order them
  match_order = match(name_sorter(par_names), par_names)

  # alpha parameters need to be transformed with exp(a), since they were estimated
  # with exp(a) * (theta - d)
  alpha_par_names = par_names[grep('alpha', par_names)]

  draws = draws[
      , (alpha_par_names) := lapply(.SD, \(x) { exp(x) }), .SDcols = alpha_par_names
    ]

  # posterior summary function
  summary_instrumentObj = function(x, probs) { 
      c(mean(x), sd(x), quantile(x, probs = probs)) 
    }

  # apply posterior summary to all parameters
  draws = draws[
      , lapply(.SD[, -c(1, 2)], summary_instrumentObj, probs = probs)
    ]

  # summary names
  draws = draws[
      , summary := c('mean', 'sd', paste0('quantile_', probs))
    ]

  # transpose data to long format (column are summary statistics)
  draws = data.table::transpose(draws, keep.names = 'parameter', 
    make.names = 'summary')

  # let's inspect draws here and see if we can't take care of parameter names
  # at this stage

  # This gives a unique list of parameters, good job.
  # unique(str_split(draws[, parameter], '\\[', simplify = TRUE)[,1])

  # reorder result for convenience
  draws = 
    draws[
        , match := order(match_order)
      ]
  
  data.table::setorder(draws, match)

  draws[
      , match := NULL
    ]

  return(draws)

}

# write a function that takes a pattern and gives the parameter summaries for filtering
# using my data.table summary function... getPar("theta", dim = 1), getPar("theta")
# just translates into the right filter function