
#' @useDynLib lmHOIRT, .registration=TRUE
#' @importFrom Rcpp evalCpp
lmHOIRT = function(model, data, start, iterations, burn, greedy_iterations) {
  start = c(rep(1, j), rep(0, j), rep(0, n), rep(0, 1))
  iterations = 1000000
  burn = 350000
  greedy_iterations = 200000
  accept = rep(0, iterations)
  iter_save = (iterations - burn)
  iter_save = 1.2 * iter_save
  x = matrix(data = 0, p, iter_save)
  indices = rep(0, p)
  amc(x = x, x_start = start, iter = iterations, burn = burn, greedy_iterations = greedy_iterations, 
      a = 0.234, data = data, lp_select = 0, accept = accept, validation_indexes = 0:(j-1), 
      validation_lower = rep(0, j), gam_correct_iter_post_burn = indices, p_reg = 1)
}