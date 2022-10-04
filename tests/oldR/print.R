#' @exportS3Method print hltObj
print.hltObj = function(x, ...) {
  
  args = list(...)
  param = args$param
  dimension = args$dimension
  
  cat('Higher-order item response theory model with ', x$nT - 1, ' first-order domains', '\n', sep = '')
  cat(' Total number of parameters: ', ncol(x$post) + nrow(x$theta), '\n', sep = '')
  
  
  
  if(!('merged' %in% names(x))) {
    cat(' Model type (1p = "Partial Credit model"; 2p = "Generalized Partial Credit Model"): ', x$args$type, '\n', sep = '')
    cat(' iterations: ', x$args$iter, '; ', 'burn-in: ', x$args$burn, '\n', sep = '')
    cat(' Proposal standard deviation: ', x$args$delta, '\n', sep = '')
    cat(' Acceptance rate: ', round(x$accept.rate, digits = 3), '\n', sep = '')
    cat('\n')
  } else {
    cat('Merged run of ', x$nchains, ' parallel chains \n', sep = '')
    cat('\n')
  }
  
  cat('Loadings: ', summary(x, param = "lambda", digits = 3)[,1], '\n', sep = ' ')
  
  if(!is.null(x$args$z)) {
    cat('Latent regression beta estimates: ', summary(x, param = "beta", '\n', digits = 3)[,1], sep = ' ')
  }
  
}
