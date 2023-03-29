#' Loadings from a thata2 model fit object
#' 
#' Get factor loadings (rotated and unrotated) from a theta2 model fit object.
#' 
#' @param x model fit object from theta2, i.e., x = theta2(data, ...)
#' @param rotate what rotation to use? Default is 'varimax'. Other option is
#' 'promax'. Only one can be specified at a time 
#' @param ... other arguments to pass down to the varimax and promax functions 
#' in the stats package. See ?varimax for more information.
#' 
#' @return a list containing the unrotated loadings, rotated loadings according
#' to the rotate argument, and the matrix of posterior estimates of the unrotated 
#' loading parameters. 
#' 
#' @example 
#' 
#' 
#'
#' 
#'
#' 
loadings = function(x, rotate = c('varimax', 'promax'), ...) {

  rotate = c('varimax', 'promax')

  # output
  out = vector('list', 3)
  names(out) = c('unrotated', 'rotated', 'posterior')
  
  # default to varimax
  rotate = rotate[1]

  # get summaries
  fit_smy = theta2::summary.theta2Obj(x)
  
  # subset loading estimates
  posterior = fit_smy[
      grep('alpha', parameter),
    ]

  # find dimensions of alpha
  posterior[['parameter']]
  
  ?stringr::str_sub()

  # order alpha into a matrix
  posterior

  # return stacked format as a separate matrix
  
  if(rotate == 'varimax') {
    lrot = varimax(l, ...)
  } else if(rotate == 'promax') {
    lrot = promax(l, ...)
  } else {
    stop('Invalid method selected. Select either varimax or promax.')
  }

  # unrotated factor loadings
  out['unrotated']

  # rotation according to rotate
  out['rotated']
  
  # posterior summaries for loadigns
  out['posterior'] = posterior

  return(out)

}