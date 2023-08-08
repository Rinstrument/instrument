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
#' @examples
#' 
#' library(theta2)
#' load('tests/case_study_results/fit_exploratory_3dim.RData')
#' 
#' ld = loadings(fit, 'varimax')
#' ld['unrotated']
#' ld['rotated']
#' ld['posterior']
#' 
#' ld = loadings(fit, 'promax')
#' ld['unrotated']
#' ld['rotated']
#' ld['posterior']
#' 
#' ld = loadings(fit, 'oblimin')
#' ld['unrotated']
#' ld['rotated']
#' ld['posterior']
#' 
#'
#' 
loadings = function(x, rotate = c('varimax', 'promax', 'oblimin'), ...) {

  # output
  out = vector('list', 3)

  names(out) = c('unrotated', 'rotated', 'posterior')
  
  # default to varimax
  rotate = rotate[1]

  # get summaries
  fit_smy = instrument::summary.instrumentObj(x)
  
  # subset loading estimates
  posterior = fit_smy[
      grep('alpha', parameter),
    ]

  # find dimensions of alpha
  par_names = posterior[['parameter']]

  # dimensions of alpha
  a_dims = as.double(
      stringr::str_extract_all(par_names[length(par_names)], '(\\d)+', 
        simplify = TRUE)
    )
  
  # unrotated factor loadings (posterior means)
  unrotated = t(
      matrix(posterior[['mean']], a_dims[1], a_dims[2], byrow = FALSE)
    )

  # select rotation
  if(rotate == 'varimax') {
    
    rotated = GPArotation::Varimax(unrotated, ...)

  } else if(rotate == 'promax') {

    rotated = promax(unrotated, ...)

  } else if(rotate == 'oblimin') {

    rotated = GPArotation::oblimin(unrotated, ...)

  } else {

    stop('Invalid method selected. Select either varimax, promax, or oblimin.')

  }

  # unrotated factor loadings
  out['unrotated'] = unrotated

  # rotation according to rotate
  out['rotated'] = rotated
  
  # posterior summaries for loadigns
  out['posterior'] = posterior

  return(out)

}