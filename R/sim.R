#' Simulate the HLT model
#' 
#' @param n number of observations
#' @param ntheta number first-level of latent dimensions
#' @param lambda latent factor coefficients
#' @param id number of questions
#' @param dL number of levels for each question
#' @param nB number of regression parameters. nB = ncol(z).
#' @param beta what value to set the regression parameters.
#' @param type type of model to simulate. `type = "1p"` for the partial
#' credit model. `type = "2p"` for the generalized partial
#' credit model. 
#' 
#' @return a `list` containing 
#' \itemize{
#' \item{x}{ - matrix of simulated item responses}
#' \item{theta}{ - matrix of true latent ability}
#' \item{id}{ - I.Ds for item membership to each dimension}
#' \item{namesx}{ - vector of column names for each item}
#' \item{s.cor}{ - true correlations between latent ability dimensions}
#' \item{s.delta}{ - true difficulty parameters}
#' \item{s.lambda}{ - true loading parameters}
#' \item{s.alpha}{ - true discrimination parameters}
#' }
#' 
#' @importFrom truncnorm rtruncnorm
#' @export
#' 
#' @examples 
# # PCM, No regression, 2 dimensions, binary response
# xdat = hltsim(type = "1p", n = 50, ntheta = 2, lambda = c(0.5, 0.5),
#               id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), dL = 2)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = xdat$id, iter = 5e5, burn = 3e5, type = "1p",
#           delta = 0.1)
# 
# # PCM, No regression, 3 dimensions, binary response
# xdat = hltsim(type = "1p", n = 200, ntheta = 3, lambda = c(0.2, 0.2, 0.1),
#               id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2), dL = 2)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = xdat$id, iter = 5e5, burn = 3e5, type = "1p",
#           delta = 0.04)
# 
# # PCM, No regression, 2 dimensions
# xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7),
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = NULL)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#           iter = 200, burn = 100, delta = 0.013)
# 
# # PCM, No regression, 4 dimensions
# xdat = hltsim(type = "1p", n = 50, ntheta = 2, lambda = c(0.5, 0.5),
#               id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), dL = 2)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = xdat$id, iter = 5e5, burn = 3e5, delta = 0.03)
# 
# 
# # GPCM, No regression, 2 dimensions
# xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7),
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = NULL)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#           iter = 200, burn = 100, delta = 0.013)

#' xdat = hltsim(n = 250, type = "2p", ntheta = 4, 
#'               lambda = c(0.5, 0.8, 0.9, 0.4), id = c(rep(0, 15),         
#'               rep(1, 15), rep(2, 15), rep(3, 15)), dL = 2)
#' mod1 = hlt(x = xdat$x, id = xdat$id, iter = 12e1, 
#'            burn = 6e1, delta = 0.023)
#' 
#' xdat = hltsim(n = 250, type = "2p", ntheta = 4, 
#'               lambda = c(0.5, 0.8, 0.9, 0.4), id = c(rep(0, 15),         
#'               rep(1, 15), rep(2, 15), rep(3, 15)), dL = 2,
#'               beta = c(0.5, -0.7), nB = 2)
#' mod2 = hlt(x = xdat$x, id = xdat$id, z = xdat$z, 
#'            iter = 12e1, burn = 6e1, delta = 0.023, nchains = 1)
#' 
#' 
hltsim = function(type, n, ntheta, lambda, id, dL, nB, beta = NULL) {
  
  if(!is.null(beta)) {
    regression = TRUE
  } else {
    regression = FALSE
  }
  
  nT = ntheta + 1
  theta = matrix(0, n, nT)
  theta[, nT] = rnorm(n, 0, sqrt(1.5))
  s.beta = beta
  
  if(regression == TRUE) {
    z = matrix(sample(c(0, 1), size = n * nB, replace = TRUE), nrow = n, 
               ncol = nB)
    theta[, nT] = theta[, nT] + z %*% s.beta + rnorm(n, 0, sqrt(1.5))
  }
  
  for(i in 1:(nT - 1)) {
    theta[, i] = lambda[i] * theta[, nT] + rnorm(n, 0, sqrt(1.5))
  }
  
  s.cor = cor(theta)
  J = length(id)
  nD = dL - 1
  s.theta = theta
  s.alpha = rtruncnorm(J, a = 0, mean = 1, sd = 0.5)
  
  s.delta = mapply(1:J, FUN = function(x) {sort(rnorm(dL, 1, 0.5))}, 
                   SIMPLIFY = "matrix")
  s.delta[1, ] = rep(0, J)
  s.lambda = lambda
  
  x = matrix(0, n, J)
  for (i in 1:n) {
    
    for (j in 1:J) {
      
      if(type == "2p") {
        exp_part = exp(cumsum((s.alpha[j] * (s.theta[i, id[j] + 1])) - s.delta[, j]))
      } else if(type == "1p") {
        exp_part = exp(cumsum((s.theta[i, id[j] + 1]) - s.delta[, j]))
      }
      
      x[i, j] = sample(1:(dL) - 1, size = 1, prob = exp_part / sum(exp_part))
    }
    
  }
  
  if(regression == TRUE) {
    ret = list(x = x, 
                z = z,
                id = id,
                theta = theta,
                namesx = paste0("x", 1:J),
                s.beta = s.beta,
                s.cor = s.cor,
                s.delta = s.delta,
                s.lambda = s.lambda)
  } else {
    ret = list(x = x, 
                theta = theta,
                id = id,
                namesx = paste0("x", 1:J),
                s.cor = s.cor,
                s.delta = s.delta,
                s.lambda = s.lambda)
  }
  
  if(type == "1p") {
    return(ret)
  } else {
    ret = append(ret, list(s.alpha))
    names(ret)[length(ret)] = "s.alpha"
    return(ret)
  }
  
}
