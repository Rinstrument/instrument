
start_val = function(x, type, isZ, n, nT, nD, J, nB) {
  
  lambda = x$lambda
  theta = x$theta
  delta = x$delta
  alpha = x$alpha
  beta = x$beta
  
  start = list(lambda = c(), theta = c(), delta = c(), 
               alpha = c(), beta = c())
  
  if(is.null(lambda)) {
    start$lambda = runif(nT - 1, -1, 1)
  } else {
    start$lambda = lambda
  }
  
  if(is.null(theta)) {
    start$theta = rep(rnorm(n, 0, 1), nT)
  } else {
    start$theta = theta
  }
  
  if(is.null(delta)) {
    start$delta = rep(rnorm(nD*J, -1, 1))
  } else {
    start$delta = delta
  }
  
  if(isZ == TRUE) {
    if(is.null(beta)) {
      start$beta = rnorm(nB, 0, 1)
    } else {
      start$beta = beta
    }
  } else {
    start$beta = c()
  }
  
  if(type == "1p") {
    start$alpha = c()
  } else if(type == "2p") {
    if(is.null(alpha)) {
      start$alpha = runif(J, 0.1, 1)
    } else {
      start$alpha = alpha
    }
  }
  
  if(type == "1p") {
    if(isZ == TRUE) {
      return(c(start$lambda, start$theta, start$delta, start$beta))
    } else {
      return(c(start$lambda, start$theta, start$delta))
    }
  } else {
    if(isZ == TRUE) {
      return(c(start$lambda, start$theta, start$delta, start$alpha, start$beta))
    } else {
      return(c(start$lambda, start$theta, start$delta, start$alpha))
    }
  }
  
}