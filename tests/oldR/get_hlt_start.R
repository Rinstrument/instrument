
#' hlt Starting Values
#' 
#' Get starting values from hlt fit object
#' 
#' @param x hlt model fit object
#' @param nchains number of chains to get starting values
#' 
#' @return a list of lists with starting values for each chain
#' 
#' @importFrom utils tail
#' @export
get_hlt_start = function(x, nchains = 1) {
  
  if(nchains > 1) {
    start_list = vector(mode = "list", length = nchains)
    
    for(i in 1:nchains) {
      start_list[[i]] = get_hlt_start(x, nchains = 1)
    }
    
    return(start = start_list)
  }
  
  post = tail(x$post, 1)
  isRegress = !is.null(x$args$z)
  is1p = x$args$type == "1p"
  nitem = ncol(x$args$x)
  
  if(isRegress) {
    nZ = ncol(x$args$z)
  }
  
  nT = x$nT
  
  if(isRegress) {
    
    if(is1p) {
      start = list(lambda = post[1, paste0("lambda", 1:(nT - 1))],
                   theta = x$theta[,1],
                   delta = post[1, grepl("^[d]", colnames(post))],
                   alpha = c(),
                   beta = post[1, paste0("beta", 1:nZ)])
    } else {
      start = list(lambda = post[1, paste0("lambda", 1:(nT - 1))],
                   theta = x$theta[,1],
                   delta = post[1, grepl("^[d]", colnames(post))],
                   alpha = post[1, paste0("a", 1:nitem)],
                   beta = post[1, paste0("beta", 1:nZ)])
    }
    
  } else {
    
    if(is1p) {
      start = list(lambda = post[1, paste0("lambda", 1:(nT - 1))],
                   theta = x$theta[,1],
                   delta = post[1, grepl("^[d]", colnames(post))],
                   alpha = c(),
                   beta = c())
    } else {
      start = list(lambda = post[1, paste0("lambda", 1:(nT - 1))],
                   theta = x$theta[,1],
                   delta = post[1, grepl("^[d]", colnames(post))],
                   alpha = post[1, paste0("a", 1:nitem)],
                   beta = c())
    }
  }
  
  return(start)
}
