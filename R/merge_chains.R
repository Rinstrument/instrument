#' Merge Chains from hlt method
#'
#' @param x object of class "hltObjList"
#' @param ... other arguments
#' 
#' @return a list of class `hltObj`. This class constructs a single `hltObj`
#' from a list of model fits by merging the chains into one matrix of draws.
#' 
#' @export
merge_chains = function (x, ...) {
  UseMethod("merge_chains", x)
}

#' @exportS3Method merge_chains hltObjList
merge_chains.hltObjList = function(x, ...) {
    
  nchains = length(x)
  post = do.call(rbind, Map(f = function(y) {y$post}, x))
  nT = x[[1]]$nT
  
  thetas = Map(f = function(y) {y$theta}, x)
  means = as.data.frame(Map(f = function(y) {y$mean}, thetas))
  sds = as.data.frame(Map(f = function(y) {y$sd}, thetas))
  
  rmeans = rowMeans(means)
  rsts = rowMeans(sds) 
  theta = data.frame(mean = rmeans, sd = rsts)
  
  result = list(post = post, 
                theta = theta,
                nT = nT,
                nchains = nchains,
                merged = TRUE,
                args = list(x = x[[1]]$args$x, z = x[[1]]$args$z, id = x[[1]]$args$id, 
                            iter = x[[1]]$args$iter, burn = x[[1]]$args$burn, 
                            delta = x[[1]]$args$delta, type = x[[1]]$args$type, 
                            start = x[[1]]$args$start,
                            progress = x[[1]]$args$progress))
  
  class(result) = c("hltObj")
  return(result)
  
}
