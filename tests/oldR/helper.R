# Calculate indices for sampling
#
# Helper to calculate indices for sampling
calc.ix = function(post_names, npar) {
  
  tot = 1
  strs = substr(post_names, 1, 1)
  ixx = c(1)
  
  for(i in 2:npar) {
    
    if(strs[i-1] != strs[i]){
      ixx[i] = tot + 1
      tot = tot + 1
    }
    
    ixx[i] = tot
  }
  
  ix = c()
  for (i in unique(ixx)) {
    ix = c(ix, which(ixx == i)[1])
  }
  
  ixe = c(ix[-1] - 1, length(ixx))
  
  return(list(ix = ix, ixe = ixe))
}

# Calculate simulated correlations
#
# Helper to calculate simulated correlations
cor.theta = function(theta) {
  theta = theta[, -ncol(theta)]
  cor(theta)
}