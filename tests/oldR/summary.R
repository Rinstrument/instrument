#' @exportS3Method summary hltObj
summary.hltObj = function(object, ...) {
  
  args = list(...)
  param = args$param
  dimension = args$dimension
  nT = object$nT
  
  if("digits" %in% names(args)) {
    digits = args$digits
  } else {
    digits = 3
  }
  
  if("transpose" %in% names(args)) {
    transpose = args$transpose
  } else {
    transpose = TRUE
  }
  
  post = object$post
  nms = colnames(post)
  
  if (param == "all") {
    smry = apply(post, 2, smy, digits = digits)
  } else if (param == "lambda") {
    smry = apply(post[, grepl("lambda", nms), drop = FALSE], 2, smy, digits = digits)
    
    #lambda = smry[1,]
    #nlambda = ncol(smry)
    #sdy = numeric(nlambda)
    # for(i in 1:nlambda) {
    #   sdy[i] = sd(summary.hltObj(object, param = "theta", dimension = i)[, 1])
    # }
    # sdx = sd(summary.hltObj(object, param = "theta", dimension = nlambda + 1)[, 1])
    # lambda_std = (lambda * sdx) / sdy
    
    # cor_mat = matrix(data = NA, nrow = nrow(object$theta), ncol = nT)
    # for(i in 1:nT) {
    #   cor_mat[, i] = summary.hltObj(object, param = "theta", dimension = i)[, 1]
    # }
    # smry = rbind(smry, std.mean = round(cor(cor_mat)[1:(nT - 1), nT], digits = digits))
    
    #smry = rbind(smry, std.mean = round(lambda_std, digits = digits))
  } else if (param == "alpha") {
    smry = apply(post[, grepl("^[a]", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "delta") {
    smry = apply(post[, grepl("^[d]", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "beta") {
    smry = apply(post[, grepl("beta", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "theta") {
    
    if("dimension" %in% names(args)) {
      dimension = args$dimension
      nT = object$nT
    } else {
      warning("Since no dimension argument was specified, summaries are returned
              for the general latent dimension.")
      dimension = nT
    }
    
    total_theta = nrow(object$theta)
    n_per_theta = total_theta / nT
    n_per_theta * dimension
    smry = t(object$theta[((n_per_theta * (dimension - 1)) + 1):(n_per_theta * dimension), ])
  } else if (param == "correlation") {
    nT = object$nT
    corr = matrix(0, nrow = nrow(object$theta), ncol = nT)
    
    for(i in 1:nT) {
      corr[, i] = summary(object, param = "theta", dimension = i)
    }
    
    smry = round(cor(corr), digits = digits)
    colnames(smry) = rownames(smry) = paste0("theta", 1:nT)
  }
  
  if(transpose == FALSE) {
    return(smry)
  } else {
    return(t(smry))
  }
}

smy = function(x, digits) {
  round(c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.975))),
        digits = digits)
}
