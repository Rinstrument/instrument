#' Explanatory and Descriptive Higher-Order Item Response Theory 
#' (Latent Trait Theory)
#' 
#' Fit a higher-order item response theory model under the generalized 
#' partial credit measurement model. The goal is to explain multiple latent dimensions
#' by a single higher-order dimension. We extend this model with an option to perform regression 
#' on the general latent dimension.
#'
#' @param x matrix of item responses. Responses must be integers where the 
#' lowest value is 0 and the highest value is the maximum possible response
#' for the item with no gaps. If a question is asked with 5 possible responses,
#' then the possible values should be c(0,1,2,3,4). For binary items, use
#' c(0,1).
#' @param z centered numeric matrix of predictors for the latent regression. 
#' Default is `z = NULL` so that no regression is performed. All columns
#' of this matrix must be numeric. For binary items, use the values c(0,1).
#' For continuous items, center the values on the mean and divide by the 
#' standard deviation (normalized). For factors with more than two levels, 
#' recode into multiple columns of c(0,1).
#' @param id I.D. vector indexing first-order latent dimension membership
#' for each of the first-order latent dimensions. We index starting from zero,
#' not one. If there are three first-order .
#' latent dimensions with 5 questions per dimension, then the vector will look
#' like c(0,0,0,0,0,1,1,1,1,1,2,2,2,2,2).
#' @param iter number of total iterations.
#' @param burn number of burn in iterations.
#' @param delta tuning parameter for Metropolis-Hanstings algorithm. Alter 
#' delta until acceptance.ratio =~ 0.234.
#' @param type type of Partial Credit Model to fit. If the partial credit model
#' is desired (i.e. all alpha parameters = 1), then choose `type = "1p"`.
#' If the Generalized Parial Credit Model is desired, then choose
#' `type = "2p"`. The default is `type = "2p"`.
#' @param start starting values for the Metropolis-Hastings algorithm. 
#' @param nchains number of independent MCMC chains. Default is `nchains = 1`.
#' @param verbose print verbose messages. Defaults to `FALSE`.
#' Provide a `list` with the following named arguments:
#' `list(lambda = c(), theta = c(), delta = c(), alpha = c(), beta = c())`
#' \itemize{
#'   \item{lambda - }{vector of starting values for the latent factor loadings.}
#'   \item{theta - }{vector of starting values for the abilities.}
#'   \item{delta - }{vector of starting values for the difficulties.}
#'   \item{alpha - }{vector of starting values for the slope parameters.}
#'   \item{beta - }{vector of starting values for the latent regression parameters}
#' }
#' If you choose specify starting values, then the lengths of the starting value
#' vectors must match the number of parameters in the model.
#' 
#' @param progress boolean, show progress bar? Defaults to TRUE.
#'
#'
#' @return A `list` containing:
#' 
#' \itemize{
#' \item{post}{ - A `matrix` of posterior estimates. Rows are the draws and columns 
#' are the named parameters.}
#' \item{accept.rate}{ - acceptance rate of MH algorithm}
#' \item{theta}{ - `matrix` of mean (first column) and standard deviation 
#' (second column) of posterior estimates of ability parameters}
#' \item{nT}{ - number of latent factors estimated}
#' \item{args}{ - returns the arguments to hlt}
#' }
#' 
#' 
#' 
#' @examples 
#' 
#' # set seed
#' set.seed(153)
#' 
#' # load the asti data set
#' data("asti")
#' 
#' # shift responses to range from 0 instead of 1
#' x = as.matrix(asti[, 1:25]) - 1
#' 
#' # subset and transform predictor data
#' z = asti[, 26:27]
#' z[, 1] = (z[, 1] == "students") * 1
#' z[, 2] = (z[, 2] == "male") * 1
#' z = as.matrix(z)
#' 
#' # specify which items from x belong to each domain
#' id = c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4)
#' 
#' # fit the model
#' mod = hlt(x, z = z, id = id, iter = 20, burn = 10, delta = 0.002)
#' 
#' mod$accept.rate # ideally 0.234
#' 
#' \donttest{
#' plot(mod, param = "lambda1", type = "trace")
#' plot(mod, param = "lambda2", type = "trace")
#' plot(mod, param = "lambda3", type = "trace")
#' plot(mod, param = "a1", type = "trace")
#' plot(mod, param = "d2_3", type = "trace")
#' plot(mod, param = "beta1", type = "trace")
#' 
#' plot(mod, item = 1, type = "icc")
#' plot(mod, item = 2, type = "icc")
#' plot(mod, item = 3, type = "icc")
#' plot(mod, item = 4, type = "icc")
#' plot(mod, item = 5, type = "icc")
#' plot(mod, item = 6, type = "icc")
#' plot(mod, item = 7, type = "icc", min = -10, max = 10)
#' 
#' summary(mod, param = "all")
#' summary(mod, param = "delta", digits = 2)
#' summary(mod, param = "lambda")
#' summary(mod, param = "alpha")
#' summary(mod, param = "delta")
#' summary(mod, param = "theta", dimension = 1)
#' summary(mod, param = "theta", dimension = 2)
#' summary(mod, param = "theta", dimension = 3)
#' summary(mod, param = "theta", dimension = 4)
#' 
#' # start from a previous run's solution
#' post = tail(mod$post, 1)
#' start = list(lambda = post[1, c("lambda1", "lambda2", "lambda3", "lambda4", "lambda5")],
#'              theta = mod$theta_mean, 
#'              delta = post[1, grepl("^[d]", colnames(post))], 
#'              alpha = post[1, paste0("a", 1:25)], 
#'              beta = post[1, c("beta1", "beta2")])
#'              
#' mod = hlt(x, z = z, id = id, start = start, iter = 20, burn = 10, delta = 0.002)
#' 
#' 
#' }
#' 
#' @importFrom stats cor quantile rnorm runif sd
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doParallel registerDoParallel
#' @useDynLib hlt, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppDist
#' @import RcppProgress
#' @export
hlt = function(x, 
               z = NULL, 
               id, 
               iter, 
               burn = iter / 2, 
               delta,
               type = "2p",
               start = list(
                 list(lambda = c(), theta = c(), delta = c(), 
                      alpha = c(), beta = c())
                 ),
               nchains = 1,
               progress = TRUE,
               verbose = FALSE) {
  
    if(nchains > 1) {
      registerDoParallel(cores = nchains)
      
      if(length(start) < nchains) {
        models = foreach (i = 1:nchains, .verbose = verbose) %dopar% {
          hlt(x = x, z = z, id = id, iter = iter, burn = burn, delta = delta, 
              type = type, progress = FALSE, nchains = 1)
        }
      } else {
        models = foreach (i = 1:nchains, .verbose = verbose) %dopar% {
          hlt(x = x, z = z, id = id, iter = iter, burn = burn, delta = delta, 
              type = type, progress = FALSE, nchains = 1, start = start[[i]])
        }
      }
      
      names(models) = paste0("chain", 1:nchains)
      class(models) = "hltObjList"
      return(models)
    }
  
    if(!is.matrix(x)) {
      x = as.matrix(x)
    }
    
    ntheta = length(unique(id))
    
    if(ntheta < 2) {
      stop("Specified ntheta < 2. Must assume at least two latent dimensions
           to perform inference.")
    }
    
    if(ntheta == 2) {
      warning("Specified ntheta == 2. The lambda parameters for the two 
              specified dimensions are set to be equal, i.e. lambda1 == lambda2.
              This constraint is only required with < 3 dimensions.")
    }
    
    if(length(id) != ncol(x)) {
      stop("The id vector must match number of columns in x.")
    }
  
    n = nrow(x)
    nT = ntheta + 1
    J = ncol(x)
    lJ = apply(x, 2, function(x) {length(unique(x))})
    nD = max(lJ) - 1
    
    isZ = !is.null(z)
    
    if(isZ) {
      
      if(!is.matrix(z)) {
        z = as.matrix(z)
      }
      
      nB = ncol(z)
      
      if(type == "1p") {
        
        npar_theta = n*nT
        npar = (nT - 1) + nD*J + nB
        npar_total = npar + npar_theta
        
        draw = matrix(nrow = 1, ncol = npar)
        draw_theta = matrix(nrow = 1, ncol = npar_theta)
        
        theta_mean = matrix(0, nrow = 1, ncol = npar_theta)
        theta_mean_sq = matrix(0, nrow = 1, ncol = npar_theta)
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                       paste0("beta", 1:nB))
        
        post_names_theta = post_names[grepl("theta", post_names)]
        post_names_no_theta = post_names[!grepl("theta", post_names)]
        
        post = matrix(0, nrow = iter - burn, ncol = npar)
        colnames(post) = post_names_no_theta
        
        start_values = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                                 nD = nD, J = J, nB = nB)
        
        draw[1, ] = start_values[post_names %in% post_names_no_theta]
        draw_theta[1, ] = start_values[post_names %in% post_names_theta]
        
      } else if(type == "2p") {
        
        npar_theta = n*nT
        npar = (nT - 1) + J + nD*J + nB
        npar_total = npar + npar_theta
        
        draw = matrix(nrow = 1, ncol = npar)
        draw_theta = matrix(nrow = 1, ncol = npar_theta)
        
        theta_mean = matrix(0, nrow = 1, ncol = npar_theta)
        theta_mean_sq = matrix(0, nrow = 1, ncol = npar_theta)
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                       paste0("a", 1:J),
                       paste0("beta", 1:nB))
        
        post_names_theta = post_names[grepl("theta", post_names)]
        post_names_no_theta = post_names[!grepl("theta", post_names)]
        
        post = matrix(0, nrow = iter - burn, ncol = npar)
        colnames(post) = post_names_no_theta
        
        start_values = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                                 nD = nD, J = J, nB = nB)
        
        draw[1, ] = start_values[post_names %in% post_names_no_theta]
        draw_theta[1, ] = start_values[post_names %in% post_names_theta]
        
      }
      
    } 
    else {
      
      if(type == "1p") {
        
        npar_theta = n*nT
        npar = (nT - 1) + nD*J
        npar_total = npar + npar_theta
        
        draw = matrix(nrow = 1, ncol = npar)
        draw_theta = matrix(nrow = 1, ncol = npar_theta)
        
        theta_mean = matrix(0, nrow = 1, ncol = npar_theta)
        theta_mean_sq = matrix(0, nrow = 1, ncol = npar_theta)
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})))
        
        post_names_theta = post_names[grepl("theta", post_names)]
        post_names_no_theta = post_names[!grepl("theta", post_names)]
        
        post = matrix(0, nrow = iter - burn, ncol = npar)
        colnames(post) = post_names_no_theta
        
        start_values = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                                 nD = nD, J = J, nB = nB)
        
        draw[1, ] = start_values[post_names %in% post_names_no_theta]
        draw_theta[1, ] = start_values[post_names %in% post_names_theta]
        
      } else if(type == "2p") {
        
        npar_theta = n*nT
        npar = (nT - 1) + J + nD*J
        npar_total = npar + npar_theta
        
        draw = matrix(nrow = 1, ncol = npar)
        draw_theta = matrix(nrow = 1, ncol = npar_theta)
        
        theta_mean = matrix(0, nrow = 1, ncol = npar_theta)
        theta_mean_sq = matrix(0, nrow = 1, ncol = npar_theta)
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                       paste0("a", 1:J))
        
        post_names_theta = post_names[grepl("theta", post_names)]
        post_names_no_theta = post_names[!grepl("theta", post_names)]
        
        post = matrix(0, nrow = iter - burn, ncol = npar)
        colnames(post) = post_names_no_theta
        
        start_values = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                                 nD = nD, J = J, nB = nB)
        
        draw[1, ] = start_values[post_names %in% post_names_no_theta]
        draw_theta[1, ] = start_values[post_names %in% post_names_theta]
        
      }
    }
    
    c.ix = calc.ix(post_names_no_theta, npar)
    ix = c.ix$ix
    ixe = c.ix$ixe
    
    accept = numeric(iter)
    accept[1] = 1
    
    if(isZ) {
      
      if(type == "1p") {
        
        if(ntheta == 2) {
          
          lt1PR2D(x = x,
                  z = z,
                  iter = iter,
                  burn = burn,
                  delta = delta,
                  post = post,
                  mean_theta = theta_mean,
                  mean_theta_sq = theta_mean_sq,
                  draw = draw,
                  draw_theta = draw_theta,
                  ix = ix,
                  ixe = ixe,
                  npar = npar,
                  ntheta = n * nT,
                  n = n,
                  nB = nB,
                  J = J,
                  nDmax = nD,
                  lJ = lJ,
                  nT = nT,
                  tJ = id,
                  accept = accept,
                  eps = .Machine$double.eps,
                  display_progress = progress)
          
        } else {
          
          lt1PR(x = x,
                z = z,
                iter = iter,
                burn = burn,
                delta = delta,
                post = post,
                mean_theta = theta_mean,
                mean_theta_sq = theta_mean_sq,
                draw = draw,
                draw_theta = draw_theta,
                ix = ix,
                ixe = ixe,
                npar = npar,
                ntheta = n * nT,
                n = n,
                nB = nB,
                J = J,
                nDmax = nD,
                lJ = lJ,
                nT = nT,
                tJ = id,
                accept = accept,
                eps = .Machine$double.eps,
                display_progress = progress)
          
        }
        
      } else if(type == "2p") {
        
        if(ntheta == 2) {
          
          lt2PR2D(x = x,
                  z = z,
                  iter = iter,
                  burn = burn,
                  delta = delta,
                  post = post,
                  mean_theta = theta_mean,
                  mean_theta_sq = theta_mean_sq,
                  draw = draw,
                  draw_theta = draw_theta,
                  ix = ix,
                  ixe = ixe,
                  npar = npar,
                  ntheta = n * nT,
                  n = n,
                  nB = nB,
                  J = J,
                  nDmax = nD,
                  lJ = lJ,
                  nT = nT,
                  tJ = id,
                  accept = accept,
                  eps = .Machine$double.eps,
                  display_progress = progress)
          
        } else {
          
          lt2PR(x = x,
                z = z,
                iter = iter,
                burn = burn,
                delta = delta,
                post = post,
                mean_theta = theta_mean,
                mean_theta_sq = theta_mean_sq,
                draw = draw,
                draw_theta = draw_theta,
                ix = ix,
                ixe = ixe,
                npar = npar,
                ntheta = n * nT,
                n = n,
                nB = nB,
                J = J,
                nDmax = nD,
                lJ = lJ,
                nT = nT,
                tJ = id,
                accept = accept,
                eps = .Machine$double.eps,
                display_progress = progress)
          
        }
        
      }
      
    } 
    else {
      
      if(type == "1p") {
        
        if(ntheta == 2) {
          
          lt1PNR2D(x = x,
                   iter = iter,
                   burn = burn,
                   delta = delta,
                   post = post,
                   mean_theta = theta_mean,
                   mean_theta_sq = theta_mean_sq,
                   draw = draw,
                   draw_theta = draw_theta,
                   ix = ix,
                   ixe = ixe,
                   npar = npar,
                   ntheta = n * nT,
                   n = n,
                   J = J,
                   nDmax = nD,
                   lJ = lJ,
                   nT = nT,
                   tJ = id,
                   accept = accept,
                   eps = .Machine$double.eps,
                   display_progress = progress)
          
        } else {
          
          lt1PNR(x = x,
                 iter = iter,
                 burn = burn,
                 delta = delta,
                 post = post,
                 mean_theta = theta_mean,
                 mean_theta_sq = theta_mean_sq,
                 draw = draw,
                 draw_theta = draw_theta,
                 ix = ix,
                 ixe = ixe,
                 npar = npar,
                 ntheta = n * nT,
                 n = n,
                 J = J,
                 nDmax = nD,
                 lJ = lJ,
                 nT = nT,
                 tJ = id,
                 accept = accept,
                 eps = .Machine$double.eps,
                 display_progress = progress)
          
        }
      } else if(type == "2p") {
        
        if(ntheta == 2) {
          lt2PNR2D(x = x,
                   iter = iter,
                   burn = burn,
                   delta = delta,
                   post = post,
                   mean_theta = theta_mean,
                   mean_theta_sq = theta_mean_sq,
                   draw = draw,
                   draw_theta = draw_theta,
                   ix = ix,
                   ixe = ixe,
                   npar = npar,
                   ntheta = n * nT,
                   n = n,
                   J = J,
                   nDmax = nD,
                   lJ = lJ,
                   nT = nT,
                   tJ = id,
                   accept = accept,
                   eps = .Machine$double.eps,
                   display_progress = progress)
          
        } else {
          
          lt2PNR(x = x,
                 iter = iter,
                 burn = burn,
                 delta = delta,
                 post = post,
                 mean_theta = theta_mean,
                 mean_theta_sq = theta_mean_sq,
                 draw = draw,
                 draw_theta = draw_theta,
                 ix = ix,
                 ixe = ixe,
                 npar = npar,
                 ntheta = n * nT,
                 n = n,
                 J = J,
                 nDmax = nD,
                 lJ = lJ,
                 nT = nT,
                 tJ = id,
                 accept = accept,
                 eps = .Machine$double.eps,
                 display_progress = progress)
          
        }
      }
    }
    
    d_end = cumsum(rep(nD, J))
    d_start = d_end - nD + 1
    d_00 = max(lJ - 1) - (lJ - 1)
    d_null = d_end[d_00 > 0] + ix[2]
    d_null_start = d_null - d_00[d_00 > 0]
    d_null = d_null - 1
    
    if(length(d_null) > 0) {
      for(j in 1:length(d_null)) {
        post[, d_null_start[j]:d_null[j]] <- 0
      }
    }
    
    theta_mean = theta_mean[1, ] / (iter - burn)
    theta_mean_sq = theta_mean_sq[1, ] / (iter - burn)
    theta_sd = sqrt(theta_mean_sq - (theta_mean^2))
    
    colnames(x) = paste0("x", 1:J)
    accept.rate = mean(accept)
    theta = data.frame(mean = theta_mean, sd = theta_sd)
    
    result = list(post = post, 
                  accept.rate = accept.rate, 
                  theta = theta,
                  nT = nT,
                  args = list(x = x, z = z, id = id, iter = iter, burn = burn, 
                              delta = delta, type = type, start = start,
                              progress = progress))
    
    class(result) = c("hltObj")
    return(result)
}
