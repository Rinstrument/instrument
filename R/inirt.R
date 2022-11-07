#' Explanatory Item Response Theory
#'
#' @export
#' @param data data frame
#' @param model text specification of the model. Not yet implemented. Instead use the other arguments
#' @param predictors list of predictors
#' @param dims dimensions
#' @param method Choose estimation method. Choices are method = 'vb' (default) for variational Bayes, or method = 'hmc' for Hamiltonian Monte Carlo.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' 
#'
inirt = function(data, model = NULL, predictors = NULL, dims, method = c("vb", "hmc"), weights = NULL, ...) {
  if(is.null(predictors)) {
    irt_data = data[, drop = FALSE]
    if(any(is.na(irt_data))) {
      model_missing_y = 1
    } else {
      model_missing_y = 0
    }
  } else {
    predictor_ulist = unlist(predictors)
    irt_data = data[, -predictor_ulist, drop = FALSE]
    reg_data = data[, predictor_ulist, drop = FALSE]
    K = ncol(reg_data)
    x = matrix(rep(t(reg_data), J), ncol = K, byrow = TRUE)
    Lbeta = length(predictor_ulist)
    if(any(is.na(irt_data))) {
      model_missing_y = 1
    } else {
      model_missing_y = 0
    }
    if(any(is.na(reg_data))) {
      model_missing_x = 1
    } else {
      model_missing_x = 0
    }
  }
  N = nrow(irt_data)
  J = ncol(irt_data)
  if(model_missin_y == 0) {
    Ncateg_max = max(irt_data)
    Ncategi = apply(irt_data, 2, max)
    names(Ncategi) = NULL
  } else if(model_missing_y == 1) {
    Ncateg_max = max(as.vector(irt_data)[!is.na(as.vector(irt_data))])
    Ncategi = apply(irt_data, 2, max) # neex to update this to account for missingness!
    names(Ncategi) = NULL
  }
  N_long = N*J
  y = as.vector(irt_data)
  nn = rep(1:N, J)
  jj = rep(1:J, each = N)
  D = dims
  L = D*(J-D) + D*(D+1)/2
  nDelta = sum(Ncategi - 1)
  # bindex = 1
  # beta_dstart = bindex
  # beta_dend = bindex + length(predictors[[1]]) - 1
  # for(d in 1:D) {
  #   beta_dstart = c(beta_dstart, bindex + length(predictors[[d]]))
  #   beta_dend = c(beta_dend, beta_dstart[d + 1] + length(predictors[[d + 1]]) - 1)
  # }
  if(is.null(predictors)) {
    regress = 0
  } else {
    regress = 1
    start_index = 1
    beta_dstart = numeric(D)
    beta_dend = numeric(D)
    for(d in 1:D) {
      beta_dstart[d] = start_index
      beta_dend[d] = start_index + length(predictors[[d]]) - 1
      start_index = start_index + length(predictors[[d]])
    }
    nobeta_dstart = numeric(D) #1000
    nobeta_dend = numeric(D) #1
    beta_dstart_rep = rep(beta_dstart, 2)
    beta_dend_rep = rep(beta_dend, 2)
    
    nobeta_dstart = beta_dstart_rep[c(-1, -length(beta_dstart_rep))]
    nobeta_dend = beta_dend_rep[c(-1, -length(beta_dend_rep))]

    beta_dstart = array(beta_dstart, dim = d)
    beta_dend = array(beta_dend, dim = d)
    nobeta_dstart = array(nobeta_dstart, dim = d)
    nobeta_dend = array(nobeta_dend, dim = d)
  }
  if(model_missing_y == 1) {
    N_miss = sum(is.na(irt_data))
    N_long_obs = N_long - N_miss
    nn = nn[!is.na(y)]
    jj = jj[!is.na(y)]
    y = y[!is.na(y)]
    if(is.null(weights)) {
      weights = rep(1, N_long_obs)
    }
  } else {
    if(is.null(weights)) {
      weights = rep(1, N_long)
    }
  }
  if(D == 1) { # Data input for unidimensional model
    if(regress == 0) { # No Regression
      if(model_missing_y == 0) {
        standata = list(N = N, J = J, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, 
                      D = D, nDelta = nDelta, L = L, weights = weights)
      } else if(model_missing_y == 1) {
        standata = list(N = N, J = J, Ncateg = Ncateg, N_long_obs = N_long_obs, nn = nn, jj = jj, 
          y = y, D = D, L = L, weights = weights)
      }
    } else { # Regression
      standata = list(N = N, J = J, K = K, Ncateg = Ncateg, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
          D = D, L = L, Lbeta = Lbeta, beta_dstart = beta_dstart, beta_dend = beta_dend)
    }
  } else { # Data input for multidimensional model
    standata = list(N = N, J = J, K = K, Ncateg = Ncateg, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
        D = D, L = L, Lbeta = Lbeta, beta_dstart = beta_dstart, beta_dend = beta_dend, 
        nobeta_dstart = nobeta_dstart, nobeta_dend = nobeta_dend)
  }
  # Select the correct inirt implementation based on input parameters
  if(D == 1) {
    if(regress == 0) {
      if(model_missing_y == 0) {
        modl = stanmodels$inirt_mirt_unidim_noRegress
      } else {
        modl = stanmodels$inirt_mirt_unidim_noRegress_ymiss
      }
    } else {
      modl = stanmodels$inirt_mirt_unidim
    }
  } else if (D > 1){
    modl = stanmodels$inirt_mirt
  }
  # Choose a model estimateion method: variational inference or HMC
  if(method[1] == "vb") {
    out = rstan::vb(modl, data = standata, ...)
  } else if(method[1] == "hmc") {
    out = rstan::sampling(modl, data = standata, ...)
  } else {
    stop("Method not yet implemeted.")
  }
  # Return the output
  # class(out) = c("inirt", class(out))
  return(out)
}