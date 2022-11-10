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
inirt = function(data, model = NULL, predictors = NULL, dims, h2_dims = 0, h2_dim_id = NULL, method = c("vb", "hmc"), weights = NULL, ...) {
  if(is.null(predictors)) {
    irt_data = data[, drop = FALSE]
    N = nrow(irt_data)
    J = ncol(irt_data)
    if(any(is.na(irt_data))) {
      model_missing_y = 1
    } else {
      model_missing_y = 0
    }
  } else {
    predictor_ulist = unlist(predictors)
    irt_data = data[, -predictor_ulist, drop = FALSE]
    reg_data = data[, predictor_ulist, drop = FALSE]
    N = nrow(irt_data)
    J = ncol(irt_data)
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
  
  if(model_missing_y == 0) {
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
  if(h2_dims == 0) {
    L = D*(J-D) + D*(D+1)/2
  } else if(h2_dims == 1){
    L = J
  } else {
    stop("h2_dims > 1 not yet implemented. Choose 0 or 1.")
  }
  nDelta = sum(Ncategi - 1)

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

    beta_dstart = array(beta_dstart, dim = D)
    beta_dend = array(beta_dend, dim = D)
    nobeta_dstart = array(nobeta_dstart, dim = D)
    nobeta_dend = array(nobeta_dend, dim = D)

    if(model_missing_x == 0) {

    } else {
      Lxmiss = sum(is.na(reg_data))
      x_miss_id = 1:Lxmiss
      reg_miss = is.na(reg_data) * 1
      reg_miss = as.vector(t(reg_miss))
      reg_miss[reg_miss == 1] = x_miss_id
      reg_miss = matrix(reg_miss, nrow = N, byrow = TRUE)
      x_miss = is.na(x) * 1
      x_in_row_is_missing = apply(x, 1, function(x) {any(is.na(x))}) * 1
      x[is.na(x)] = 0
    }
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

  if(h2_dims > 0) {
    # h2_dim_id = list(c(1,2,4,5,7),c(3,6,9,11),c(8,10,12,13,14,15))
    h2_dim_id_ulist = unlist(h2_dim_id)
    d_lengths = sapply(h2_dim_id, function(x) {length(x)})
    d_seq_start = 1
    new_id_list = vector(mode = "list", length = D)
    for(i in 1:D) {
      new_id_list[[i]] = d_seq_start:(d_lengths[i] + d_seq_start - 1)
      d_seq_start = max(new_id_list[[i]]) + 1
    }
    new_id_list = lapply(new_id_list, function(x) {c(x[1], x[length(x)])})
    irt_data = irt_data[, h2_dim_id_ulist]
    alpha_dstart = array(sapply(new_id_list, function(x) {x[1]}), dim = D)
    alpha_dend = array(sapply(new_id_list, function(x) {x[length(x)]}), dim = D)
    lambda_ind = c()
    for(i in 1:D) {
      lambda_ind = c(lambda_ind, rep(i, d_lengths[i]))
    }
    lambda_ind = as.vector(t(replicate(N, lambda_ind)))
  }

  if(D == 1) { # Data input for unidimensional model
    if(regress == 0) { # No Regression
      if(model_missing_y == 0) {
        if(model_missing_x == 0) {

        } else {

        }
        standata = list(N = N, J = J, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, 
                      D = D, nDelta = nDelta, L = L, weights = weights)
      } else if(model_missing_y == 1) {
        standata = list(N = N, J = J, Ncateg = Ncateg, N_long_obs = N_long_obs, nn = nn, jj = jj, 
          y = y, D = D, L = L, weights = weights)
      }
    } else { # Regression
      if(model_missing_x == 0) {
        standata = list(N = N, J = J, K = K, Ncateg = Ncateg, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
          D = D, L = L, Lbeta = Lbeta, beta_dstart = beta_dstart, beta_dend = beta_dend)
      } else {
        standata = list(N = N, J = J, K = K, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
          D = D, L = L, Lbeta = Lbeta, beta_dstart = beta_dstart, beta_dend = beta_dend, Lxmiss = Lxmiss,
          x_miss = x_miss, reg_miss = reg_miss, x_in_row_is_missing = x_in_row_is_missing)
      }
    }
  } else if(D > 1 & h2_dims == 0) { # Data input for multidimensional model
    standata = list(N = N, J = J, Ncateg = Ncateg, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
        D = D, L = L, Lbeta = Lbeta, beta_dstart = beta_dstart, beta_dend = beta_dend, 
        nobeta_dstart = nobeta_dstart, nobeta_dend = nobeta_dend)
  } else if(D > 1 & h2_dims > 0) { # second-order HO-IRT model (single second-order dimension)
    standata = list(N = N, J = J, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, 
        y = y, D = D, nDelta = nDelta, L = L, alpha_dstart = alpha_dstart, alpha_dend = alpha_dend, 
        lambda_ind = lambda_ind, weights = weights)
  }
  # Select the correct inirt implementation based on input parameters
  if(D == 1) { # unidimensional analysis
    if(regress == 0) {
      if(model_missing_y == 0) {
        modl = stanmodels$inirt_mirt_unidim_no_tregress
      } else {
        modl = stanmodels$inirt_mirt_unidim_no_tregress_ymiss
      }
    } else if(regress == 1) {
      if(model_missing_y == 0) {
        if(model_missing_x == 0) {
          modl = stanmodels$inirt_mirt_unidim_tregress
        } else {
          modl = stanmodels$inirt_mirt_unidim_tregress_xmiss
        }
      } else {
        if(model_missing_x == 0) {
          modl = stanmodels$inirt_mirt_unidim_tregress_ymiss
        } else {
          modl = stanmodels$inirt_mirt_unidim_tregress_xmiss_ymiss
        }
      }
    } 
  } else if(D > 1 & h2_dims == 0) { # multidimensional analysis with no higher-order trait
    if(regress == 0) {
      if(model_missing_y == 0) {
        modl = stanmodels$inirt_mirt_no_tregress
      } else {
        modl = stanmodels$inirt_mirt_no_tregress_ymiss
      }
    } else {
      if(model_missing_y == 0) {
        if(model_missing_x == 0) {
          modl = stanmodels$inirt_mirt_tregress
        } else {
          modl = stanmodels$inirt_mirt_tregress_xmiss
        }
      } else {
        if(model_missing_x == 0) {
          modl = stanmodels$inirt_mirt_tregress_ymiss
        } else {
          modl = stanmodels$inirt_mirt_tregress_xmiss_ymiss
        }
      }
    }
  } else if (D > 1 & h2_dims > 0){ # second order analysis
    if(regress == 0) {
      if(model_missing_y == 0) {
        modl = stanmodels$inirt_hoirt_2ord_unidim_no_tregress
      } else {

      }
    } else {
      if(model_missing_y == 0) {
        if(model_missing_x == 0) {

        } else {

        }
      } else {
        if(model_missing_x == 0) {

        } else {
          
        }
      }
    }
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