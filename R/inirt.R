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
inirt = function(data, item_id, model = NULL, predictors = NULL, predictors_ranef = NULL, ranef_id = NULL, 
    predictors_ranef_corr = NULL, n_pranef_cor = NULL, dims, h2_dims = 0, h2_dim_id = NULL, structural_design = NULL, 
    method = c("vb", "hmc"), weights = NULL, ...) {
  irt_data = data[, item_id, drop = FALSE]
  N = nrow(irt_data)
  J = ncol(irt_data)
  if(any(is.na(irt_data))) {
    model_missing_y = 1
  } else {
    model_missing_y = 0
  }
  has_treg = 0
  any_rand = 0
  any_rand_ind = 0
  any_rand_cor = 0
  model_missing_x = 0
  model_missing_x_rand = 0
  if(!is.null(predictors)) {
    has_treg = 1
    predictor_ulist = unlist(predictors)
    reg_data = data[, predictor_ulist, drop = FALSE]
    K = ncol(reg_data)
    x = matrix(rep(t(reg_data), J), ncol = K, byrow = TRUE)
    if(any(is.na(reg_data))) {
        model_missing_x = 1
    }
  }
  if(!is.null(predictors_ranef)) {
    has_treg = 1
    any_rand = 1
    any_rand_ind = 1
    predictor_ranef_ulist = unlist(predictors_ranef)
    reg_data_ranef = data[, predictor_ranef_ulist, drop = FALSE]
    Lzeta = ncol(reg_data_ranef)
    z = matrix(rep(t(reg_data_ranef), J), ncol = Lzeta, byrow = TRUE)
    if(any(is.na(reg_data_ranef))) {
        model_missing_x_rand = 1
    }
  }
  if(!is.null(predictors_ranef_corr)) {
    has_treg = 1
    any_rand = 1
    any_rand_cor = 1
    predictors_ranef_corr_ulist = unlist(predictors_ranef_corr)
    reg_data_ranef_cor = data[, predictors_ranef_corr_ulist, drop = FALSE]
    Lzeta_cor = ncol(reg_data_ranef_cor)
    z_c = matrix(rep(t(reg_data_ranef_cor), J), ncol = Lzeta_cor, byrow = TRUE)
    l_Lzeta_cor = n_pranef_cor
    u_Lzeta_cor = Lzeta_cor / l_Lzeta_cor
    cor_z_item_ind = rep(1:u_Lzeta_cor, l_Lzeta_cor)
    cor_z_item_elem_ind = rep(1:l_Lzeta_cor, each = u_Lzeta_cor)
  }

  # if(is.null(predictors) & is.null(predictors_ranef) & is.null(predictors_ranef_corr)) {
  #   model_missing_x = 0
  #   has_treg = 0
  #   any_rand = 0
  #   any_rand_ind = 0
  #   any_rand_cor = 0
  #   Lzeta_cor = 0
  # } else {
  #   if(!is.null(predictors) & !is.null(predictors_ranef) & !is.null(predictors_ranef_corr)) {

  #   } else if(

  #   )
  #   if(!is.null(predictors) & !is.null(predictors_ranef)) {
  #     has_treg = 1
  #     any_rand = 1
  #     any_rand_ind = 1
  #     any_rand_cor = 0
  #     predictor_ulist = unlist(predictors)
  #     predictor_ranef_ulist = unlist(predictors_ranef)
  #     reg_data = data[, predictor_ulist, drop = FALSE]
  #     reg_data_ranef = data[, predictor_ranef_ulist, drop = FALSE]
  #     K = ncol(reg_data)
  #     Lzeta = ncol(reg_data_ranef)
  #     x = matrix(rep(t(reg_data), J), ncol = K, byrow = TRUE)
  #     z = matrix(rep(t(reg_data_ranef), J), ncol = Lzeta, byrow = TRUE)
  #     if(any(is.na(reg_data)) | any(is.na(reg_data_ranef))) {
  #       model_missing_x = 1
  #     } else {
  #       model_missing_x = 0
  #     }
  #   } else if(!is.null(predictors) & !is.null(predictors_ranef_corr)) {
  #     has_treg = 1
  #     any_rand = 1
  #     any_rand_ind = 0
  #     any_rand_cor = 1
  #     predictor_ulist = unlist(predictors)
  #     predictors_ranef_corr_ulist = unlist(predictors_ranef_corr)
  #     reg_data = data[, predictor_ulist, drop = FALSE]
  #     reg_data_ranef_cor = data[, predictors_ranef_corr_ulist, drop = FALSE]
  #     K = ncol(reg_data)
  #     Lzeta_cor = ncol(reg_data_ranef_cor)
  #     x = matrix(rep(t(reg_data), J), ncol = K, byrow = TRUE)
  #     z_c = matrix(rep(t(reg_data_ranef_cor), J), ncol = Lzeta_cor, byrow = TRUE)
  #     l_Lzeta_cor = n_pranef_cor
  #     u_Lzeta_cor = Lzeta_cor / l_Lzeta_cor
  #     cor_z_item_ind = rep(1:u_Lzeta_cor, l_Lzeta_cor)
  #     cor_z_item_elem_ind = rep(1:l_Lzeta_cor, each = u_Lzeta_cor)
  #     # if(any(is.na(reg_data)) | any(is.na(reg_data_ranef))) {
  #     #   model_missing_x = 1
  #     # } else {
  #     #   model_missing_x = 0
  #     # } missing data in the random effect design?
  #     model_missing_x = 0

  #   } else if(!is.null(predictors_ranef) & is.null(predictors)) {
  #     any_rand = 1
  #     any_rand_ind = 1
  #     any_rand_cor = 0
  #     predictor_ranef_ulist = unlist(predictors_ranef)
  #     reg_data = data[, predictor_ulist, drop = FALSE]
  #   } else {
  #     any_rand = 0
  #     any_rand_ind = 0
  #     any_rand_cor = 0
  #   }
  # }
  
  if(model_missing_y == 0) {
    Ncateg_max = max(irt_data)
    Ncategi = as.integer(apply(irt_data, 2, max))
    names(Ncategi) = NULL
  } else {
    Ncateg_max = max(as.vector(irt_data)[!is.na(as.vector(irt_data))])
    Ncategi = as.integer(apply(irt_data, 2, max)) # neex to update this to account for missingness!
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

  regress = 0
  if(!is.null(predictors)) {
    regress = 1
    start_index = 1
    if(h2_dims > 0) {
      beta_dstart = numeric(1)
      beta_dend = numeric(1)
      beta_dstart[1] = start_index
      beta_dend[1] = start_index + length(predictors[[1]]) - 1
      start_index = start_index + length(predictors[[1]])
      # beta_dstart = array(beta_dstart, dim = 1)
      # beta_dend = array(beta_dend, dim = 1)
    }
    if(h2_dims == 0) {
      start_index = 1
      beta_dstart = numeric(D)
      beta_dend = numeric(D)
      for(d in 1:D) {
        beta_dstart[d] = start_index
        beta_dend[d] = start_index + length(predictors[[d]]) - 1
        start_index = start_index + length(predictors[[d]])
      }
      beta_dstart = array(beta_dstart, dim = D)
      beta_dend = array(beta_dend, dim = D)
    }

    if(model_missing_x == 0) {
      reg_miss = is.na(reg_data) * 1
      x_miss = matrix(rep(t(reg_miss), J), ncol = K, byrow = TRUE)
      x_in_row_is_missing = apply(x, 1, function(x) {any(is.na(x))}) * 1
    }
    if(model_missing_x > 0) {
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
  if(!is.null(predictors_ranef)) {
    regress = 1
    start_index = 1
    zeta_dstart = numeric(D)
    zeta_dend = numeric(D)
    for(d in 1:D) {
      zeta_dstart[d] = start_index
      zeta_dend[d] = start_index + length(predictors_ranef[[d]]) - 1
      start_index = start_index + length(predictors_ranef[[d]])
    }
    zeta_dstart = array(zeta_dstart, dim = D)
    zeta_dend = array(zeta_dend, dim = D)
    Lzeta_sd = length(unique(ranef_id))
    zeta_sd_ind = ranef_id
  } else {
    Lzeta = 1
    Lzeta_sd = 0
    zeta_sd_ind = array(0, dim = 1)
    z = matrix(0, nrow = N_long, ncol = Lzeta)
    zeta_dstart = array(1, dim = D)
    zeta_dend = array(1, dim = D)
  }
  if(!is.null(predictors_ranef_corr)) {
    regress = 1
  }

  # if(is.null(predictors) & is.null(predictors_ranef) & is.null(predictors_ranef_corr)) {
  #   regress = 0
  # } else {
  #   if(!is.null(predictors) & !is.null(predictors_ranef)) {
  #     regress = 1
  #     start_index = 1
  #     if(h2_dims > 0) {
  #       beta_dstart = numeric(1)
  #       beta_dend = numeric(1)
  #       beta_dstart[1] = start_index
  #       beta_dend[1] = start_index + length(predictors[[1]]) - 1
  #       start_index = start_index + length(predictors[[1]])
  #       # beta_dstart = array(beta_dstart, dim = 1)
  #       # beta_dend = array(beta_dend, dim = 1)
  #     } else {
  #       start_index = 1
  #       beta_dstart = numeric(D)
  #       beta_dend = numeric(D)
  #       for(d in 1:D) {
  #         beta_dstart[d] = start_index
  #         beta_dend[d] = start_index + length(predictors[[d]]) - 1
  #         start_index = start_index + length(predictors[[d]])
  #       }
  #       beta_dstart = array(beta_dstart, dim = D)
  #       beta_dend = array(beta_dend, dim = D)

  #       start_index = 1
  #       zeta_dstart = numeric(D)
  #       zeta_dend = numeric(D)
  #       for(d in 1:D) {
  #         zeta_dstart[d] = start_index
  #         zeta_dend[d] = start_index + length(predictors_ranef[[d]]) - 1
  #         start_index = start_index + length(predictors_ranef[[d]])
  #       }
  #       zeta_dstart = array(zeta_dstart, dim = D)
  #       zeta_dend = array(zeta_dend, dim = D)

  #       Lzeta_sd = length(unique(ranef_id))
  #       zeta_sd_ind = ranef_id
  #     }

  #     if(model_missing_x == 0) {
  #       reg_miss = is.na(reg_data) * 1
  #       x_miss = matrix(rep(t(reg_miss), J), ncol = K, byrow = TRUE)
  #       x_in_row_is_missing = apply(x, 1, function(x) {any(is.na(x))}) * 1
  #     } else {
  #       Lxmiss = sum(is.na(reg_data))
  #       x_miss_id = 1:Lxmiss
  #       reg_miss = is.na(reg_data) * 1
  #       reg_miss = as.vector(t(reg_miss))
  #       reg_miss[reg_miss == 1] = x_miss_id
  #       reg_miss = matrix(reg_miss, nrow = N, byrow = TRUE)
  #       x_miss = is.na(x) * 1
  #       x_in_row_is_missing = apply(x, 1, function(x) {any(is.na(x))}) * 1
  #       x[is.na(x)] = 0
  #     }
  #   } else if(!is.null(predictors) & !is.null(predictors_ranef_corr)) {
  #     regress = 1
  #     Lzeta = 1
  #     Lzeta_sd = 0
  #     zeta_sd_ind = array(0, dim = 1)
  #     z = matrix(0, nrow = N_long, ncol = Lzeta)
  #     start_index = 1
  #     if(h2_dims > 0) {
  #       # beta_dstart = numeric(1)
  #       # beta_dend = numeric(1)
  #       # beta_dstart[1] = start_index
  #       # beta_dend[1] = start_index + length(predictors[[1]]) - 1
  #       # start_index = start_index + length(predictors[[1]])
  #       # # beta_dstart = array(beta_dstart, dim = 1)
  #       # # beta_dend = array(beta_dend, dim = 1)
  #     } else {
  #       beta_dstart = numeric(D)
  #       beta_dend = numeric(D)
  #       for(d in 1:D) {
  #         beta_dstart[d] = start_index
  #         beta_dend[d] = start_index + length(predictors[[d]]) - 1
  #         start_index = start_index + length(predictors[[d]])
  #       }
  #       beta_dstart = array(beta_dstart, dim = D)
  #       beta_dend = array(beta_dend, dim = D)

  #       # for multiple dimensions
  #       start_index = 1
  #       # zeta_dstart = numeric(D)
  #       # zeta_dend = numeric(D)
  #       # for(d in 1:D) {
  #         # zeta_dstart[d] = start_index
  #         # zeta_dend[d] = start_index + length(predictors_ranef[[d]]) - 1
  #         # start_index = start_index + length(predictors_ranef[[d]])
  #       # }
  #       # zeta_dstart = array(zeta_dstart, dim = D)
  #       # zeta_dend = array(zeta_dend, dim = D)
  #       zeta_dstart = array(1, dim = D)
  #       zeta_dend = array(1, dim = D)
  #     }

  #     # if(model_missing_x == 0) {
  #     #   reg_miss = is.na(reg_data) * 1
  #     #   x_miss = matrix(rep(t(reg_miss), J), ncol = K, byrow = TRUE)
  #     #   x_in_row_is_missing = apply(x, 1, function(x) {any(is.na(x))}) * 1
  #     # } else {
  #     #   Lxmiss = sum(is.na(reg_data))
  #     #   x_miss_id = 1:Lxmiss
  #     #   reg_miss = is.na(reg_data) * 1
  #     #   reg_miss = as.vector(t(reg_miss))
  #     #   reg_miss[reg_miss == 1] = x_miss_id
  #     #   reg_miss = matrix(reg_miss, nrow = N, byrow = TRUE)
  #     #   x_miss = is.na(x) * 1
  #     #   x_in_row_is_missing = apply(x, 1, function(x) {any(is.na(x))}) * 1
  #     #   x[is.na(x)] = 0
  #     # }
  #   } 
  #   # else if(!is.null(predictors_ranef) & is.null(predictors)) {
      
  #   # } else {
      
  #   # }
  # }

  # sample weights
  if(is.null(weights)) {
    if(model_missing_y == 1) {
      weights = rep(1, N_long_obs)
    } else {
      weights = rep(1, N_long)
    }
  }

  # missing values in the item observation matrix
  if(model_missing_y == 1) {
    N_miss = sum(is.na(irt_data))
    N_long_obs = N_long - N_miss
    nn = nn[!is.na(y)]
    jj = jj[!is.na(y)]
    y = y[!is.na(y)]
  }

  if(h2_dims > 0) {
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

  # Set up the structural regressions (item regression)
  if(!is.null(structural_design)) {
    if(any(c("alpha") %in% names(structural_design))) {
      a_design = as.matrix(structural_design[["alpha"]])
      nAlpha_r = ncol(a_design)
    }
    if(any(c("delta") %in% names(structural_design))) {
      d_design = as.matrix(structural_design[["delta"]])
      nDelta_r = ncol(d_design)
    }
  }

  if(D == 1) {
    standata = list(N = N, J = J, K = K, any_rand = any_rand, any_rand_ind = any_rand_ind, any_rand_cor = any_rand_cor, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
        D = D, nDelta = nDelta, L = L, has_treg = has_treg, beta_dstart = beta_dstart, beta_dend = beta_dend, zeta_dstart = zeta_dstart, zeta_dend = zeta_dend, 
        weights = weights, x_miss = x_miss, reg_miss = reg_miss, x_in_row_is_missing = x_in_row_is_missing, nDelta_r = nDelta_r, nAlpha_r = nAlpha_r,
        d_design = d_design, a_design = a_design, Lzeta = Lzeta, u_Lzeta_cor = u_Lzeta_cor, l_Lzeta_cor = l_Lzeta_cor, Lzeta_cor = Lzeta_cor, z = z, 
        Lzeta_sd = Lzeta_sd, zeta_sd_ind = zeta_sd_ind, cor_z_item_ind = cor_z_item_ind, cor_z_item_elem_ind = cor_z_item_elem_ind, z_c = z_c)
  } else if(D > 1 & h2_dims == 0) {
    standata = list(N = N, J = J, K = K, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
        D = D, nDelta = nDelta, L = L, has_treg = has_treg, beta_dstart = beta_dstart, beta_dend = beta_dend, weights = weights, x_miss = x_miss, 
        reg_miss = reg_miss, x_in_row_is_missing = x_in_row_is_missing, nDelta_r = nDelta_r, nAlpha_r = nAlpha_r, d_design = d_design, 
        a_design = a_design)
  } else {    # D > 1 & h2_dims > 0
    standata = list(N = N, J = J, K = K, Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
        D = D, nDelta = nDelta, L = L, has_treg = has_treg, alpha_dstart = alpha_dstart, alpha_dend = alpha_dend, lambda_ind = lambda_ind,
        beta_dstart = beta_dstart, beta_dend = beta_dend, weights = weights, x_miss = x_miss, reg_miss = reg_miss, 
        x_in_row_is_missing = x_in_row_is_missing, nDelta_r = nDelta_r, nAlpha_r = nAlpha_r, d_design = d_design, a_design = a_design)
  }

  # Select the correct inirt implementation based on input parameters
  if(D == 1) {
    modl = stanmodels$inirt_unidim
  } else if(D > 1 & h2_dims == 0) {
    modl = stanmodels$inirt_mirt
  } else {    # D > 1 & h2_dims > 0
    modl = stanmodels$inirt_hoirt
  }

  # Choose a model estimateion method: variational inference or HMC
  if(method[1] == "vb") {
    out = rstan::vb(modl, data = standata, ...)
  } else if(method[1] == "hmc") {
    out = rstan::sampling(modl, data = standata, ...)
  } else {
    stop("Method not yet implemeted.")
  }

  return(out)
}