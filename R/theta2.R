#' Item Response Theory
#'
#' @export
#' @param data data frame
#' @param model text specification of the model. Not yet implemented. Instead use the other arguments
#' @param predictors list of predictors for each dimension. Each element of the list gives the predictors for the given dimension.
#' @param predictors_ranef list of random effect parameters for each dimension.
#' @param ranef_id id for which items belong to the current random effect.
#' @param predictors_ranef_corr list of correlated random effect parameters for each dimension.
#' @param 
#' @param dims dimensions
#' @param method Choose estimation method. Choices are method = 'vb' (default) for variational Bayes, or method = 'hmc' for Hamiltonian Monte Carlo.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' 
#' @export 
theta2 = function(data, model, exploratory = FALSE, method = c("vb", "hmc"), 
  weights = NULL, ...) {

  # item_id = NULL
  # predictors = NULL
  # predictors_ranef = NULL
  # ranef_id = NULL
  # predictors_ranef_corr = NULL
  # n_pranef_cor = NULL
  # dims = 1
  # h2_dims = 0
  # h2_dim_id = NULL
  
  # parse input model
  model_data = parse_model(model = model, data = data, exploratory = exploratory)
  dim_model_data = length(model_data)
  names_model_data = unlist(purrr::map(1:dim_model_data, \(x) {model_data[[x]]$type}))

  # next, translate input into each of the stan input variables

  # What type of IRT model are we fitting? - parsed, translated to input parameters
  irt_model = model_data[[1]]
  item_id = irt_model$item_id
  dims = irt_model$dims
  dim_id = irt_model$dim_id
  h2_dims = irt_model$h2_dims
  h2_dim_id = irt_model$h2_dim_id

  model_data = model_data[-1]

  # What type of IRT regression are we fitting? - parsed, translated to input parameters
  regr_alpha_delta = str_detect(names_model_data, "alpha|delta")
  regr_theta = model_data[[which(!regr_alpha_delta)]]
  regr_alpha_delta = model_data[-which(!regr_alpha_delta)]
  regr_alpha_data = regr_alpha_delta[str_detect(names_model_data[-1], "alpha")][[1]]
  regr_delta_data = regr_alpha_delta[str_detect(names_model_data[-1], "delta")][[1]]

  # For a single theta
  predictors = regr_theta$predictors
  if(!is.null(predictors)) {
    predictors = list(regr_theta$predictors)
  }
  predictors_ranef = regr_theta$predictors_ranef
  ranef_id = regr_theta$ranef_id
  predictors_ranef_corr = regr_theta$predictors_ranef_cor
  n_pranef_cor = regr_theta$n_pranef_cor

  irt_data = as.matrix(data[, item_id, drop = FALSE])
  N = nrow(irt_data)
  J = ncol(irt_data)
  N_long = N*J

  # alpha and delta
  a_design = matrix(1, nrow = N, ncol = 1)
  nAlpha_r = 1
  d_design = matrix(1, nrow = N, ncol = 1)
  nDelta_r = 1

  if(length(regr_alpha_data$predictors)) {
    a_design = cbind(a_design, data[, regr_alpha_data$predictors])
    nAlpha_r = ncol(a_design)
  }

  if(length(regr_delta_data$predictors)) {
    d_design = cbind(d_design, data[, regr_delta_data$predictors])
    nDelta_r = ncol(d_design)
  }

  # # Set up the structural regressions (item regression)
  # if(!is.null(structural_design)) {
  #   if(any(c("alpha") %in% names(structural_design))) {
  #     a_design = as.matrix(structural_design[["alpha"]])
  #     nAlpha_r = ncol(a_design)
  #   }
  #   if(any(c("delta") %in% names(structural_design))) {
  #     d_design = as.matrix(structural_design[["delta"]])
  #     nDelta_r = ncol(d_design)
  #   }
  # } else {
  #   a_design = matrix(1, nrow = N, ncol = 1)
  #   nAlpha_r = 1
  #   d_design = matrix(1, nrow = N, ncol = 1)
  #   nDelta_r = 1
  # }

  structural_design_ranef = 
    list(a_predictors = NULL, 
         a_predictors_ranef = NULL, 
         a_ranef_id = NULL, 
         a_predictors_ranef_corr = NULL, 
         a_n_pranef_cor = NULL,
         d_predictors = NULL, 
         d_predictors_ranef = NULL, 
         d_ranef_id = NULL, 
         d_predictors_ranef_corr = NULL, 
         d_n_pranef_cor = NULL)

  # if(!is.null(regr_alpha_data$predictors_ranef)) {
  structural_design_ranef$a_predictors_ranef = regr_alpha_data$new_reg_data[, regr_alpha_data$predictors_ranef]
  structural_design_ranef$a_ranef_id = regr_alpha_data$ranef_id

  structural_design_ranef$d_predictors_ranef = regr_delta_data$new_reg_data[, regr_alpha_data$predictors_ranef]
  structural_design_ranef$d_ranef_id = regr_delta_data$ranef_id
  # }

  structural_design_ranef$a_predictors_ranef_corr = regr_alpha_data$new_reg_data[, regr_alpha_data$predictors_ranef_cor]
  structural_design_ranef$a_n_pranef_cor = regr_alpha_data$n_pranef_cor

  structural_design_ranef$d_predictors_ranef_corr = regr_delta_data$new_reg_data[, regr_delta_data$predictors_ranef_cor]
  structural_design_ranef$d_n_pranef_cor = regr_delta_data$n_pranef_cor

  # Model input is parsed and converted into numeric variables. The rest
  # if inirt::inirt sets up the model for stan and then calls the pre-compiled
  # stan programs

  if(any(is.na(irt_data))) {
    model_missing_y = 1
  } else {
    model_missing_y = 0
  }
  
  has_treg = 0
  K = 0

  x = array(0, dim = c(N, 0))
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
    x = reg_data
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
    z = reg_data_ranef
    if(any(is.na(reg_data_ranef))) {
        model_missing_x_rand = 1
    }
  }

  u_Lzeta_cor = 0
  l_Lzeta_cor = 0
  Lzeta_cor = 0
  cor_z_item_ind = array(0, dim = c(0))
  cor_z_item_elem_ind = array(0, dim = c(0))
  z_c = array(0, dim = c(N, 0))

  if(!is.null(predictors_ranef_corr)) {
    has_treg = 1
    any_rand = 1
    any_rand_cor = 1
    predictors_ranef_corr_ulist = unlist(predictors_ranef_corr)
    reg_data_ranef_cor = data[, predictors_ranef_corr_ulist, drop = FALSE]
    Lzeta_cor = ncol(reg_data_ranef_cor)
    z_c = reg_data_ranef_cor
    l_Lzeta_cor = n_pranef_cor
    u_Lzeta_cor = Lzeta_cor / l_Lzeta_cor
    cor_z_item_ind = rep(1:u_Lzeta_cor, l_Lzeta_cor)
    cor_z_item_elem_ind = rep(1:l_Lzeta_cor, each = u_Lzeta_cor)
  }

  # ar = array(0, dim = c(N, 0))
  # dr = array(0, dim = c(N, 0))

  # uncorrelate random effect input quantities for the alpha, delta regression models
  alindex = array(0, dim = c(0))
  any_rand_ind_a = 0
  Laeta = 0
  Laeta_sd = 0
  aeta_sd_ind = array(0, dim = c(0))
  ar = array(0, dim = c(N, Laeta))
  
  dlindex = array(0, dim = c(0))
  any_rand_ind_d = 0
  Ldeta = 0
  Ldeta_sd = 0
  deta_sd_ind = array(0, dim = c(0))
  dr = array(0, dim = c(N, Ldeta))

  # uncorrelated random effects on alpha, delta parameters
  if(dim(structural_design_ranef$a_predictors_ranef)[2] | dim(structural_design_ranef$d_predictors_ranef)[2]) {
    
    if(dim(structural_design_ranef$a_predictors_ranef)[2]) {

      any_rand_ind_a = 1
      ar = structural_design_ranef$a_predictors_ranef
      Laeta = ncol(ar)
      a_ranef_id = structural_design_ranef$a_ranef_id
      Laeta_sd = length(unique(a_ranef_id))
      aeta_sd_ind = a_ranef_id
      alindex = 1:Laeta

    }

    if(dim(structural_design_ranef$d_predictors_ranef)[2]) {
      
      any_rand_ind_d = 1
      dr = structural_design_ranef$d_predictors_ranef
      Ldeta = ncol(dr)
      d_ranef_id = structural_design_ranef$d_ranef_id
      Ldeta_sd = length(unique(d_ranef_id))
      deta_sd_ind = d_ranef_id
      dlindex = 1:Ldeta

    }
    
  }

  # alpha set-up
  any_rand_cor_a = 0
  a_c = array(0, dim = c(N, 0))
  cor_a_item_ind = array(0, dim = c(0))
  cor_a_item_elem_ind = array(0, dim = c(0))
  Laeta_cor = 0
  u_Laeta_cor = 0
  l_Laeta_cor = 0

  # delta set-up
  any_rand_cor_d = 0
  d_c = array(0, dim = c(N, 0))
  cor_d_item_ind = array(0, dim = c(0))
  cor_d_item_elem_ind = array(0, dim = c(0))
  Ldeta_cor = 0
  u_Ldeta_cor = 0
  l_Ldeta_cor = 0

  if(dim(structural_design_ranef$a_predictors_ranef_corr)[2] | dim(structural_design_ranef$d_predictors_ranef)[2]) {

    if(dim(structural_design_ranef$a_predictors_ranef_corr)[2]) {
      any_rand = 1
      any_rand_cor_a = 1
      a_c = structural_design_ranef$a_predictors_ranef_corr
      Laeta_cor = ncol(structural_design_ranef$a_predictors_ranef_corr)
      l_Laeta_cor = structural_design_ranef$a_n_pranef_cor
      u_Laeta_cor = Laeta_cor / l_Laeta_cor
      cor_a_item_ind = rep(1:u_Laeta_cor, l_Laeta_cor)
      cor_a_item_elem_ind = rep(1:l_Laeta_cor, each = u_Laeta_cor)
    }

    if(dim(structural_design_ranef$d_predictors_ranef)[2]) {
      any_rand = 1
      any_rand_cor_d = 1
      d_c = structural_design_ranef$d_predictors_ranef_corr
      Ldeta_cor = ncol(structural_design_ranef$d_predictors_ranef_corr)
      l_Ldeta_cor = structural_design_ranef$d_n_pranef_cor
      u_Ldeta_cor = Ldeta_cor / l_Ldeta_cor
      cor_d_item_ind = rep(1:u_Ldeta_cor, l_Ldeta_cor)
      cor_d_item_elem_ind = rep(1:l_Ldeta_cor, each = u_Ldeta_cor)
    }

  }

  if(model_missing_y == 0) {
    Ncateg_max = max(irt_data)
    Ncategi = as.integer(apply(irt_data, 2, max))
    names(Ncategi) = NULL
  } else {
    Ncateg_max = max(as.vector(irt_data)[!is.na(as.vector(irt_data))])
    Ncategi = as.integer(apply(irt_data, 2, max)) # neex to update this to account for missingness!
    names(Ncategi) = NULL
  }
  
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
  x_miss = array(0, dim = c(N_long, 0))
  # reg_miss = 0
  # x_miss = 0
  # x_in_row_is_missing = 0
  beta_dstart = array(0, dim = c(0))
  beta_dend = array(0, dim = c(0))

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

  # Independent random effects in the theta regression (theta ~ rand1 + rand2 + ....)
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
    Lzeta = 0
    Lzeta_sd = 0
    zeta_sd_ind = array(0, dim = c(0))
    z = matrix(0, nrow = N, ncol = Lzeta)
    zeta_dstart = array(0, dim = c(0))
    zeta_dend = array(0, dim = c(0))
  }
  
  if(!is.null(predictors_ranef_corr)) {
    regress = 1
  }

  # sample weights
  if(is.null(weights)) {
    if(model_missing_y == 1) {
      weights = rep(1, N_long_obs)
    } else {
      weights = rep(1, N_long)
    }
  }

  # missing values in the item observation matrix
  N_miss = 0
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

  if(D == 1) {
    standata = list(N = N, J = J, K = K, any_rand = any_rand, any_rand_ind = any_rand_ind, any_rand_cor = any_rand_cor, 
        any_rand_ind_a = any_rand_ind_a, any_rand_cor_a = any_rand_cor_a, any_rand_ind_d = any_rand_ind_d, any_rand_cor_d = any_rand_cor_d,
        Ncateg_max = Ncateg_max, Ncategi = Ncategi, N_long = N_long, nn = nn, jj = jj, y = y, x = x, 
        D = D, nDelta = nDelta, L = L, has_treg = has_treg, beta_dstart = beta_dstart, beta_dend = beta_dend, zeta_dstart = zeta_dstart, zeta_dend = zeta_dend, 
        weights = weights, x_miss = x_miss, nDelta_r = nDelta_r, nAlpha_r = nAlpha_r,
        d_design = d_design, a_design = a_design, Lzeta = Lzeta, Laeta = Laeta, Ldeta = Ldeta, 
        u_Lzeta_cor = u_Lzeta_cor, l_Lzeta_cor = l_Lzeta_cor, 
        u_Laeta_cor = u_Laeta_cor, l_Laeta_cor = l_Laeta_cor, 
        u_Ldeta_cor = u_Ldeta_cor, l_Ldeta_cor = l_Ldeta_cor, 
        Lzeta_cor = Lzeta_cor, 
        Laeta_cor = Laeta_cor, 
        Ldeta_cor = Ldeta_cor, 
        z = z, 
        ar = ar,
        dr = dr,
        Lzeta_sd = Lzeta_sd, zeta_sd_ind = zeta_sd_ind, cor_z_item_ind = cor_z_item_ind, cor_z_item_elem_ind = cor_z_item_elem_ind, 
        Laeta_sd = Laeta_sd, alindex = alindex, aeta_sd_ind = aeta_sd_ind, cor_a_item_ind = cor_a_item_ind, cor_a_item_elem_ind = cor_a_item_elem_ind,
        Ldeta_sd = Ldeta_sd, dlindex = dlindex, deta_sd_ind = deta_sd_ind, cor_d_item_ind = cor_d_item_ind, cor_d_item_elem_ind = cor_d_item_elem_ind,
        z_c = z_c,
        a_c = a_c, 
        d_c = d_c)
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
