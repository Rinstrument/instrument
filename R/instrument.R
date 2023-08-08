#' Bayesian Explanatory Multidimensional Item Response Theory
#' 
#' Fit a multidimensional item response theory model using Stan.
#' See github.com/theta2pack/doc for examples, tutorials, and documentation.
#'
#' @param data a named `data.frame` or `matrix`
#' @param model text specification of the model. See examples or website.
#' @param itype item type of the items. Options are '1pl', '2pl', '3pl' 
#' @param ranef_id id for which items belong to the current random effect.
#' @param exploratory fit exploratory MIRT model? Only used with multidimensional
#' model that is not higher-order or bifactor
#' @param method Choose estimation method. Use method = 'hmc' for Hamiltonian 
#' Monte Carlo. Other option is method = 'vb' for Variational Bayes. HMC should 
#' be the default for published results. VB is experimental.
#' @param fweights frequency weights. These weights represent how many times
#' each row of the data should be repeated in the analysis. fweights weight
#' the likelihood function.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' 
#' @export 
instrument = function(data, model, itype, exploratory = FALSE, method = c("vb", "hmc"), 
  fweights = NULL, ...) {

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

  # Given IRT model, what item types are being fitted?
  if(length(itype) == 1) {
    itype = rep(which(itype == c("1pl", "2pl", "3pl")), length(item_id))
  } else {
    itype = match(itype, c("1pl", "2pl", "3pl"))
  }

  # Check for 3parameter item types
  all_3pl = itype == 3
  any_eta3pl = any(all_3pl) * 1
  nEta3pl = sum(all_3pl)

  # What type of IRT regression are we fitting? - parsed, translated to input parameters
  regr_alpha_delta = stringr::str_detect(names_model_data, "alpha|delta")
  regr_theta = model_data[which(!regr_alpha_delta)] #[[]] ???
  names_regr_theta = unlist(sapply(regr_theta, \(x) {x["type"]}))
  regr_alpha_delta = model_data[which(regr_alpha_delta)]
  names_model_data = names_model_data[!(names_model_data %in% names_regr_theta)] # changed this line from "theta" to generic name
  regr_alpha_data = regr_alpha_delta[stringr::str_detect(names_model_data, "alpha")][[1]] # [-1] ???
  regr_delta_data = regr_alpha_delta[stringr::str_detect(names_model_data, "delta")][[1]]

  predictors = NULL

  predictors = lapply(regr_theta, \(x) {x$predictors})

  # does this need to be protected agains empty names_regr_theta
  find_dims = which(irt_model[['dim_names']] %in% names_regr_theta)

  # predictors = regr_theta$predictors
  #   if(!is.null(predictors)) {
  #     predictors = list(regr_theta$predictors)
  #   }

  if(all(sapply(regr_theta, \(x) { is.null(x$predictors_ranef_cor) }))) {
    # base case
    which_dim_cor_reg = rep(0, dims)

  } else if(length(regr_theta) == 1) { 
    # For a single theta
    which_dim_cor_reg = c(which(names_regr_theta == irt_model$dim_names), rep(0, dims - 1))

  } else {
    # all or multiple thetas?
    which_dim_cor_reg = c(which(names_regr_theta == irt_model$dim_names), rep(0, dims - length(names_regr_theta)))

  }

  which_dim_ind_reg = array(0, dim = dims)
  which_dim_ind_reg_sort = array(0, dim = dims)

  if(any(sapply(regr_theta, \(x) { !is.null(x[['predictors_ranef']]) }))) {

    which_dim_ind_reg = fill_match(find_dims, dims)

    which_dim_ind_reg = c(
        which_dim_ind_reg[which_dim_ind_reg > 0],
        which_dim_ind_reg[which_dim_ind_reg == 0]
      )
    which_dim_ind_reg_sort = sort(which_dim_ind_reg, TRUE)
    # which_dim_ind_reg = which_dim_ind_reg_sort
    which_dim_ind_reg_sort = c(
        (1:dims)[which_dim_ind_reg_sort > 0],
        (rep(0, dims))[which_dim_ind_reg_sort == 0]
      )
    
  }

  

  # if(all(sapply(regr_theta, \(x) { is.null(x$predictors_ranef) }))) {
  #   # base case
  #   which_dim_ind_reg = rep(0, dims)

  # } else if(length(regr_theta) == 1) { 
  #   # For a single theta
  #   which_dim_ind_reg = c(which(names_regr_theta == irt_model$dim_names), rep(0, dims - 1))

    


  # } else {
  #   # all or multiple thetas?
  #   which_dim_ind_reg = c(which(names_regr_theta == irt_model$dim_names), rep(0, dims - length(names_regr_theta)))

  # }
  
  predictors_ranef = lapply(regr_theta, \(x) {x$predictors_ranef})
  ranef_id = lapply(regr_theta, \(x) {x$ranef_id})
  predictors_ranef_corr = lapply(regr_theta, \(x) {x$predictors_ranef_cor})
  n_pranef_cor = lapply(regr_theta, \(x) {x$n_pranef_cor})

  irt_data = as.matrix(data[, item_id, drop = FALSE])
  irt_data = collapse_categories(irt_data)  # collapse missing categories in IRT data set
  N = nrow(irt_data)
  J = ncol(irt_data)
  N_long = N*J

  y = as.vector(irt_data)
  
  # alpha and delta
  if(any(itype >= 2)) {
    if(is.null(regr_alpha_data$predictors)) {
      a_design = matrix(1, nrow = N, ncol = 0)
      nAlpha_r = 0
      LMean = 0
      DAlpha = 1
    } else {
      a_design = matrix(1, nrow = N, ncol = 1)
      nAlpha_r = 1
      LMean = 1
      DAlpha = 1
    }
  } else {
    a_design = matrix(1, nrow = N, ncol = 0)
    nAlpha_r = 0
    LMean = 0
    DAlpha = 0
  }
  
  if(is.null(regr_delta_data$predictors)) {
    d_design = matrix(1, nrow = N, ncol = 0)
    deltaMean = 0
    nDelta_r = 0
  } else {
    d_design = matrix(1, nrow = N, ncol = 1)
    deltaMean = 1
    nDelta_r = 1
  }
  
  if(length(regr_alpha_data$predictors)) {
    if(all(itype == 1)) {
      stop("Invalid model. Cannot specify an alpha ~ regression when itype == 1pl.
            Alpha is a parameter in the 2pl model.")
    }
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

  structural_design_ranef$d_predictors_ranef = regr_delta_data$new_reg_data[, regr_delta_data$predictors_ranef]
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

  l_y_nonMiss = length(y)
  y_nonMiss = !is.na(y)
  if(any(!y_nonMiss)) {
    l_y_nonMiss = sum(y_nonMiss)
  }
  
  has_treg = 0
  K = 0

  x = array(0, dim = c(N, 0))
  xLong = array(0, dim = c(N_long, 0))
  any_rand = 0
  any_rand_ind = 0
  any_rand_cor = 0
  model_missing_x = 0
  model_missing_x_rand = 0

  if(any(!is.null(unlist(predictors)))) {
    has_treg = 1
    predictor_ulist = unlist(predictors)
    reg_data = data[, predictor_ulist, drop = FALSE]
    K = ncol(reg_data)
    x = reg_data
    if(any(is.na(reg_data))) {
        model_missing_x = 1
    }
    xLong = matrix(rep(t(x), J), ncol = ncol(x), byrow = TRUE)
  }

  # extra memory slots used if more than one set of correlated random effects
  # is to be estimates
  extra_mem_slots = 32

  # if(length(regr_theta) > 1) {
  for(i in 2:extra_mem_slots) {
    eval(parse(text = paste0("rand_ind_g", i - 1, " = 0")))
    eval(parse(text = paste0("Lzeta_", i, " = 0")))
    eval(parse(text = paste0("Lzeta_sd_", i, " = 0")))
    eval(parse(text = paste0("zeta_sd_ind_", i, " = array(0, dim = 0)")))
    eval(parse(text = paste0("z_", i, " = array(0, dim = c(N, 0))")))
    eval(parse(text = paste0("zLong_", i, " = array(0, dim = c(l_y_nonMiss, 0))")))# ??? N_long
    eval(parse(text = paste0("zeta_sd_ind_diag_", i, " = array(0, dim = c(0, 0))")))
  }
  # }

  Lzeta = 0
  z = array(0, dim = c(N, 0))
  Lzeta_sd = 0
  zeta_sd_ind = array(0, dim = c(0))
  zeta_sd_ind_diag = array(0, dim = c(0, 0))
  zeta_sd_ind_ones = array(0, dim = c(0))
    # z = matrix(0, nrow = N, ncol = Lzeta)
  if(!is.null(unlist(predictors_ranef))) {
    # has_treg = 1
    any_rand = 1
    any_rand_ind = 1

    # how to distringuish
    reg_data_ranef = lapply(regr_theta, \(x) { x$new_reg_data })
    Lzeta_sd_list = lapply(regr_theta, \(x) { length(unique(x[['ranef_id']])) })
    zeta_sd_ind_list = lapply(regr_theta, \(x) { x[['ranef_id']] })

    Lzeta = ncol(reg_data_ranef[[1]])
    z = reg_data_ranef[[1]]
    Lzeta_sd = Lzeta_sd_list[[1]]
    zeta_sd_ind = zeta_sd_ind_list[[1]]

    # diagonalize
    zeta_sd_ind_diag = matrix(0, Lzeta, Lzeta_sd)
    for(ll in 1:Lzeta) {
      for(jl in 1:Lzeta_sd) {
        zeta_sd_ind_diag[ll, jl] = ifelse(zeta_sd_ind[ll] == jl, 1, 0)
      }
    }
    
    if(length(regr_theta) > 1) {
      for(i in 2:length(regr_theta)) {
        eval(parse(text = paste0("rand_ind_g", i - 1, " = 1")))
        eval(parse(text = paste0("Lzeta_", i, " = ncol(reg_data_ranef[[", i, "]])")))
        eval(parse(text = paste0("Lzeta_sd_", i, " = Lzeta_sd_list[[", i, "]]")))
        eval(parse(text = paste0("zeta_sd_ind_", i, " = zeta_sd_ind_list[[", i, "]]")))
        eval(parse(text = paste0("z_", i, " = reg_data_ranef[[", i, "]]")))
        eval(parse(text = paste0("zLong_", i, " = array(0, dim = c(l_y_nonMiss, 0))")))
        eval(parse(text = paste0("zeta_sd_ind_diag_", i, " = array(0, dim = c(Lzeta_", i, ", Lzeta_sd_", i, "))")))
        eval(
          parse(
            text = paste0(
              "for(ll in 1:Lzeta_", i, ") {
                for(jl in 1:Lzeta_sd_", i, ") {
                  zeta_sd_ind_diag_", i, "[ll, jl] = ifelse(zeta_sd_ind_", i, "[ll] == jl, 1, 0)
                }
              }"
            )
          )
        )
      }
    }
    
  }

  u_Lzeta_cor = 0 # number of correlated random effect params
  l_Lzeta_cor = 0 # length of random effect parameter vector

  Lzeta_cor = 0 # total number of first set of correlated random effects

  cor_z_item_ind = array(0, dim = c(0)) # index random effect vector
  cor_z_item_elem_ind = array(0, dim = c(0)) # index the position within random effect vectors
  z_c = array(0, dim = c(N, 0)) # random effect design matrix

 

  # if(length(regr_theta) > 1) {
  for(i in 2:extra_mem_slots) {
    eval(parse(text = paste0("rand_cor_g", i - 1, " = 0")))
    eval(parse(text = paste0("u_Lzeta_cor_", i, " = 0")))
    eval(parse(text = paste0("l_Lzeta_cor_", i, " = 0")))
    eval(parse(text = paste0("Lzeta_cor_", i, " = 0")))
    eval(parse(text = paste0("cor_z_item_ind_", i, " = array(0, dim = c(0))")))
    eval(parse(text = paste0("cor_z_item_elem_ind_", i, " = array(0, dim = c(0))")))
    eval(parse(text = paste0("z_c_", i, " = array(0, dim = c(N, 0))")))
    eval(parse(text = paste0("z_cLong_", i, " = array(0, dim = c(l_y_nonMiss, 0))")))
  }
  # }
  

  if(!is.null(unlist(predictors_ranef_corr))) {
    # has_treg = 1
    any_rand = 1
    any_rand_cor = 1
    # predictors_ranef_corr_ulist = unlist(predictors_ranef_corr)
    # reg_data_ranef_cor = data[, predictors_ranef_corr_ulist, drop = FALSE]

    # reg_data_ranef_cor = lapply(predictors_ranef_corr, \(x, d) {d[, x, drop = FALSE]}, d = data)

    reg_data_ranef_cor = lapply(regr_theta, \(x) {x$new_reg_data})

    Lzeta_cor = ncol(reg_data_ranef_cor[[1]])
    z_c = reg_data_ranef_cor[[1]]
    l_Lzeta_cor = n_pranef_cor[[1]]
    u_Lzeta_cor = Lzeta_cor / l_Lzeta_cor
    cor_z_item_ind = rep(1:u_Lzeta_cor, l_Lzeta_cor)
    cor_z_item_elem_ind = rep(1:l_Lzeta_cor, each = u_Lzeta_cor)

    if(length(regr_theta) > 1) {
      for(i in 2:length(regr_theta)) {
        eval(parse(text = paste0("rand_cor_g", i - 1, " = 1")))
        eval(parse(text = paste0("Lzeta_cor_", i, paste0(" = ncol(reg_data_ranef_cor[[", i, "]])"))))
        eval(parse(text = paste0("z_c_", i, paste0(" = reg_data_ranef_cor[[", i, "]]"))))
        eval(parse(text = paste0("l_Lzeta_cor_", i, paste0(" = n_pranef_cor[[", i, "]]"))))
        eval(parse(text = paste0("u_Lzeta_cor_", i, paste0(" = Lzeta_cor_", i, " / ", "l_Lzeta_cor_", i))))
        eval(parse(text = paste0("cor_z_item_ind_", i, paste0(" = rep(1:u_Lzeta_cor_", i, ", l_Lzeta_cor_", i, ")"))))
        eval(parse(text = paste0("cor_z_item_elem_ind_", i, paste0(" = rep(1:l_Lzeta_cor_", i, ", each = u_Lzeta_cor_", i, ")"))))
      }
    }

    # long version of z_c, z_cLong
    z_cLong = matrix(rep(t(z_c), J), 
                     ncol = ncol(z_c), 
                     byrow = TRUE
                     )[y_nonMiss, , drop = FALSE]
    if(length(regr_theta) > 1) {
      for(i in 2:length(regr_theta)) {
        eval(parse(text = paste0("z_cLong_", i, " = matrix(rep(t(z_c_", i, "), J), 
                                                           ncol = ncol(z_c_", i, "), 
                                                           byrow = TRUE
                                                           )[y_nonMiss, , drop = FALSE]")))
      }
    }

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

      if(all(itype == 1)) {
      stop("Invalid model. Cannot specify an alpha ~ regression when itype == 1pl.
            Alpha is a parameter in the 2pl model.")
      }

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

      if(all(itype == 1)) {
      stop("Invalid model. Cannot specify an alpha ~ regression when itype == 1pl.
            Alpha is a parameter in the 2pl model.")
      }

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

  # if(model_missing_y == 0) {
  Ncategi = as.integer(apply(irt_data, 2, \(x) {max(x, na.rm = TRUE)}))
  Ncateg_max = max(Ncategi)
  names(Ncategi) = NULL
  # } else {
    # Ncategi = as.integer(apply(irt_data, 2, \(x) {max(x, na.rm = TRUE)}))
    # Ncateg_max = max(Ncategi)
    # names(Ncategi) = NULL
  # }
  
  Ncategi_jj = rep(Ncategi, each = N)
  nn = rep(1:N, J)
  jj = rep(1:J, each = N)

  # missing values in the item observation matrix
  N_miss = 0
  if(model_missing_y == 1) {
    Ncategi_jj = Ncategi_jj[y_nonMiss]
    N_miss = sum(is.na(irt_data))
    N_long = N_long - N_miss
    nn = nn[y_nonMiss]
    jj = jj[y_nonMiss]
    y = y[y_nonMiss]
    xLong = matrix(rep(t(x), J), 
                   ncol = ncol(x), 
                   byrow = TRUE
                   )[y_nonMiss, , drop = FALSE]

    if(!is.null(unlist(predictors_ranef))) {
      zLong = matrix(rep(t(z), J), 
                     ncol = ncol(z), 
                     byrow = TRUE
                     )[y_nonMiss, , drop = FALSE]
      for(i in 2:length(regr_theta)) {
        eval(parse(text = paste0("zLong_", i, " = matrix(rep(t(z_", i, "), J), 
                                                         ncol = ncol(z_", i, "), 
                                                         byrow = TRUE
                                                         )[y_nonMiss, , drop = FALSE]")))
      }
    } else { # this could probably be moved to the other pre-instantiations up above
      zLong = array(0, dim = c(N_long, Lzeta))
    }

  }

  D = dims
  match_eta3pl = match(itype, 3, nomatch = 0)
  find_eta3pl = rep(cumsum(match_eta3pl) * match_eta3pl, each = N)

  if(model_missing_y == 1) {
    find_eta3pl = find_eta3pl[1:N_long]
  }

  itype = rep(itype, each = N)

  if(model_missing_y == 1) {
    itype = itype[1:N_long]
  }
  
  if(h2_dims == 0) {
    if(all(itype >= 2)) {
      L = D*(J-D) + D*(D+1)/2
    } else if(all(itype == 1)) {
      L = 0
    }
  } else if(h2_dims == 1){
    if(all(itype == 2)) {
      L = J
    } else if(all(itype == 1)) {
      L = 0
    }
  } else {
    stop("h2_dims > 1 not yet implemented. Choose 0 or 1.")
  }

  nDelta = sum(Ncategi - 1)

  regress = 0
  x_miss = array(0, dim = c(N, 0))
  # reg_miss = 0
  # x_miss = 0
  # x_in_row_is_missing = 0
  beta_dstart = array(0, dim = c(0))
  beta_dend = array(0, dim = c(0))

  which_dim_fixed_reg = array(0, dim = dims)
  which_dim_fixed_reg_sort = array(0, dim = dims)

  if(any(!is.null(unlist(predictors)))) {
    regress = 1
    start_index = 1
    if(h2_dims > 0) {
      beta_dstart = numeric(1)
      beta_dend = numeric(1)
      beta_dstart[1] = start_index
      beta_dend[1] = start_index + length(predictors[[1]]) - 1
      start_index = start_index + length(predictors[[1]])
      beta_dstart = array(beta_dstart, dim = 1)
      beta_dend = array(beta_dend, dim = 1)
    }

    if(h2_dims == 0) {
      # lapply(predictors, \(x) {is.null(x)})
      # sapply(regr_theta, \(x) {x[['type']]}) ### here ###
      # irt_model$dim_names
      # names_regr_theta

      which_dim_fixed_reg = fill_match(find_dims, dims)

      which_dim_fixed_reg = c(
          which_dim_fixed_reg[which_dim_fixed_reg > 0],
          which_dim_fixed_reg[which_dim_fixed_reg == 0]
        )
      which_dim_fixed_reg_sort = sort(which_dim_fixed_reg, TRUE)
      # which_dim_ind_reg = which_dim_ind_reg_sort
      which_dim_fixed_reg_sort = c(
          (1:dims)[which_dim_fixed_reg_sort > 0],
          (rep(0, dims))[which_dim_fixed_reg_sort == 0]
        )

      predictors_by_dim = irt_model$dim_names %in% names_regr_theta
      start_index = 1
      beta_dstart = numeric(D) + 1
      beta_dend = numeric(D)
      p_ind = 1
      for(d in 1:D) {
        if(which_dim_fixed_reg[d] > 0) {
          beta_dstart[d] = start_index
          beta_dend[d] = start_index + length(predictors[[p_ind]]) - 1
          start_index = start_index + length(predictors[[p_ind]])
          p_ind = p_ind + 1
        }
      }
      beta_dstart = array(beta_dstart, dim = D)
      beta_dend = array(beta_dend, dim = D)
    }

    if(model_missing_x == 0) {
      reg_miss = is.na(reg_data) * 1
      x_miss = reg_miss
      # x_miss = matrix(rep(t(reg_miss), J), ncol = K, byrow = TRUE)
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
  if(!is.null(unlist(predictors_ranef))) {
    regress = 1
    start_index = 1
    zeta_dstart = rep(1, D)
    zeta_dend = numeric(D)
    p_ind = 1
    for(d in 1:D) { # insert which_ thing here (equivalent to correlated case)
      if(which_dim_ind_reg[d] > 0) {
        zeta_dstart[d] = 1
        zeta_dend[d] = 1 + length(predictors_ranef[[p_ind]]) - 1
        # start_index = start_index + length(predictors_ranef[[p_ind]])
        p_ind = p_ind + 1
      }
    }
    zeta_dstart = array(zeta_dstart, dim = D)
    zeta_dend = array(zeta_dend, dim = D)
    # Lzeta_sd = length(unique(ranef_id))
    # zeta_sd_ind = ranef_id
  } else {
    Lzeta = 0
    # Lzeta_sd = 0
    # zeta_sd_ind = array(0, dim = c(0))
    # z = matrix(0, nrow = N, ncol = Lzeta)
    zeta_dstart = array(0, dim = c(0))
    zeta_dend = array(0, dim = c(0))
  }
  
  if(!is.null(unlist(predictors_ranef_corr))) {
    regress = 1
  }

  # # missing values in the item observation matrix
  # N_miss = 0
  # if(model_missing_y == 1) {
  #   N_miss = sum(is.na(irt_data))
  #   N_long = N_long - N_miss
  #   nn = nn[!is.na(y)]
  #   jj = jj[!is.na(y)]
  #   y = y[!is.na(y)]
  # }

  # frequency weights
  if(is.null(fweights)) {
    if(model_missing_y == 1) {
      fweights = rep(1, N_long)
    } else {
      fweights = rep(1, N_long)
    }
  }

  if(h2_dims > 0) {
    # h2_dim_id_ulist = unlist(h2_dim_id)
    # d_lengths = sapply(h2_dim_id, function(x) {length(x)})
    d_lengths = as.vector(table(h2_dim_id))
    d_seq_start = 1
    new_id_list = vector(mode = "list", length = D)
    for(i in 1:D) {
      new_id_list[[i]] = d_seq_start:(d_lengths[i] + d_seq_start - 1)
      d_seq_start = max(new_id_list[[i]]) + 1
    }
    new_id_list = lapply(new_id_list, function(x) {c(x[1], x[length(x)])})
    # irt_data = irt_data[, h2_dim_id_ulist]
    alpha_dstart = array(sapply(new_id_list, function(x) {x[1]}), dim = D)
    alpha_dend = array(sapply(new_id_list, function(x) {x[length(x)]}), dim = D)
    # lambda_ind = c()
    # for(i in 1:D) {
    #   lambda_ind = c(lambda_ind, rep(i, d_lengths[i]))
    # }
    # lambda_ind = as.vector(t(replicate(N, lambda_ind)))
    lambda_ind = rep(h2_dim_id, each = N)
  } else if(h2_dims == 0 & exploratory == TRUE){
    alpha_dstart = array(rep(1, D), dim = D)
    alpha_dend = array(rep(max(item_id), D), dim = D)
  } else if(h2_dims == 0 & exploratory == FALSE) {
    # update alpha_dstart and alpha_dstart
  }

  if(D == 1) {

    standata = list(
        N = N, 
        J = J, 
        K = K, 
        any_rand = any_rand, 
        any_rand_ind = any_rand_ind, 
        any_rand_cor = any_rand_cor, 
        any_rand_ind_a = any_rand_ind_a, 
        any_rand_cor_a = any_rand_cor_a, 
        any_rand_ind_d = any_rand_ind_d, 
        any_rand_cor_d = any_rand_cor_d,
        Ncateg_max = Ncateg_max, 
        Ncategi = Ncategi, 
        N_long = N_long, 
        nn = nn, 
        jj = jj, 
        y = y, 
        itype = itype, 
        any_eta3pl = any_eta3pl, 
        nEta3pl = nEta3pl, 
        find_eta3pl = find_eta3pl, 
        x = x, 
        D = D, 
        DAlpha = DAlpha, 
        nDelta = nDelta, 
        L = L, 
        has_treg = has_treg, 
        beta_dstart = beta_dstart, 
        beta_dend = beta_dend, 
        zeta_dstart = zeta_dstart, 
        zeta_dend = zeta_dend, 
        fweights = fweights, 
        # x_miss = x_miss, 
        nDelta_r = nDelta_r, 
        nAlpha_r = nAlpha_r, 
        LMean = LMean, 
        deltaMean = deltaMean, 
        d_design = d_design, 
        a_design = a_design, 
        Lzeta = Lzeta, 
        Laeta = Laeta, 
        Ldeta = Ldeta, 
        u_Lzeta_cor = u_Lzeta_cor, 
        l_Lzeta_cor = l_Lzeta_cor, 
        u_Laeta_cor = u_Laeta_cor, 
        l_Laeta_cor = l_Laeta_cor, 
        u_Ldeta_cor = u_Ldeta_cor, 
        l_Ldeta_cor = l_Ldeta_cor, 
        Lzeta_cor = Lzeta_cor, 
        Laeta_cor = Laeta_cor, 
        Ldeta_cor = Ldeta_cor, 
        z = z, 
        ar = ar,
        dr = dr,
        Lzeta_sd = Lzeta_sd, 
        zeta_sd_ind = zeta_sd_ind, 
        cor_z_item_ind = cor_z_item_ind, 
        cor_z_item_elem_ind = cor_z_item_elem_ind, 
        Laeta_sd = Laeta_sd, 
        alindex = alindex, 
        aeta_sd_ind = aeta_sd_ind, 
        cor_a_item_ind = cor_a_item_ind, 
        cor_a_item_elem_ind = cor_a_item_elem_ind,
        Ldeta_sd = Ldeta_sd, 
        dlindex = dlindex, 
        deta_sd_ind = deta_sd_ind, 
        cor_d_item_ind = cor_d_item_ind, 
        cor_d_item_elem_ind = cor_d_item_elem_ind,
        z_c = z_c,
        a_c = a_c, 
        d_c = d_c
      )

  } else if(D > 1 & h2_dims == 0) {

         
        # Lzeta_2  = NULL
        # Lzeta_3  = NULL
        # Lzeta_4  = NULL
        # Lzeta_5  = NULL
        # Lzeta_6  = NULL
        # Lzeta_7  = NULL
        # Lzeta_8  = NULL
        # Lzeta_9  = NULL
        # Lzeta_10 = NULL
        # Lzeta_11 = NULL
        # Lzeta_12 = NULL
        # Lzeta_13 = NULL
        # Lzeta_14 = NULL
        # Lzeta_15 = NULL
        # Lzeta_16 = NULL
        # Lzeta_17 = NULL
        # Lzeta_18 = NULL
        # Lzeta_19 = NULL
        # Lzeta_20 = NULL
        # Lzeta_21 = NULL
        # Lzeta_22 = NULL
        # Lzeta_23 = NULL
        # Lzeta_24 = NULL
        # Lzeta_25 = NULL
        # Lzeta_26 = NULL
        # Lzeta_27 = NULL
        # Lzeta_28 = NULL
        # Lzeta_29 = NULL
        # Lzeta_30 = NULL
        # Lzeta_31 = NULL
        # Lzeta_32 = NULL
        # zLong_2  = NULL
        # zLong_3  = NULL
        # zLong_4  = NULL
        # zLong_5  = NULL
        # zLong_6  = NULL
        # zLong_7  = NULL
        # zLong_8  = NULL
        # zLong_9  = NULL
        # zLong_10 = NULL
        # zLong_11 = NULL
        # zLong_12 = NULL
        # zLong_13 = NULL
        # zLong_14 = NULL
        # zLong_15 = NULL
        # zLong_16 = NULL
        # zLong_17 = NULL
        # zLong_18 = NULL
        # zLong_19 = NULL
        # zLong_20 = NULL
        # zLong_21 = NULL
        # zLong_22 = NULL
        # zLong_23 = NULL
        # zLong_24 = NULL
        # zLong_25 = NULL
        # zLong_26 = NULL
        # zLong_27 = NULL
        # zLong_28 = NULL
        # zLong_29 = NULL
        # zLong_30 = NULL
        # zLong_31 = NULL
        # zLong_32 = NULL
        # Lzeta_sd_2  = NULL
        # Lzeta_sd_3  = NULL
        # Lzeta_sd_4  = NULL
        # Lzeta_sd_5  = NULL
        # Lzeta_sd_6  = NULL
        # Lzeta_sd_7  = NULL
        # Lzeta_sd_8  = NULL
        # Lzeta_sd_9  = NULL
        # Lzeta_sd_10 = NULL
        # Lzeta_sd_11 = NULL
        # Lzeta_sd_12 = NULL
        # Lzeta_sd_13 = NULL
        # Lzeta_sd_14 = NULL
        # Lzeta_sd_15 = NULL
        # Lzeta_sd_16 = NULL
        # Lzeta_sd_17 = NULL
        # Lzeta_sd_18 = NULL
        # Lzeta_sd_19 = NULL
        # Lzeta_sd_20 = NULL
        # Lzeta_sd_21 = NULL
        # Lzeta_sd_22 = NULL
        # Lzeta_sd_23 = NULL
        # Lzeta_sd_24 = NULL
        # Lzeta_sd_25 = NULL
        # Lzeta_sd_26 = NULL
        # Lzeta_sd_27 = NULL
        # Lzeta_sd_28 = NULL
        # Lzeta_sd_29 = NULL
        # Lzeta_sd_30 = NULL
        # Lzeta_sd_31 = NULL
        # Lzeta_sd_32 = NULL

        # zeta_sd_ind    = NULL
        # zeta_sd_ind_2  = NULL
        # zeta_sd_ind_3  = NULL
        # zeta_sd_ind_4  = NULL
        # zeta_sd_ind_5  = NULL
        # zeta_sd_ind_6  = NULL
        # zeta_sd_ind_7  = NULL
        # zeta_sd_ind_8  = NULL
        # zeta_sd_ind_9  = NULL
        # zeta_sd_ind_10 = NULL
        # zeta_sd_ind_11 = NULL
        # zeta_sd_ind_12 = NULL
        # zeta_sd_ind_13 = NULL
        # zeta_sd_ind_14 = NULL
        # zeta_sd_ind_15 = NULL
        # zeta_sd_ind_16 = NULL
        # zeta_sd_ind_17 = NULL
        # zeta_sd_ind_18 = NULL
        # zeta_sd_ind_19 = NULL
        # zeta_sd_ind_20 = NULL
        # zeta_sd_ind_21 = NULL
        # zeta_sd_ind_22 = NULL
        # zeta_sd_ind_23 = NULL
        # zeta_sd_ind_24 = NULL
        # zeta_sd_ind_25 = NULL
        # zeta_sd_ind_26 = NULL
        # zeta_sd_ind_27 = NULL
        # zeta_sd_ind_28 = NULL
        # zeta_sd_ind_29 = NULL
        # zeta_sd_ind_30 = NULL
        # zeta_sd_ind_31 = NULL
        # zeta_sd_ind_32 = NULL

        # zeta_sd_ind_diag    = NULL
        # zeta_sd_ind_diag_2  = NULL
        # zeta_sd_ind_diag_3  = NULL
        # zeta_sd_ind_diag_4  = NULL
        # zeta_sd_ind_diag_5  = NULL
        # zeta_sd_ind_diag_6  = NULL
        # zeta_sd_ind_diag_7  = NULL
        # zeta_sd_ind_diag_8  = NULL
        # zeta_sd_ind_diag_9  = NULL
        # zeta_sd_ind_diag_10 = NULL
        # zeta_sd_ind_diag_11 = NULL
        # zeta_sd_ind_diag_12 = NULL
        # zeta_sd_ind_diag_13 = NULL
        # zeta_sd_ind_diag_14 = NULL
        # zeta_sd_ind_diag_15 = NULL
        # zeta_sd_ind_diag_16 = NULL
        # zeta_sd_ind_diag_17 = NULL
        # zeta_sd_ind_diag_18 = NULL
        # zeta_sd_ind_diag_19 = NULL
        # zeta_sd_ind_diag_20 = NULL
        # zeta_sd_ind_diag_21 = NULL
        # zeta_sd_ind_diag_22 = NULL
        # zeta_sd_ind_diag_23 = NULL
        # zeta_sd_ind_diag_24 = NULL
        # zeta_sd_ind_diag_25 = NULL
        # zeta_sd_ind_diag_26 = NULL
        # zeta_sd_ind_diag_27 = NULL
        # zeta_sd_ind_diag_28 = NULL
        # zeta_sd_ind_diag_29 = NULL
        # zeta_sd_ind_diag_30 = NULL
        # zeta_sd_ind_diag_31 = NULL
        # zeta_sd_ind_diag_32 = NULL


    standata = list(
        N = N, 
        J = J, 
        K = K, 
        any_rand = any_rand, 
        any_rand_ind = any_rand_ind, 
        any_rand_cor = any_rand_cor, 
        any_rand_ind_a = any_rand_ind_a, 
        any_rand_cor_a = any_rand_cor_a, 
        any_rand_ind_d = any_rand_ind_d, 
        any_rand_cor_d = any_rand_cor_d,
        Ncateg_max = Ncateg_max, 
        Ncategi = Ncategi, 
        N_long = N_long, 
        nn = nn, 
        jj = jj, 
        y = y, 
        itype = itype, 
        any_eta3pl = any_eta3pl, 
        nEta3pl = nEta3pl, 
        find_eta3pl = find_eta3pl, 
        x = x, 
        xLong = xLong,
        D = D, 
        DAlpha = DAlpha, 
        alpha_dstart = alpha_dstart,
        alpha_dend = alpha_dend,
        nDelta = nDelta, 
        L = L, 
        has_treg = has_treg, 
        beta_dstart = beta_dstart, 
        beta_dend = beta_dend, 
        zeta_dstart = zeta_dstart, 
        zeta_dend = zeta_dend, 
        fweights = fweights, 
        # x_miss = x_miss, 
        nDelta_r = nDelta_r, 
        nAlpha_r = nAlpha_r, 
        LMean = LMean, 
        deltaMean = deltaMean, 
        d_design = d_design, 
        a_design = a_design, 
        which_dim_fixed_reg = which_dim_fixed_reg,
        which_dim_fixed_reg_sort = which_dim_fixed_reg_sort,
        which_dim_ind_reg = which_dim_ind_reg,
        which_dim_ind_reg_sort = which_dim_ind_reg_sort,

        Lzeta    = Lzeta, 
        Lzeta_2  = Lzeta_2,
        Lzeta_3  = Lzeta_3,
        Lzeta_4  = Lzeta_4,
        Lzeta_5  = Lzeta_5,
        Lzeta_6  = Lzeta_6,
        Lzeta_7  = Lzeta_7,
        Lzeta_8  = Lzeta_8,
        Lzeta_9  = Lzeta_9,
        Lzeta_10 = Lzeta_10,
        Lzeta_11 = Lzeta_11,
        Lzeta_12 = Lzeta_12,
        Lzeta_13 = Lzeta_13,
        Lzeta_14 = Lzeta_14,
        Lzeta_15 = Lzeta_15,
        Lzeta_16 = Lzeta_16,
        Lzeta_17 = Lzeta_17,
        Lzeta_18 = Lzeta_18,
        Lzeta_19 = Lzeta_19,
        Lzeta_20 = Lzeta_20,
        Lzeta_21 = Lzeta_21,
        Lzeta_22 = Lzeta_22,
        Lzeta_23 = Lzeta_23,
        Lzeta_24 = Lzeta_24,
        Lzeta_25 = Lzeta_25,
        Lzeta_26 = Lzeta_26,
        Lzeta_27 = Lzeta_27,
        Lzeta_28 = Lzeta_28,
        Lzeta_29 = Lzeta_29,
        Lzeta_30 = Lzeta_30,
        Lzeta_31 = Lzeta_31,
        Lzeta_32 = Lzeta_32,

        Laeta = Laeta, 
        Ldeta = Ldeta, 
        which_dim_cor_reg = which_dim_cor_reg,
        rand_cor_g1 = rand_cor_g1,
        rand_cor_g2 = rand_cor_g2,
        u_Lzeta_cor = u_Lzeta_cor, 
        u_Lzeta_cor_2 = u_Lzeta_cor_2, 
        u_Lzeta_cor_3 = u_Lzeta_cor_3, 
        l_Lzeta_cor = l_Lzeta_cor, 
        l_Lzeta_cor_2 = l_Lzeta_cor_2,
        l_Lzeta_cor_3 = l_Lzeta_cor_3,
        u_Laeta_cor = u_Laeta_cor, 
        l_Laeta_cor = l_Laeta_cor, 
        u_Ldeta_cor = u_Ldeta_cor, 
        l_Ldeta_cor = l_Ldeta_cor, 
        Lzeta_cor = Lzeta_cor,
        Lzeta_cor_2 = Lzeta_cor_2,
        Lzeta_cor_3 = Lzeta_cor_3,
        Laeta_cor = Laeta_cor, 
        Ldeta_cor = Ldeta_cor, 
        z = z, 

        zLong    = zLong,
        zLong_2  = zLong_2,
        zLong_3  = zLong_3,
        zLong_4  = zLong_4,
        zLong_5  = zLong_5,
        zLong_6  = zLong_6,
        zLong_7  = zLong_7,
        zLong_8  = zLong_8,
        zLong_9  = zLong_9,
        zLong_10 = zLong_10,
        zLong_11 = zLong_11,
        zLong_12 = zLong_12,
        zLong_13 = zLong_13,
        zLong_14 = zLong_14,
        zLong_15 = zLong_15,
        zLong_16 = zLong_16,
        zLong_17 = zLong_17,
        zLong_18 = zLong_18,
        zLong_19 = zLong_19,
        zLong_20 = zLong_20,
        zLong_21 = zLong_21,
        zLong_22 = zLong_22,
        zLong_23 = zLong_23,
        zLong_24 = zLong_24,
        zLong_25 = zLong_25,
        zLong_26 = zLong_26,
        zLong_27 = zLong_27,
        zLong_28 = zLong_28,
        zLong_29 = zLong_29,
        zLong_30 = zLong_30,
        zLong_31 = zLong_31,
        zLong_32 = zLong_32,

        z_2 = z_2,
        z_3 = z_3,
        ar = ar,
        dr = dr,

        Lzeta_sd    = Lzeta_sd, 
        Lzeta_sd_2  = Lzeta_sd_2, 
        Lzeta_sd_3  = Lzeta_sd_3, 
        Lzeta_sd_4  = Lzeta_sd_4,
        Lzeta_sd_5  = Lzeta_sd_5,
        Lzeta_sd_6  = Lzeta_sd_6,
        Lzeta_sd_7  = Lzeta_sd_7,
        Lzeta_sd_8  = Lzeta_sd_8,
        Lzeta_sd_9  = Lzeta_sd_9,
        Lzeta_sd_10 = Lzeta_sd_10,
        Lzeta_sd_11 = Lzeta_sd_11,
        Lzeta_sd_12 = Lzeta_sd_12,
        Lzeta_sd_13 = Lzeta_sd_13, 
        Lzeta_sd_14 = Lzeta_sd_14,
        Lzeta_sd_15 = Lzeta_sd_15,
        Lzeta_sd_16 = Lzeta_sd_16,
        Lzeta_sd_17 = Lzeta_sd_17,
        Lzeta_sd_18 = Lzeta_sd_18,
        Lzeta_sd_19 = Lzeta_sd_19,
        Lzeta_sd_20 = Lzeta_sd_20,
        Lzeta_sd_21 = Lzeta_sd_21,
        Lzeta_sd_22 = Lzeta_sd_22,
        Lzeta_sd_23 = Lzeta_sd_23, 
        Lzeta_sd_24 = Lzeta_sd_24,
        Lzeta_sd_25 = Lzeta_sd_25,
        Lzeta_sd_26 = Lzeta_sd_26,
        Lzeta_sd_27 = Lzeta_sd_27,
        Lzeta_sd_28 = Lzeta_sd_28,
        Lzeta_sd_29 = Lzeta_sd_29,
        Lzeta_sd_30 = Lzeta_sd_30,
        Lzeta_sd_31 = Lzeta_sd_31,
        Lzeta_sd_32 = Lzeta_sd_32,

        zeta_sd_ind    = zeta_sd_ind, 
        zeta_sd_ind_2  = zeta_sd_ind_2, 
        zeta_sd_ind_3  = zeta_sd_ind_3, 
        zeta_sd_ind_4  = zeta_sd_ind_4, 
        zeta_sd_ind_5  = zeta_sd_ind_5,
        zeta_sd_ind_6  = zeta_sd_ind_6,
        zeta_sd_ind_7  = zeta_sd_ind_7,
        zeta_sd_ind_8  = zeta_sd_ind_8,
        zeta_sd_ind_9  = zeta_sd_ind_9,
        zeta_sd_ind_10 = zeta_sd_ind_10,
        zeta_sd_ind_11 = zeta_sd_ind_11,
        zeta_sd_ind_12 = zeta_sd_ind_12,
        zeta_sd_ind_13 = zeta_sd_ind_13,
        zeta_sd_ind_14 = zeta_sd_ind_14, 
        zeta_sd_ind_15 = zeta_sd_ind_15,
        zeta_sd_ind_16 = zeta_sd_ind_16,
        zeta_sd_ind_17 = zeta_sd_ind_17,
        zeta_sd_ind_18 = zeta_sd_ind_18,
        zeta_sd_ind_19 = zeta_sd_ind_19,
        zeta_sd_ind_20 = zeta_sd_ind_20,
        zeta_sd_ind_21 = zeta_sd_ind_21,
        zeta_sd_ind_22 = zeta_sd_ind_22,
        zeta_sd_ind_23 = zeta_sd_ind_23,
        zeta_sd_ind_24 = zeta_sd_ind_24, 
        zeta_sd_ind_25 = zeta_sd_ind_25,
        zeta_sd_ind_26 = zeta_sd_ind_26,
        zeta_sd_ind_27 = zeta_sd_ind_27,
        zeta_sd_ind_28 = zeta_sd_ind_28,
        zeta_sd_ind_29 = zeta_sd_ind_29,
        zeta_sd_ind_30 = zeta_sd_ind_30,
        zeta_sd_ind_31 = zeta_sd_ind_31,
        zeta_sd_ind_32 = zeta_sd_ind_32,

        zeta_sd_ind_diag    = zeta_sd_ind_diag,
        zeta_sd_ind_diag_2  = zeta_sd_ind_diag_2,
        zeta_sd_ind_diag_3  = zeta_sd_ind_diag_3,
        zeta_sd_ind_diag_4  = zeta_sd_ind_diag_4,
        zeta_sd_ind_diag_5  = zeta_sd_ind_diag_5,
        zeta_sd_ind_diag_6  = zeta_sd_ind_diag_6,
        zeta_sd_ind_diag_7  = zeta_sd_ind_diag_7,
        zeta_sd_ind_diag_8  = zeta_sd_ind_diag_8,
        zeta_sd_ind_diag_9  = zeta_sd_ind_diag_9,
        zeta_sd_ind_diag_10 = zeta_sd_ind_diag_10,
        zeta_sd_ind_diag_11 = zeta_sd_ind_diag_11,
        zeta_sd_ind_diag_12 = zeta_sd_ind_diag_12,
        zeta_sd_ind_diag_13 = zeta_sd_ind_diag_13,
        zeta_sd_ind_diag_14 = zeta_sd_ind_diag_14,
        zeta_sd_ind_diag_15 = zeta_sd_ind_diag_15,
        zeta_sd_ind_diag_16 = zeta_sd_ind_diag_16,
        zeta_sd_ind_diag_17 = zeta_sd_ind_diag_17,
        zeta_sd_ind_diag_18 = zeta_sd_ind_diag_18,
        zeta_sd_ind_diag_19 = zeta_sd_ind_diag_19,
        zeta_sd_ind_diag_20 = zeta_sd_ind_diag_20,
        zeta_sd_ind_diag_21 = zeta_sd_ind_diag_21,
        zeta_sd_ind_diag_22 = zeta_sd_ind_diag_22,
        zeta_sd_ind_diag_23 = zeta_sd_ind_diag_23,
        zeta_sd_ind_diag_24 = zeta_sd_ind_diag_24,
        zeta_sd_ind_diag_25 = zeta_sd_ind_diag_25,
        zeta_sd_ind_diag_26 = zeta_sd_ind_diag_26,
        zeta_sd_ind_diag_27 = zeta_sd_ind_diag_27,
        zeta_sd_ind_diag_28 = zeta_sd_ind_diag_28,
        zeta_sd_ind_diag_29 = zeta_sd_ind_diag_29,
        zeta_sd_ind_diag_30 = zeta_sd_ind_diag_30,
        zeta_sd_ind_diag_31 = zeta_sd_ind_diag_31,
        zeta_sd_ind_diag_32 = zeta_sd_ind_diag_32,

        cor_z_item_ind = cor_z_item_ind, 
        cor_z_item_elem_ind = cor_z_item_elem_ind, 
        cor_z_item_ind_2 = cor_z_item_ind_2, 
        cor_z_item_elem_ind_2 = cor_z_item_elem_ind_2, 
        cor_z_item_ind_3 = cor_z_item_ind_3, 
        cor_z_item_elem_ind_3 = cor_z_item_elem_ind_3, 
        Laeta_sd = Laeta_sd, 
        alindex = alindex, 
        aeta_sd_ind = aeta_sd_ind, 
        cor_a_item_ind = cor_a_item_ind, 
        cor_a_item_elem_ind = cor_a_item_elem_ind,
        Ldeta_sd = Ldeta_sd, 
        dlindex = dlindex, 
        deta_sd_ind = deta_sd_ind, 
        cor_d_item_ind = cor_d_item_ind, 
        cor_d_item_elem_ind = cor_d_item_elem_ind,
        z_c = z_c,
        z_c_2 = z_c_2,
        z_c_3 = z_c_3,
        z_cLong = z_cLong,
        z_cLong_2 = z_cLong_2,
        z_cLong_3 = z_cLong_3,
        a_c = a_c, 
        d_c = d_c
      )

  } else {    # D > 1 & h2_dims > 0

    standata = list(
        N = N, 
        J = J, 
        K = K, 
        any_rand = any_rand, 
        any_rand_ind = any_rand_ind, 
        any_rand_cor = any_rand_cor, 
        any_rand_ind_a = any_rand_ind_a, 
        any_rand_cor_a = any_rand_cor_a, 
        any_rand_ind_d = any_rand_ind_d, 
        any_rand_cor_d = any_rand_cor_d,
        Ncateg_max = Ncateg_max, 
        Ncategi = Ncategi, 
        N_long = N_long, 
        nn = nn, 
        jj = jj, 
        y = y, 
        itype = itype, 
        any_eta3pl = any_eta3pl, 
        nEta3pl = nEta3pl, 
        find_eta3pl = find_eta3pl, 
        x = x, 
        D = D, 
        DAlpha = DAlpha, 
        alpha_dstart = alpha_dstart,
        alpha_dend = alpha_dend,
        nDelta = nDelta, 
        lambda_ind = lambda_ind,
        L = L, 
        has_treg = has_treg, 
        beta_dstart = beta_dstart, 
        beta_dend = beta_dend, 
        zeta_dstart = zeta_dstart, 
        zeta_dend = zeta_dend, 
        fweights = fweights, 
        # x_miss = x_miss, 
        nDelta_r = nDelta_r, 
        nAlpha_r = nAlpha_r, 
        LMean = LMean, 
        deltaMean = deltaMean, 
        d_design = d_design, 
        a_design = a_design, 
        Lzeta = Lzeta, 
        Laeta = Laeta, 
        Ldeta = Ldeta, 
        u_Lzeta_cor = u_Lzeta_cor, 
        l_Lzeta_cor = l_Lzeta_cor, 
        u_Laeta_cor = u_Laeta_cor, 
        l_Laeta_cor = l_Laeta_cor, 
        u_Ldeta_cor = u_Ldeta_cor, 
        l_Ldeta_cor = l_Ldeta_cor, 
        Lzeta_cor = Lzeta_cor, 
        Laeta_cor = Laeta_cor, 
        Ldeta_cor = Ldeta_cor, 
        z = z, 
        ar = ar,
        dr = dr,
        Lzeta_sd = Lzeta_sd, 
        zeta_sd_ind = zeta_sd_ind, 
        cor_z_item_ind = cor_z_item_ind, 
        cor_z_item_elem_ind = cor_z_item_elem_ind, 
        Laeta_sd = Laeta_sd, 
        alindex = alindex, 
        aeta_sd_ind = aeta_sd_ind, 
        cor_a_item_ind = cor_a_item_ind, 
        cor_a_item_elem_ind = cor_a_item_elem_ind,
        Ldeta_sd = Ldeta_sd, 
        dlindex = dlindex, 
        deta_sd_ind = deta_sd_ind, 
        cor_d_item_ind = cor_d_item_ind, 
        cor_d_item_elem_ind = cor_d_item_elem_ind,
        z_c = z_c,
        a_c = a_c, 
        d_c = d_c
      )

  }

  # Select the correct inirt implementation based on input parameters
  if(D == 1) {

    modl = stanmodels$instrument_unidim
    mtype = "unidim"

  } else if(D > 1 & h2_dims == 0) {

    modl = stanmodels$instrument_mirt
    mtype = "mdim"

  } else {    # D > 1 & h2_dims > 0

    modl = stanmodels$instrument_soirt
    mtype = "soirt"

  }

  # Choose a model estimateion method: variational inference or HMC
  if(method[1] == "vb") {

    out = rstan::vb(modl, data = standata, ...)

  } else if(method[1] == "hmc") {

    out = rstan::sampling(modl, data = standata, ...)

  } else {

    stop("Method not yet implemeted.")

  }

  output = list(
    stanfit = out,
    args = list(match.call()),
    mtype = mtype
  )

  class(output) = 'instrumentObj'

  return(output)

}
