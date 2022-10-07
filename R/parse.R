#' Parse lmHOIRT model
#' 
#' 
parse.lmHOIRT = function(model, col_names_regr, col_names_theta) {
  
  col_names_regr = c(paste0("z", 1:10))
  col_names_theta = c(paste0("x", 1:15))
  
  model = 
  "theta_g ~ theta1 + theta2 + theta3 + theta3
   theta1 = [1,2,3,4,5]
   theta2 = [5,6,7]
   theta3 = [8,9,10,11,12]
   theta4 = [13,14,15,16,17]
  "
  
  model = 
  "theta_g ~ theta1 + theta2 + theta3 + theta3
   theta_g ~ z1 + z2 + z3 + z4
   theta1 = [x1,x2,x3,x4,x5]
   theta2 = [x5,x6,x7]
   theta3 = [x8,x9,x10,x11,x12]
   theta4 = [x13,x14,x15,x16,x17]
  "
  
  model = 
  "theta_g1 ~ theta1 + theta2 + theta3 + theta4
   theta_g2 ~ theta5 + theta6 + theta7 + theta8
   theta1 = [1,2,3,4,5]
   theta2 = [5,6,7]
   theta3 = [8,9,10,11,12]
   theta4 = [13,14,15,16,17]
   theta5 = [1,2,3,4,5]
   theta6 = [5,6,7]
   theta7 = [8,9,10,11,12]
   theta8 = [13,14,15,16,17]
  "
  
  model = 
  "theta_v1 ~ theta_g1 + theta_g2
   theta_g1 ~ theta1 + theta2 + theta3 + theta4
   theta_g2 ~ theta5 + theta6 + theta7 + theta8
   theta1 = [1,2,3,4,5]
   theta2 = [5,6,7]
   theta3 = [8,9,10,11,12]
   theta4 = [13,14,15,16,17]
   theta5 = [1,2,3,4,5]
   theta6 = [5,6,7]
   theta7 = [8,9,10,11,12]
   theta8 = [13,14,15,16,17]
  "
  
  model = 
  "theta_v1 ~ theta_g1 + theta_g2
   theta_g1 ~ theta1 + theta2 + theta3 + theta4
   theta_g2 ~ theta5 + theta6 + theta7 + theta8
   theta_v1 ~ z1 + z2 + z3 + z4
   theta1 = [1,2,3,4,5]
   theta2 = [5,6,7]
   theta3 = [8,9,10,11,12]
   theta4 = [13,14,15,16,17]
   theta5 = [18,19,20,21,22]
   theta6 = [23,24,25]
   theta7 = [26,27,28,29,30]
   theta8 = [31,32,33,34,35]
  "

  model = 
  "theta_g1 = [1,2,3,4,5,6,7,8,9,10]
   theta_g2 ~ theta1 + theta2 + theta3 + theta4
   theta_g3 = [11,12,13,14,15,16,17,18]
   theta1 = [1,2,3,4,5]
   theta2 = [5,6,7]
   theta3 = [8,9,10,11,12]
   theta4 = [13,14,15,16,17]
  "
  
  model = unlist(strsplit(model, "\n"))
  model = model[!grepl("^\\s*$", model)]
  model = trimws(model)
  
  dim_distinct = grepl("[[].*[]]", model)
  model_thetas = model[!dim_distinct] 
  model_dims = model[dim_distinct]
  
  model_thetas
  model_dims

  # ? is there a regression model?
  theta_eqs_splt = strsplit(model_thetas, "~")
  
  theta_eqs_splt[[1]][2]

  child_dims  = lapply(strsplit(theta_eqs_splt[[1]][2], "+", fixed = TRUE), \(x) {gsub("[[:space:]]", "", x)})

  any_parent_dims = function(model_dims, model_thetas) {
    a = c()
    for(i in 1:length(model_dims)) {
      a[i] = grepl(gsub("[[:space:]]", "", sub("\\=.*", "", model_dims[i])), model_thetas, fixed=TRUE)
    }
    return(model_dims[a])
  }

  dims_w_parents = any_parent_dims(model_dims, model_thetas)
dims_w_parents

  no_parents = setdiff(model_dims, dims_w_parents)
  no_parents

  col_names_regr_pattern = paste(col_names_regr, sep = " ", collapse = "|")
  
  grepl(col_names_regr_pattern, theta_eqs_splt[[2]][2])
  
  lat_reg_id = grepl(col_names_regr_pattern, lapply(1:length(theta_eqs_splt), \(x) {theta_eqs_splt[[x]][2]}))
  
  lat_regr_models = model_thetas[lat_reg_id]
  lat_regr_models
  
  model_thetas

  col_names_theta_pattern = paste(c(col_names_regr, col_names_theta), sep = " ", collapse = "|")
  
  factor_m_id = !grepl(col_names_theta_pattern, lapply(1:length(theta_eqs_splt), \(x) {theta_eqs_splt[[x]][2]}))
  
  factor_models = model_thetas[factor_m_id]
  factor_models
  
  
  
  
  
  lat_regr_models
  factor_models
  model_dims
  n_mho = length(factor_models)
  n_orders = 2
  n_first_orders = length(model_dims)

  get_diminfo = function(model_dims) {
    lapply(model_dims, FUN = \(x) {as.numeric(unlist(strsplit(sub(".*\\[([^][]+)].*", "\\1",x), ",")))})
  }
  
  dim_info = get_diminfo(model_dims)
  dim_lengths = unlist(lapply(dim_info, \(x) {length(x)}))
  
  lp_data = c(n_mho, n_orders, n_first_orders, dim_lengths)
  lp_data
  #
  return(list(lp_data = lp_data))
}
