#' Parse Theta Equation for inirt Model
#' 
#' Parse an inirt model theta equation. For internal use only.
#' This function is called by the inirt::inirt() function to parse
#' a user specified theta equation. See examples on this page.
#' 
#' @param model string specifying regression formula
#' @param data data frame 
#' @param exploratory do exploratory multidimensional IRT? TRUE/FALSE
#' 
#' @return a list with the named elements:
#' 
#' @examples 
#' 
#' # data = as.data.frame(matrix(0, 10, 50))
#' names(data) = paste0("x", 1:50)
#' data$School = paste0("s", rep(1:5, each = 2))
#' data$age = runif(10, 10, 20)
#' mod_theta = "c(1:10, 20:35) + x12 + x13 + x14"
#'   mod_theta = "x12 + x13 + x14"
#'   mod_theta = "c(1:10, 20:35)"
#' mod = "t1 = c(1:10, 20:35) + x12 + x13 + x14"
#' mod = "t1 = c(1:10, 20:35)"
#' mod = "t1 = c(1:10)"
#' mod = "t1 = x12 + x13 + x14"
#' mod = "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8"
#' mod = "t1 = x1 + c(20:35)"
#' mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8")
#' mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", 
#'   "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8", "t4 = t1 + t2 + t3")
#' parse_regression_eq(mod, data = data)
#' mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8")
#' parse_regression_eq(mod, data = data)
#' mod = c("t1 = c(1:50)", "t2 = c(1:50)", "t3 = c(1:50)", "t4 = c(1:50)")
#' parse_regression_eq(mod, data = data, exploratory = TRUE)
#' 
#' 
#'
#'
#'
parse_theta_eq = function(model, data, exploratory = FALSE) {

  # Check if user defined the latent dimensions using correct syntax
  if(!all(str_detect(model, "="))) {
    stop("latent factor equations require a '=' between factor name and items. E.g. t1 = c(1:25).")
  }

  # remove all spaces from model definition right away
  model = str_replace_all(model, pattern=" ", repl="")

  # if equation breaks to a new line after a '+', e.g., t1 = x1 + 
  #                                                       x2 + x3
  # remove the line break \n
  model = str_replace_all(model, "[\r\n]" , "")

  # split left hand side anr right hand side; lhs defines names; rhs defines model
  split_model = str_split(model, "=")

  # store parts of model
  dim_names = purrr::map_chr(split_model, function(x) {x[1]})
  models = purrr::map_chr(split_model, function(x) {x[2]})
  model_order = str_detect(models, paste0(dim_names, collapse = "|"))

  # instantiate return values. these values store the quantities which define
  # the model for the inirt::inirt() internals
  first_order = models[!model_order]
  item_id = NULL
  dims = length(first_order)
  second_order = NULL
  h2_dims = NULL
  h2_dim_id = NULL
  dim_id = NULL

  # sort through model types: multidimesional, higher-order etc.
  if(any(model_order)) { # second-order model
    second_order = models[model_order]
    h2_dims = (!is.null(second_order)) * 1
    for(i in 1:length(first_order)) {
      model = str_split(first_order[i], c("\\+|\\,"))[[1]]
      model = unlist(str_split(str_remove_all(model, "^c\\(|\\)"), ","))

      presence = str_detect(model, ":") | str_detect(model, "\\D", negate = TRUE)

      for(j in 1:length(model)) {
        if(presence[j]) {
          new_items = eval(parse(text = model[j]))
          item_id = c(item_id, new_items)
          h2_dim_id = c(h2_dim_id, rep(i, length(new_items)))
        } else {
          name_id = which(colnames(data) %in% model[j])
          item_id = c(item_id, name_id)
          h2_dim_id = c(h2_dim_id, rep(i, 1))
        }
      }
    }
  } else {  # first-order multidimensional model
    h2_dims = 0
    # fit an exploratory muldimensional model?
    if(exploratory == TRUE) {
      model = str_split(first_order[1], c("\\+|\\,"))[[1]]
      model = unlist(str_split(str_remove_all(model, "^c\\(|\\)"), ","))

      presence = str_detect(model, ":") | str_detect(model, "\\D", negate = TRUE)

      for(j in 1:length(model)) {
        if(presence[j]) {
          new_items = eval(parse(text = model[j]))
          item_id = c(item_id, new_items)
        } else {
          name_id = which(colnames(data) %in% model[j])
          item_id = c(item_id, name_id)
        }
      }
    } else {
      for(i in 1:length(first_order)) {
        model = str_split(first_order[i], c("\\+|\\,"))[[1]]
        model = unlist(str_split(str_remove_all(model, "^c\\(|\\)"), ","))

        presence = str_detect(model, ":") | str_detect(model, "\\D", negate = TRUE)

        for(j in 1:length(model)) {
          if(presence[j]) {
            new_items = eval(parse(text = model[j]))
            item_id = c(item_id, new_items)
            dim_id = c(dim_id, rep(i, length(new_items)))
          } else {
            name_id = which(colnames(data) %in% model[j])
            item_id = c(item_id, name_id)
            dim_id = c(dim_id, rep(i, 1))
          }
        }
      }
    }
  }

  # return quantities for inirt::inirt() to set up and fit model given the 
  # parsed model definition
  return(list(item_id = item_id, dims = dims, dim_id = dim_id, h2_dims = h2_dims, 
    h2_dim_id = h2_dim_id))

}


# data = as.data.frame(matrix(0, 10, 50))
# names(data) = paste0("x", 1:50)
# data$School = paste0("s", rep(1:5, each = 2))
# data$age = runif(10, 10, 20)

# mod_theta = "c(1:10, 20:35) + x12 + x13 + x14"
#   mod_theta = "x12 + x13 + x14"
#   mod_theta = "c(1:10, 20:35)"
# mod = "t1 = c(1:10, 20:35) + x12 + x13 + x14"
# mod = "t1 = c(1:10, 20:35)"
# mod = "t1 = c(1:10)"
# mod = "t1 = x12 + x13 + x14"
# mod = "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8"
# mod = "t1 = x1 + c(20:35)"

# mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8")
# mod = c("t1 = x1 + 
#              c(20:35)", "t2 = c(1:10, 20, 35)", "t3 = x1 + x2 + x3 + x4 + x5 + 
#                                                       x6 + x7 + x8")
# parse_theta_eq(mod, data = data)

# mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", 
#   "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8", "t4 = t1 + t2 + t3")
# parse_theta_eq(mod, data = data)

# mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8")
# parse_theta_eq(mod, data = data)

# mod = c("t1 = c(1:50)", "t2 = c(1:50)", "t3 = c(1:50)", "t4 = c(1:50)")
# parse_theta_eq(mod, data = data, exploratory = TRUE)
