

parse_regression_eq = function(model, data) {

  model  = mod
  if(!all(str_detect(model, "="))) {
    stop("latent factor equations require a = between factor name and items. E.g. t1 = c(1:25).")
  }

  split_model = lapply(str_split(model, "="), function(x) {str_squish(x)})

  dim_names = purrr::map_chr(split_model, function(x) {x[1]})
  models = purrr::map_chr(split_model, function(x) {x[2]})
  model_order = str_detect(models, paste0(dim_names, collapse = "|"))

  first_order = models[!model_order]
  item_id = NULL
  dims = length(first_order)
  second_order = NULL
  h2_dim_id = NULL
  if(any(model_order)) {
    second_order = models[model_order]
    h2_dims = (!is.null(second_order)) * 1
    for(i in 1:length(first_order)) {
      # model = str_squish(unlist(str_split(model, "="))[2])
      model = str_squish(str_split(first_order[i], c("\\+|\\,"))[[1]])
      model = str_squish(unlist(str_split(str_remove_all(model, "^c\\(|\\)"), ",")))

      presence = str_detect(model, ":")

      for(j in 1:length(model)) {
        if(presence[j]) {
          new_items = eval(parse(text = model[j]))
          item_id = c(item_id, new_items)
          h2_dim_id = c(h2_dim_id, rep(i, length(new_items)))
        } else {
          name_id = which(names(data) %in% model[j])
          item_id = c(item_id, name_id)
          h2_dim_id = c(h2_dim_id, rep(i, 1))
        }
      }
    }
  } else {

    model = str_squish(unlist(str_split(model, "="))[2])
    model = str_squish(str_split(model, c("\\+|\\,"))[[1]])
    model = str_squish(unlist(str_split(str_remove_all(model, "^c\\(|\\)"), ",")))

    presence = str_detect(model, ":")

    for(i in 1:length(model)) {
      if(presence[i]) {
        item_id = c(item_id, eval(parse(text = model[i])))
      } else {
        name_id = which(names(data) %in% model[i])
        item_id = c(item_id, name_id)
      }
    }

  }

  return(list(item_id = item_id, dims = dims, h2_dims = h2_dims, 
    h2_dim_id = h2_dim_id))

}


data = as.data.frame(matrix(0, 10, 50))
names(data) = paste0("x", 1:50)
data$School = paste0("s", rep(1:5, each = 2))
data$age = runif(10, 10, 20)

mod_theta = "c(1:10, 20:35) + x12 + x13 + x14"
  mod_theta = "x12 + x13 + x14"
  mod_theta = "c(1:10, 20:35)"
mod = "t1 = c(1:10, 20:35) + x12 + x13 + x14"
mod = "t1 = c(1:10, 20:35)"
mod = "t1 = c(1:10)"
mod = "t1 = x12 + x13 + x14"
mod = "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8"
mod = "t1 = x1 + c(20:35)"

mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8")

mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", 
  "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8", "t4 = t1 + t2 + t3")
parse_regression_eq(mod, data = data)

mod = c("t1 = x1 + c(20:35)", "t2 = c(1:10, 20:35)", "t3 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8")
parse_regression_eq(mod, data = data)
