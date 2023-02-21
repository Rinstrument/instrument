
insert_default_regression = function(model_regression) {
  
  # start with one, but eventually this will need to be a loop

  # # split left hand side and right had side of model
  # model_regression = unlist(str_split(model_regression, "~"))
  # type = model_regression[1]

  split_model = str_split(model_regression, "=")
  # dim_names = purrr::map_chr(split_model, function(x) {x[1]})

  spec_type = c()
  for(i in 1:length(model_regression)) {
    type = unlist(str_split(model_regression[i], "~"))[1]
    spec_type = c(spec_type, type)
  }

  for(i in 1:length(spec_type)) {
    if(!("theta" %in% spec_type)) {
      model_regression = c(model_regression, "theta~0")
    }
    if(!("alpha" %in% spec_type)) {
      model_regression = c(model_regression, "alpha~0")
    }
    if(!("delta" %in% spec_type)) {
      model_regression = c(model_regression, "delta~0")
    }
  }

  return(model_regression)

}
