#' Parse inirt model
#' 
#' Parse inirt model. For internal use only.
#' 
#' @param model text string which gives model definition
#' @param data data frame passed from user
#' @param exploratory do exploratory multidimensional IRT? TRUE/FALSE
#' 
#' This function takes a user inputted model fron inirt::inirt and the data,
#' translates the input and returns the data which stan and inirt::inirt require to 
#' fit the appropriate model.
#' 
#' @return list containing translated data to set up and fit the model
#' 
parse_model = function(model, data, exploratory = FALSE){

  # Remove all whitespace from model definition
  model = str_replace_all(model, pattern=" ", repl="")

  # remove all '+\n' where model definition is broked to a new line by user
  # e.g., y ~ x1 + x2 + 
  #          x3 + x4
  model = str_replace_all(model, "\\+\n", "+")

  # Splot the sub-models into a vector; each sub-model is separated by a new line 
  # by the user; no semi-colors are permitted to break lines
  model = str_split(model, "\n")[[1]]

  # remove extra newlines if the user inserted them; usually occur at the 
  # beginning or the end of the model definition
  model = model[model != ""]

  # separate latent dimension definitions from regression equations, collect
  # latent dimension definitions into one vector
  latent_models = str_detect(model, "=")
  # separate the two problems
  model_latent = model[latent_models]
  model_regression = model[!latent_models]

  # How many sub-models were defined?
  model_length = length(model_regression) + 1

  # list to store processed data for each sub-model
  input_data = vector(mode = "list", length = model_length)

  # parse the latent definition of the model, store as first element of input_data
  input_data[[1]] = parse_theta_eq(model_latent, data = data, exploratory = exploratory)
  
  # loop through the vector of models
  for(i in 1:(model_length - 1)) {
    # if detect a regression model, else detect & process a dimension definition
    # if(str_detect(model[i], "~")) {
    input_data[[i + 1]] = parse_regression_eq(model = model_regression[i], data = data)
    # } else if(str_detect(model[i], "=")) {
      # input_data[[i]] = parse_theta_eq(model = model[i], data = data, exploratory = exploratory)
    # }
  }

  # return list of input data for each sub-model given in the input
  return(input_data)

}

# item_id = str_squish(unlist(str_split(unlist(str_split(mod_theta, "="))[2], "\\+")))
# data = as.data.frame(matrix(0, 10, 50))
# names(data) = paste0("x", 1:50)
# data$School = paste0("s", rep(1:5, each = 2))
# data$age = runif(10, 10, 20)

# model_data = parse_model(
#   "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
#    t2 = x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18
#    t1 ~ x12 + x13 + x15
#    t2 ~ 0
#    alpha ~ x22 + x20 + x19
#    delta ~ x40 + x41 + x42 + x50",
#   data = data, 
#   exploratory = FALSE)

# model_data = parse_model(
#   "theta = c(1:50)
#    theta ~ 0
#    alpha ~ 0
#    delta ~ 0",
#   data = data, 
#   exploratory = FALSE)

# str(model_data)
# model_data[[2]]
