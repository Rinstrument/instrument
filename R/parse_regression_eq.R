#' Parse Regression Equation for inirt Model
#' 
#' Parse an inirt model regression equation. For internal use only.
#' This function is called by the inirt::inirt() function to parse
#' a user specified regression equation. See examples on this page.
#' 
#' @param model string specifying regression formula
#' @param data data frame 
#' 
#' @return a list with the named elements:
#' data = new data frame coded for fitting the regression model given the formula
#' predictors = numeric column number for the fixed effects
#' predictors_ranef = numeric column number for the independent random effects
#' ranef_id = numeric ID for which random effect variance each random effect belongs
#' predictors_ranef_cor = numeric column number for the correlated random effects
#' n_pranef_cor = vector length of correlated random effects
#' 
#' @examples 
#' data = as.data.frame(matrix(0, 10, 20))
#' names(data) = paste0("x", 1:20)
#' data$School = paste0("s", rep(1:5, each = 2))
#' data$age = runif(10, 10, 20)
#' 
#' reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + x12 + x13 + x15", data = data)
#' reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School) + x12 + x13 + x15", data = data)
#' reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School)", data = data)
#' reg_data = parse_regression_eq(model = "t1 ~ 0", data = data)
#' reg_data = parse_regression_eq(model = "t1 ~ x12 + x13 + x15", data = data)
#'
#'
#'
parse_regression_eq = function(model, data) {

  # Check if user defined the regression equations using correct syntax
  if(!str_detect(model, "~")) {
    stop("regression equations require a ~ between response and predictors. E.g. a ~ b.")
  }

  # remove all spaces from model definition right away
  model = str_replace_all(model, pattern=" ", repl="")

  # if equation breaks to a new line, e.g., t1 ~ x1 + 
  #                                              x2 + x3
  # remove the line break \n
  model = str_replace_all(model, "[\r\n]" , "")

  # instantiate return values. these values store the quantities which define
  # the model for the inirt::inirt() internals
  predictors = NULL
  predictors_ranef = NULL
  ranef_id = NULL
  predictors_ranef_cor = NULL
  n_pranef_cor = NULL

  # split left hand side and right had side of model
  model = unlist(str_split(model, "~"))
  type = model[1]
  mod_rhs = model[2]
  # fixed effect portion of model definition
  mod_fixed = str_replace_all(mod_rhs, "\\s*\\([^\\)]+\\)", "")
  mod_fixed = str_subset(unlist(str_split(mod_fixed, "\\+")), ".+")

  if(identical(mod_fixed, character(0)) || mod_fixed[1] == "0") {
    mod_fixed = NULL
  } else {
    predictors = which(colnames(data) %in% mod_fixed)
  }

  # if(length(mod_fixed[1]) > 0 & mod_fixed[1] != "0") {
    
  # }

  # parse and store the random effect portion of the regression equation
  mod_ranef = unlist(str_match_all(mod_rhs, "(?<=\\().*?(?=\\))"))
  mod_ranef = str_split(mod_ranef, "\\|")

  # Random effects take some sorting out. They can be either uncorrelated or 
  # correlated and there can be multiple separate definitions in one formula
  # Therefore, a for loop is needed to work through something like
  # t1 ~ (1|School) + (1 + age|Cohort), etc.
  new_reg_data = matrix(0, nrow(data), 0)
  if(length(mod_ranef) > 0) {
    ranef_id = c()
    predictors_ranef = c()
    for(i in 1:length(mod_ranef)) {
      current = mod_ranef[[i]]
      var = current[2]
      M = table(stack(setNames(strsplit(paste0(var, data[[var]]), "/"), 1:10))[2:1])
      M = matrix(M, ncol = ncol(M), dimnames = dimnames(M))
      lhs = str_split(current[1], "\\+")[[1]]
      if(length(lhs) > 1) {
        n_pranef_cor = length(lhs)
        for(j in 2:n_pranef_cor) {
          New_M = M
          New_M[M == 1] = data[[lhs[j]]]
          M = cbind(M, New_M)
        }
        #predictors_ranef_cor = c(predictors_ranef_cor, (ncol(data) + 1):(ncol(data) + ncol(M)))
        predictors_ranef_cor = c(predictors_ranef_cor, 1:ncol(M))
        new_reg_data = cbind(new_reg_data, M)
      } else {
        ranef_id = c(ranef_id, rep(i, ncol(M)))
        #predictors_ranef = c(predictors_ranef, (ncol(data) + 1):(ncol(data) + ncol(M)))
        predictors_ranef = c(predictors_ranef, 1:ncol(M))
        if(current[1] != "1") {
          var_val = data[[current[1]]]
          colnames(M) = paste0(colnames(M), "_", current[1])
          new_reg_data = cbind(new_reg_data, M * var_val)
        } else {
          new_reg_data = cbind(new_reg_data, M)
        }
      }
    }
  }

  # return quantities for inirt::inirt() to set up and fit model given the 
  # parsed model definition
  return(list(type = type, new_reg_data = new_reg_data, predictors = predictors, 
    predictors_ranef = predictors_ranef, ranef_id = ranef_id, 
    predictors_ranef_cor = predictors_ranef_cor, n_pranef_cor = n_pranef_cor))
}

# data = as.data.frame(matrix(0, 10, 20))
# names(data) = paste0("x", 1:20)
# data$School = paste0("s", rep(1:5, each = 2))
# data$age = runif(10, 10, 20)

# reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + x12 + x13 + x15", data = data)
# reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School) + x12 + x13 + x15", data = data)
# reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School)", data = data)
# reg_data = parse_regression_eq(model = "t1 ~ 0", data = data)
# reg_data = parse_regression_eq(model = "t1 ~ x12 + x13 + x15", data = data)

# reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + c(12:20)", data = data)
# reg_data = parse_regression_eq(model = "t1 ~ c(12:18) + c(22, 24)", data = data)

# # Also, need to add multiline equation functionality
# reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + 
#                                             x12 + x13 + x15", data = data)
# reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School) + x12 + x13 + x15", data = data)
# reg_data
# reg_data = parse_regression_eq(model = "t1 ~ c(12:18) + 
#                                           c(22, 24)", data = data)

