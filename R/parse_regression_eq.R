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

  if(!str_detect(model, "~")) {
    stop("regression equations require a ~ between response and predictors. E.g. a ~ b.")
  }

  predictors = NULL
  predictors_ranef = NULL
  ranef_id = NULL
  predictors_ranef_cor = NULL
  n_pranef_cor = NULL

  model = unlist(str_split(model, "~"))
  mod_lhs = str_squish(model[1])
  mod_rhs = str_squish(model[2])
  mod_fixed = str_replace_all(mod_rhs, "\\s*\\([^\\)]+\\)", "")
  mod_fixed = str_subset(str_squish(unlist(str_split(mod_fixed, "\\+"))), ".+")

  if(length(mod_fixed[1]) > 0 & mod_fixed[1] != "0") {
    predictors = which(names(data) %in% mod_fixed)
  }

  mod_ranef = unlist(str_match_all(mod_rhs, "(?<=\\().*?(?=\\))"))
  mod_ranef = str_split(mod_ranef, "\\|")

  if(length(mod_ranef) > 0) {
    ranef_id = c()
    predictors_ranef = c()
    for(i in 1:length(mod_ranef)) {
      current = mod_ranef[[i]]
      var = current[2]
      M = table(stack(setNames(strsplit(paste0(var, data[[var]]), "/"), 1:10))[2:1])
      M = matrix(M, ncol = ncol(M), dimnames = dimnames(M))
      lhs = str_squish(str_split(current[1], "\\+")[[1]])
      if(length(lhs) > 1) {
        n_pranef_cor = length(lhs)
        for(j in 2:n_pranef_cor) {
          New_M = M
          New_M[M == 1] = data[[lhs[j]]]
          M = cbind(M, New_M)
        }
        predictors_ranef_cor = c(predictors_ranef_cor, (ncol(data) + 1):(ncol(data) + ncol(M)))
        data = cbind(data, M)
      } else {
        ranef_id = c(ranef_id, rep(i, ncol(M)))
        predictors_ranef = c(predictors_ranef, (ncol(data) + 1):(ncol(data) + ncol(M)))
        if(current[1] != "1") {
          var_val = data[[current[1]]]
          colnames(M) = paste0(colnames(M), "_", current[1])
          data = cbind(data, M * var_val)
        } else {
          data = cbind(data, M)
        }
      }
    }
  }

  return(list(data = data, predictors = predictors, predictors_ranef = predictors_ranef, 
    ranef_id = ranef_id, predictors_ranef_cor = predictors_ranef_cor, 
    n_pranef_cor = n_pranef_cor))
}

data = as.data.frame(matrix(0, 10, 20))
names(data) = paste0("x", 1:20)
data$School = paste0("s", rep(1:5, each = 2))
data$age = runif(10, 10, 20)

reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + x12 + x13 + x15", data = data)
reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School) + x12 + x13 + x15", data = data)
reg_data = parse_regression_eq(model = "t1 ~ (1 + age|School)", data = data)
reg_data = parse_regression_eq(model = "t1 ~ 0", data = data)
reg_data = parse_regression_eq(model = "t1 ~ x12 + x13 + x15", data = data)

reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + c(12:20)", data = data)
reg_data = parse_regression_eq(model = "t1 ~ c(12:18) + c(22, 24)", data = data)

# Also, need to add multiline equation functionality
reg_data = parse_regression_eq(model = "t1 ~ (1|School) + (age|School) + 
                                            x12 + x13 + x15", data = data)
reg_data = parse_regression_eq(model = "t1 ~ c(12:18) + 
                                          c(22, 24)", data = data)

model = "t1 ~ c(12:18) + 
                                          c(22, 24)"
model
model = str_replace_all(model, pattern=" ", repl="") # smarter: remove all whitespace right away, then deal with other issues
model
model = str_replace_all(model, "[\r\n]" , "") # remove +\n combination
model # left with correct model
