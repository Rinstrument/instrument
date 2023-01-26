
mod = "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
       t1 ~ (1 + age|School) + x12 + x13 + x15
       alpha ~ a1 + a2
       delta ~ d1 + d2 + d3 + d4"
# how about definitions that carry to the next line??? - need to parse those differently
mod = "t1 = c(1:50)
       t1 ~ (1|School) + (age|School) + x12 + x13 + x15
       alpha ~ a1 + a2
       delta ~ d1 + d2 + d3 + d4"
       #t1 = c(1:10, 20:35)
       #t1 = c(1,2,5,6,10,12,14,20:30)

mod = "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
       t1 ~ (1|School) + (age|School) + x12 + x13 + x15
       alpha ~ a1 + a2
       delta ~ d1 + d2 + d3 + d4"

mod = "t1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
       t1 ~ x12 + x13 + x15
       alpha ~ a1 + a2
       delta ~ d1 + d2 + d3 + d4"

       #   t1 ~ x12 + x14
       # t2 ~ x12 + x13 + x15
       # t3 ~ x12 + x14
# t2 = x11 + x12 + x13 + x14 + x15
#        t3 = x20 + x22 + x23 + x24 + x25
#        tg = t1 + t2 + t3
library(stringr)

mod = unlist(str_split(mod, "\\n"))
mod

model = "    

t1 = x1 + x2 + x3 + x4 + x5 + 
            x6 + x7 + x8

            t2 = c(1:5, 2,4) + 
            x6 + x7 + c(1:20)
      
       t1 ~ x12 + 
           x13 + x15

       alpha ~ x1 + 
          x2

       delta ~ x4 + 
         x6
         "
model = str_replace_all(model, pattern=" ", repl="")
model
# model = str_remove(model, "^(+\n)$")
model = str_replace_all(model, "\\+\n", "+")
model
cat(model)
model = str_split(model, "\n")[[1]]
model = model[model != ""]
model

for(i in 1:length(mod)) {
  input_data = vector(mode = "list", length = length(mod))
  if(str_detect("~", mod[i])) {
    input_data[[i]] = parse_regression_eq(model = mod[i], data = data)
  } else if(str_detect("=", mod[i])) {
    input_data[[i]] = parse_theta_eq(model = mod[i], data = data)
  }
}

mod_alpha_reg = str_detect(mod, "alpha")
mod_alpha_reg = mod[mod_alpha_reg]
mod_alpha_reg

mod_delta_reg = str_detect(mod, "delta")
mod_delta_reg = mod[mod_delta_reg]
mod_delta_reg

mod_theta = str_detect(mod, " = ")
mod_theta = mod[mod_theta]
mod_theta

mod = "t1 = c(1:10, 20:35) + x12 + x13 + x14
       t1 ~ (1|School) + (age|School) + x12 + x13 + x15
       alpha ~ a1 + a2
       delta ~ d1 + d2 + d3 + d4"
       #t1 = c(1:10, 20:35)
       #t1 = c(1,2,5,6,10,12,14,20:30)

mod = unlist(str_split(mod, "\\n"))
mod

mod_theta = str_detect(mod, " = ")
mod_theta = mod[mod_theta]
mod_theta

mod_theta = str_squish(unlist(str_split(mod_theta, "="))[2])
mod_theta

# mod_theta_diff_range = str_remove_all(mod_theta, "(?<=c\\().*?(?=\\))")
# mod_theta = paste(Reduce(setdiff, strsplit(c(mod_theta, mod_theta_diff_range), split = "")), collapse = "")
# use str_remove instead of line above
mod_theta
str_replace_all("c(1,2,3)", mod_theta_diff_range, "")
str_remove_all("c(1,2,3)", "^c\\(|\\)$")
str_remove_all(mod_theta, "^c\\(|\\)$")


# grab the [:] notation for specifying a model!
mod_theta = "c(1:10, 20:35) + x12 + x13 + x14"
str_remove_all(mod_theta, "^c\\(|\\)")



item_id = str_squish(unlist(str_split(unlist(str_split(mod_theta, "="))[2], "\\+")))
data = as.data.frame(matrix(0, 10, 20))
names(data) = paste0("x", 1:20)
data$School = paste0("s", rep(1:5, each = 2))
data$age = runif(10, 10, 20)

# new theta parser!
mod_theta = "c(1:10, 20:35) + x12 + x13 + x14"
mod_theta = "x12 + x13 + x14"
mod_theta = "c(1:10, 20:35)"
mod_theta = str_squish(str_split(mod_theta, c("\\+|\\,"))[[1]])
mod_theta = str_squish(unlist(str_split(str_remove_all(mod_theta, "^c\\(|\\)"), ",")))
mod_theta
mod_theta[1]
data[, eval(parse(text = mod_theta[2]))]



mod_theta_reg = str_detect(mod, " ~ ") & (!str_detect(mod, "alpha")) & (!str_detect(mod, "delta"))
mod_theta_reg = mod[mod_theta_reg]
mod_theta_reg


item_id = str_squish(unlist(str_split(unlist(str_split(mod_theta, "="))[2], "\\+")))
data = as.data.frame(matrix(0, 10, 20))
names(data) = paste0("x", 1:20)
data$School = paste0("s", rep(1:5, each = 2))
data$age = runif(10, 10, 20)


mod_theta_reg_fixed = unlist(str_split(mod_theta_reg, "~"))[2]
mod_theta_reg_fixed = str_squish(unlist(str_split(mod_theta_reg_fixed, "\\+")))
mod_theta_reg_fixed = str_remove_all(mod_theta_reg_fixed, "(?<=\\().*?(?=\\))")
mod_theta_reg_fixed = mod_theta_reg_fixed[!str_detect(mod_theta_reg_fixed, "\\(")]
mod_theta_reg_fixed

predictors = which(names(data) %in% mod_theta_reg_fixed)
predictors


mod_theta_reg_ranef = unlist(str_match_all(mod_theta_reg, "(?<=\\().*?(?=\\))"))
mod_theta_reg_ranef = str_split(mod_theta_reg_ranef, "\\|")
mod_theta_reg_ranef

if(length(mod_theta_reg_ranef) > 0) {
  ranef_id = c()
  predictors_ranef = c()
  for(i in 1:length(mod_theta_reg_ranef)) {
    current = mod_theta_reg_ranef[[i]]
    var = current[2]
    M = table(stack(setNames(strsplit(paste0(var, data[[var]]), "/"), 1:10))[2:1])
    M = matrix(M, ncol = ncol(M), dimnames = dimnames(M))
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

item_id = which(names(data) %in% item_id)
item_id

# mod_theta_reg = str_squish(unlist(str_split(unlist(str_split(mod_theta_reg, "~"))[2], "\\+")))
# mod_theta_reg = which(names(data) %in% mod_theta_reg)
# mod_theta_reg


# predictors_ranef = which(str_detect(colnames(data), "School"))
predictors_ranef
ranef_id

# Get the parenthesis and what is inside
k = str_extract_all(" (1 + age|School) + (1| City) + x12 + x13 + x15", "\\([^()]+\\)")[[1]]
k
# Remove parenthesis
k = substring(k, 2, nchar(k) - 1)
k

str_replace_all(" (1 + age|School) + (1| City) + x12 + x13 + x15", " \\s*\\([^\\)]+\\)", "")

model = "t1 ~ (1 + age|School) + (1| City) + x12 + x13 + x15"

parse_regression_eq = function(model, data) {

  predictors = NULL
  predictors_ranef = NULL
  ranef_id = NULL
  predictors_ranef_cor = NULL
  n_pranef_cor = NULL

  if(!str_detect(model, "~")) {
    stop("regression equations require a ~ between response and predictors. E.g. a ~ b.")
  }

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

data, item_id, model = NULL, predictors = NULL, predictors_ranef = NULL, ranef_id = NULL, 
    predictors_ranef_corr = NULL, n_pranef_cor = NULL, dims = 1, h2_dims = 0, h2_dim_id = NULL, structural_design = NULL, 
    structural_design_ranef = list(a_predictors = NULL, a_predictors_ranef = NULL, a_ranef_id = NULL, a_predictors_ranef_corr = NULL, a_n_pranef_cor = NULL,
                                   d_predictors = NULL, d_predictors_ranef = NULL, d_ranef_id = NULL, d_predictors_ranef_corr = NULL, d_n_pranef_cor = NULL),
    method = c("vb", "hmc"), weights = NULL, vb_algorithm = "fullrank", ...
