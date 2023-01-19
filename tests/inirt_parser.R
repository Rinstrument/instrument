
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

mod_alpha_reg = str_detect(mod, "alpha")
mod_alpha_reg = mod[mod_alpha_reg]
mod_alpha_reg

mod_delta_reg = str_detect(mod, "delta")
mod_delta_reg = mod[mod_delta_reg]
mod_delta_reg

mod_theta = str_detect(mod, " = ")
mod_theta = mod[mod_theta]
mod_theta

mod_theta_reg = str_detect(mod, " ~ ") & (!str_detect(mod, "alpha")) & (!str_detect(mod, "delta"))
mod_theta_reg = mod[mod_theta_reg]
mod_theta_reg

item_id = str_squish(unlist(str_split(unlist(str_split(mod_theta, "="))[2], "\\+")))
data = as.data.frame(matrix(0, 10, 20))
names(data) = paste0("x", 1:20)
item_id = which(names(data) %in% item_id)
item_id

mod_theta_reg = str_squish(unlist(str_split(unlist(str_split(mod_theta_reg, "~"))[2], "\\+")))
mod_theta_reg = which(names(data) %in% mod_theta_reg)
mod_theta_reg

predictors = list(mod_theta_reg)
predictors



data, item_id, model = NULL, predictors = NULL, predictors_ranef = NULL, ranef_id = NULL, 
    predictors_ranef_corr = NULL, n_pranef_cor = NULL, dims = 1, h2_dims = 0, h2_dim_id = NULL, structural_design = NULL, 
    structural_design_ranef = list(a_predictors = NULL, a_predictors_ranef = NULL, a_ranef_id = NULL, a_predictors_ranef_corr = NULL, a_n_pranef_cor = NULL,
                                   d_predictors = NULL, d_predictors_ranef = NULL, d_ranef_id = NULL, d_predictors_ranef_corr = NULL, d_n_pranef_cor = NULL),
    method = c("vb", "hmc"), weights = NULL, vb_algorithm = "fullrank", ...