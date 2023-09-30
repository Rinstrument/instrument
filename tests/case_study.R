# ------------------------------------------------------------------------------
# libraries
library(instrument)

# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(devtools)
install(dependencies = FALSE)

# ------------------------------------------------------------------------------
# data
data(familyrisk)

# ------------------------------------------------------------------------------
# fit model
# get this working on the cluster tomorrow
# model = 'theta1 = c(3:16)
#          theta1 ~ (1 | id) + wave'

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         theta1 ~ (1 | id) + wave
         theta2 ~ (1 | id) + wave
         theta3 ~ (1 | id) + wave'

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         theta1 ~ (1 | id) + wave'
fit = instrument::instrument(data = familyrisk, model = model, iter = 10)

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         theta1 ~ wave + (1 + wave | id)
         theta2 ~ wave + (1 + wave | id)
         theta3 ~ wave + (1 + wave | id)'

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         thetag = theta1 + theta2 + theta3
         thetag ~ wave + (1 + wave | id)'

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         theta1 ~ wave + (1 | wave) + (1 | id)'

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         theta3 ~ wave + (1 + wave | id)'

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)
         theta1 ~ wave + (1 | id)'

# model = 'tg = theta1 + theta2 + theta3
#          theta1 = c(3:16)
#          theta2 = c(3:16)
#          theta3 = c(3:16)
#          tg ~ (1 | id) + wave'


data = familyrisk
model = model
itype = '2pl'
exploratory = TRUE
method = 'hmc'
iter = 500
chains = 1
fweights = NULL
library(devtools)
load_all()

library(tidyverse)
standata %>% names()
str(standata, list.len = length(standata))

data(familyrisk)

fit = instrument::instrument(data = familyrisk, model = model, itype = '2pl',
  exploratory = TRUE, method = 'hmc', iter = 10, chains = 1)

fit = theta2::theta2(data = familyrisk, model = model, itype = '2pl',
  exploratory = TRUE, method = 'vb', tol_rel_obj = 0.001, iter = 10000)

save(fit, file = "tests/case_study_results/fit_exploratory_3dim_wave_random_intercept_dim1.RData")

object = fit



data = ds; model = model; itype = "2pl";
exploratory = TRUE; method = "hmc"; iter = 500; chains = 1
fweights = NULL
library(devtools)
load_all()
# ------------------------------------------------------------------------------

# summarize results

library(theta2)
library(stringr)
library(tidyverse)

rm(list = ls())

# load('./tests/case_study_results/mirtRandInterceptLong.rda')
# load('./tests/case_study_results/mirtCorRanef1plusWave_theta1_version2.rda')
# load('./tests/case_study_results/mirtCorRanef1plusWave_theta2_version2.rda')
# load('./tests/case_study_results/mirtCorRanef1plusWave_theta3_version2.rda')
load('./tests/case_study_results/mirtCorRanef1plusWave_all_3_version2.rda')

library(rstan)

matrix_of_draws = as.matrix(fit$stanfit)
dim(matrix_of_draws)

Rhat(matrix_of_draws[,1])
ess_bulk(matrix_of_draws[,1])
ess_tail(matrix_of_draws[,1])


out = summary.theta2Obj(fit)

(out %>%
  pull(parameter) %>%
  str_split(., '\\[', simplify = TRUE))[, 1] %>%
  unique()

omega = matrix(out[grep('Omega', parameter), 'mean']$mean, nrow = 2)[,1:2]
tau = diag(out[grep('tau', parameter), 'mean']$mean[1:2])
vcov = tau %*% omega %*% tau
vcov

omega = matrix(out[grep('Omega', parameter), 'mean']$mean, nrow = 2)[,3:4]
tau = diag(out[grep('tau', parameter), 'mean']$mean[3:4])
vcov = tau %*% omega %*% tau
vcov

omega = matrix(out[grep('Omega', parameter), 'mean']$mean, nrow = 2)[,5:6]
tau = diag(out[grep('tau', parameter), 'mean']$mean[5:6])
vcov = tau %*% omega %*% tau
vcov

vcov_p = function(p) {
  omega = matrix(out[grep('Omega', parameter), ..p][[1]], nrow = 2)
  tau = diag(out[grep('tau', parameter), ..p][[1]])
  vcov = tau %*% omega %*% tau
  vcov
}
vcov_p(p = 'quantile_0.025')
vcov_p(p = 'quantile_0.975')

out[
    grep('alpha', parameter),
  ][
    , mean
   ] %>%
  matrix(., byrow = FALSE, nrow = 14) %>%
  as.data.frame() %>%
  mutate_all(\(x) { ifelse(x < 0.4, '--', as.character(round(x, digits = 1)))})

out[
    grep('betat', parameter),
  ]

out[
    grep('zeta', parameter),
  ]

out[
    grep('zeta\\[', parameter),
  ][
    , sd(mean)
  ]

out[
    grep('zeta_2\\[', parameter),
  ][
    , sd(mean)
  ]

out[
    grep('zeta_3\\[', parameter),
  ][
    , sd(mean)
  ]



object = fit
pars = 'default'
probs = c(0.025, 0.50, 0.975)
out = summary.theta2Obj(fit)






out[grep('betat\\[1,1\\]', parameter), ]
out[grep('Omega', parameter), ]
omega = matrix(out[grep('Omega', parameter), 'mean']$mean, nrow = 2)
tau = diag(out[grep('tau', parameter), 'mean']$mean)
vcov = tau %*% omega %*% tau
vcov

vcov[1,2] / (sqrt(vcov[1,1])*sqrt(vcov[2,2]))

out$param

reffs = cbind(
  int = out[grep('zeta_c', parameter), ][grep(',1\\]', parameter), ]$mean,
  slope = out[grep('zeta_c', parameter), ][grep(',2\\]', parameter), ]$mean
)

hist(reffs[,2])

bdat = data.frame(x = seq(-2, 2, by = 0.1))

bdat = list()
x = seq(-2, 2, by = 0.1)
for(i in 1:nrow(reffs)) {
  bdat = rbind(bdat, data.frame(ln = paste0(i), x = x,
    y = reffs[i, 1] + reffs[i, 2]*x))
}

ggplot(bdat) + aes(x = x, y = y, group = ln) +
    geom_line()


c(1:10)[rep(c(T, F), 5)]
# ------------------------------------------------------------------------------

# work out the loading piece
library(theta2)
load('tests/case_study_results/fit_exploratory_3dim.RData')

x = fit
object = fit
pars = "default"
probs = c(0.025, 0.50, 0.975)
