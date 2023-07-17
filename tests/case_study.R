# ------------------------------------------------------------------------------
# libraries
library(theta2)

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
         theta1 ~ (1 + wave | id)
         theta2 ~ (1 + wave | id)
         theta3 ~ (1 + wave | id)'

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

fit = theta2::theta2(data = familyrisk, model = model, itype = '2pl', 
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
