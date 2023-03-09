# INIRT higher-order IRT tests
library(devtools)
library(Rcpp)
compileAttributes()
document()
devtools::install(dependencies = FALSE)
rstan::stanc(file = "./inst/stan/theta2_mirt.stan", verbose = TRUE)

library(theta2)
library(foreach)

# source simulation function
# source("tests/sim_mirt_data.R")
#---
# itype = "2pl"
# method = "hmc"
# iter = 500
# warmup = 300
# chains = 1
# fweights = NULL
# cores = 1
# exploratory = TRUE

# make a cluster and run simulation in parallel
make_parallel_compute = function(n_sim = 100, n_cores = parallel::detectCores()) {

  # locally sourced data simulation function
  source("tests/sim_mirt_data.R", local = TRUE)

  # replace this with parallely
  # Construct cluster
  cl = parallelly::makeClusterPSOCK(n_cores, autoStop = TRUE)

  # After the function is run, shutdown the cluster.
  on.exit(parallel::stopCluster(cl))

  # Register parallel backend
  doParallel::registerDoParallel(cl)   # Modify with any do*::registerDo*()

  # Execute for loop in parallel to get simulation results
  estimates = foreach::foreach(i = iterators::icount(n_sim),
                               .packages = "theta2") %dopar% {
    # seed with formula 
    set.seed(i * 10 / pi)

    # simulate mirt data
    mirt_data = sim_mirt_data(type = "mirt")
    fit_data = mirt_data$fit_data  # data for fitting a model
    sim_data = mirt_data$sim_data  # underlying truth

    # truth
    true = sim_data$true

    # simulated data set
    data = fit_data$data

    # model
    model = "theta1 = c(1:80)
             theta2 = c(1:80)
             theta3 = c(1:80)
             theta4 = c(1:80)"
    
    # fit model
    fit = theta2::theta2(data = data, model = model, itype = "2pl", exploratory = TRUE, 
      method = "hmc", iter = 10, warmup = 5, chains = 1, cores = 1)

    # produce summary
    fit_smy = theta2::summary.theta2Obj(fit)

    # merge summary output with truth
    fit_smy = fit_smy[true, on = 'parameter']

    # return summary
    fit_smy
  }

  # Release results
  return(estimates)
}

# run simulation
res = make_parallel_compute(n_sim = 2, n_cores = 6)

# assess model with metrics: mse, bias, etc.
mod_assessment = evaluate_model(res)









#array(c(0.0), dim = 1)
library(rstan)

sim_data$theta
sim_data$theta_g
cor(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta_g)
plot(rstan::summary(fit, pars = c("theta"))$summary[,1], sim_data$theta_g)
theta_resid = matrix(rstan::summary(fit, pars = c("theta_resid"))$summary[,1], ncol = 4, byrow = TRUE)
tg = rstan::summary(fit, pars = c("theta"))$summary[,1]
lam = rstan::summary(fit, pars = c("lambda_identify"))$summary[,1]

rstan::summary(fit, pars = c("sig_sq_thetag_reg"))

plot(tg, sim_data$theta_g)

plot(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])

plot(theta_resid[,1], sim_data$theta[,1])
plot(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
plot(theta_resid[,2], sim_data$theta[,2])


plot(tg*lam[4] + theta_resid[,4], sim_data$theta_g*sim_data$lambda[4] + sim_data$theta[,4])


hist(tg*lam[1] + theta_resid[,1])

# plot(theta_resid[,1], sim_data$theta[,1])
cor(tg, sim_data$theta_g)
plot(tg, sim_data$theta_g)
plot(theta_resid[,1], sim_data$theta[,1])
cor(theta_resid[,1], sim_data$theta[,1])
cor(theta_resid[,2], sim_data$theta[,2])
cor(theta_resid[,3], sim_data$theta[,3])
cor(theta_resid[,4], sim_data$theta[,4])

df = data.frame(x1 = tg*lam[1] + theta_resid[,1], 
                x2 = sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
ggplot(df) + 
  aes(x = x1, y = x2) + 
  geom_point() + 
  xlim(-1.5, 1.5)

plot(tg*lam[1] + theta_resid[,1], sim_data$theta_g*sim_data$lambda[1] + sim_data$theta[,1])
cor(tg*lam[2] + theta_resid[,2], sim_data$theta_g*sim_data$lambda[2] + sim_data$theta[,2])
cor(tg*lam[3] + theta_resid[,3], sim_data$theta_g*sim_data$lambda[3] + sim_data$theta[,3])
cor(tg*lam[4] + theta_resid[,4], sim_data$theta_g*sim_data$lambda[4] + sim_data$theta[,4])

sim_data$lambda
rstan::summary(fit, pars = c("lambda"))$summary[,1]


rstan::summary(fit, pars = c("sig_thetag_reg"))$summary[,1]

rstan::summary(fit, pars = c("lambda_identify"))$summary[,1]

rstan::summary(fit, pars = c("beta_l"))$summary[,1]

library(rstan)
traceplot(fit, pars = c("lambda_identify"))
traceplot(fit, pars = c("sig_thetag_reg"))
traceplot(fit, pars = c("theta_resid[1,1]", "theta_resid[2,1]"))

rstan::summary(fit, pars = c("sig_thetag_reg"))

# fit = inirt::inirt(data = fit_data$data, model = fit_data$model, predictors = fit_data$predictors, dims = fit_data$dims, 
#     h2_dims = fit_data$h2_dims, h2_dim_id = fit_data$h2_dim_id, structural_design = fit_data$structural_design, 
#     method = fit_data$method, weights = fit_data$weights, iter = fit_data$iter, 
#     init = fit_data$init)


# summary(fit, pars = "alpha")$summary[,"mean"]
aest = matrix(rstan::summary(fit, pars = "alpha")$summary[,1], nrow = 4, byrow = TRUE)
aest
plot(exp(aest[aest != 0.0]), exp(sim_data$alpha[sim_data$alpha != 0.0]))
cor(aest[1,1:5], sim_data$alpha[1,1:5])
plot(exp(aest[1,1:10]), exp(sim_data$alpha[1,1:10]))
cor(aest[2,21:40], alpha[2,21:40])
cor(aest[3,41:60], alpha[3,41:60])

dest = matrix(rstan::summary(fit, pars = c("delta_trans"))$summary[,1], nrow = 20, byrow = TRUE)
cor(dest[,1], sim_data$delta[,2])
cor(dest[,2], delta[,3])
cor(dest[,3], delta[,4])

tgest = summary(fit, pars = c("theta_g"))$summary[,"mean"]
cor(tgest, theta_g)
plot(tgest, theta_g)
test = matrix(summary(fit, pars = c("theta_resid"))$summary[,"mean"], nrow = n, byrow = TRUE)
cor(tgest + test[,1], theta[,1])
cor(test[,1], theta[,1])
cor(test[,2], theta[,2])
cor(test[,3], theta[,3])

cor(tgest + test[,2], theta[,2])
cor(tgest + test[,3], theta[,3])

summary(fit, pars = c("lambda"))$summary[,"mean"]
