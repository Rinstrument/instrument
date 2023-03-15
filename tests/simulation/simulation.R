#-------------------------------------------------------------------------------
# Author: Michael J Kleinsasser
# Title: Simulation Study of theta2 for the manuscript ___
# Date: March 9, 2023
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# import packages
library(theta2)
library(foreach)
library(tidyverse)

#-------------------------------------------------------------------------------
# Globally sources files
source("tests/simulation/evaluate_model.R")
source("tests/simulation/sim_mirt_data.R")
source("tests/simulation/tables.R")

#-------------------------------------------------------------------------------
# study level seed
seed_study_level = 422578599
set.seed(seed_study_level)

#-------------------------------------------------------------------------------
# simulate one set of truth values
true_values = sim_mirt_pars("mirt")

#-------------------------------------------------------------------------------
# make a cluster and run simulation in parallel
make_parallel_compute = function(n_sim = 100, n_cores = parallel::detectCores()) {

  # locally sourced data simulation function
  source("tests/simulation/sim_mirt_data.R", local = TRUE)

  # study level seed
  
  # seed_study_level = tail(.Random.seed, 1)

  # replace this with parallely
  # Construct cluster
  cl = parallelly::makeClusterPSOCK(n_cores, autoStop = TRUE)

  # After the function is run, shutdown the cluster.
  on.exit(parallel::stopCluster(cl))

  # Register parallel backend
  doParallel::registerDoParallel(cl)   # Modify with any do*::registerDo*()

  # Execute for loop in parallel to get simulation results
  estimates = foreach::foreach(i = iterators::icount(n_sim),
                               .packages = "theta2",
                               .export = "true_values") %dopar% {

    # seed with formula (for random sampling using true values sampled earlier)
    seed_replication_level = i * 10000 / pi
    set.seed(seed_replication_level)

    # simulate mirt data
    mirt_data = sim_mirt_data(type = "mirt", pars = true_values)
    fit_data = mirt_data$fit_data  # data for fitting a model
    sim_data = mirt_data$sim_data  # underlying truth

    # truth
    true = sim_data$true

    # simulated data set
    data = fit_data$data

    # model: first test this
    model = "theta1 = c(1:80)
             theta2 = c(1:80)
             theta3 = c(1:80)
             theta4 = c(1:80)
             theta1 ~ z1 + z2 + z3
             theta2 ~ z1 + z2 + z3
             theta3 ~ z1 + z2 + z3
             theta4 ~ z1 + z2 + z3
             alpha  ~ z1 + z2 + z3
             delta  ~ z1 + z2 + z3"

    library(devtools)
    load_all()
    itype = "2pl"
    method = "hmc"
    iter = 500
    warmup = 300
    chains = 1
    fweights = NULL
    cores = 1
    exploratory = TRUE

    # model: base simulation on this
    # model = "theta1 = c(1:80)
    #          theta2 = c(1:80)
    #          theta3 = c(1:80)
    #          theta4 = c(1:80)
    #          theta1 ~ z1 + z2 + z3 + (1 + age | school)
    #          theta2 ~ z1 + z2 + z3 + (1 + age | school)
    #          theta3 ~ z1 + z2 + z3 + (1 + age | school)
    #          theta4 ~ z1 + z2 + z3 + (1 + age | school)
    #          alpha  ~ z1 + z2 + z3
    #          delta  ~ z1 + z2 + z3"
    
    # fit model
    fit = theta2::theta2(data = data, model = model, itype = "2pl", exploratory = TRUE, 
      method = "hmc", iter = 500, warmup = 300, chains = 1, cores = 1)

    # produce summary
    fit_smy = theta2::summary.theta2Obj(fit)

    # merge summary output with truth
    fit_smy = fit_smy[true, on = 'parameter']

    # return summary from foreach iteration i
    fit_smy

  }

  # Return results
  return(estimates)

}

#-------------------------------------------------------------------------------
# run simulation
res = make_parallel_compute(n_sim = 5, n_cores = 5)
system.time(make_parallel_compute(n_sim = 5, n_cores = 5))

#-------------------------------------------------------------------------------
# assess model with metrics: mse, bias, etc.
mod_assessment = evaluate_model(res)

#-------------------------------------------------------------------------------
# Tables
tables = make_gt_tables(mod_assessment)

# Then end
#-------------------------------------------------------------------------------
