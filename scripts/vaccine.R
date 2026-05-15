# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

# Define parameter values and derived quantities -------------------------

sim_parameters <- vaccine_scenario_defaults() |>
  within({
    rm(beta_A2)
  }) |>
  do.call(params_scenarios_grid, args=_) |>
  merge(
    data.frame(beta_A2 = log(1 - c(0.7, 0, 0.3, 0.5, 0.9))),
    by=NULL
  ) |>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize() |>
  within({
    VE = 1-rr_ps
    scenario_nr = seq_along(VE)
  })

# Constants for simulation -----------------------------------------------

N_sim <- 10000
alpha <- 0.05

# List of analysis functions ---------------------------------------------

my_analyse <- list(
  iv       = analyse_vaccine_ivreg(ci_level = 1-alpha, VE_margin = 0.3),
  ps_cov   = analyse_vaccine_ps(ci_level = 1-alpha, VE_margin = 0.3, covariates_in_outcomes_model = TRUE),
  ps_nocov = analyse_vaccine_ps(ci_level = 1-alpha, VE_margin = 0.3, covariates_in_outcomes_model = FALSE),
  pp       = analyse_vaccine_pp(ci_level = 1-alpha, VE_margin = 0.3)
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  iv       = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_cov   = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_nocov = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  pp       = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  iv       = summarise_test(alpha, name="test"),
  ps_cov   = summarise_test(alpha, name="test"),
  ps_nocov = summarise_test(alpha, name="test"),
  pp       = summarise_test(alpha, name="test")
)

# Run the simulations ----------------------------------------------------

cl <- makeCluster(detectCores(logical=FALSE)-1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

clusterExport(cl = cl, varlist = c("alpha"))

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_vaccine,
  analyse = my_analyse,
  summarise = my_summarise,
  fixed_objects = list(include_unobserved=FALSE),
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Save results -----------------------------------------------------------

save(results, main_sessioninfo, nodes_sessioninfo, file=format(Sys.time(), paste0("results_vaccine_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")))
