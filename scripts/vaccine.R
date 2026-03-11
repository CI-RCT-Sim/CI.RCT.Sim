# devtools::install()
# renv::deactivate()
# update.packages(ask=FALSE)
library(CI.RCT.Sim)
library(parallel)

# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- vaccine_scenario() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize() |>
  within({
    VE = 1-rr_ps
  })

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 5000
alpha <- 0.05

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

my_analyse <- list(
  iv = analyse_vaccine_ivreg(ci_level = 1-alpha, VE_margin = 0.3),
  ps = analyse_vaccine_ps(ci_level = 1-alpha, VE_margin = 0.3),
  pp = analyse_vaccine_pp(ci_level = 1-alpha, VE_margin = 0.3)
)

# -------------------------------------------------------------------
# List of summarisation functions
# -------------------------------------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  iv = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0, name="iv_est"),
  ps = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0, name="ps_est"),
  pp = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0, name="pp_est"),
  iv = summarise_test(alpha, name="iv_test"),
  ps = summarise_test(alpha, name="ps_test"),
  pp = summarise_test(alpha, name="pp_test")
)

# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

clusterExport(cl = cl, varlist = c("alpha"))

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

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

save(results, file=format(Sys.Date(), "results_test_%Y-%m%-%d_%H%M.Rdata"))
