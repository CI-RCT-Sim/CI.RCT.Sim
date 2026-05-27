# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

source("scripts/vaccine_scenario_classes.R")

# Define parameter values and derived quantities -------------------------

scenario <- Sys.getenv("scenario")

message(paste("Scenario", scenario, "selected"))

selected_scenario <- switch(
  scenario,
  A1 = vaccine_scenario_A1,
  A2 = vaccine_scenario_A2,
  B1 = vaccine_scenario_B1,
  C1 = vaccine_scenario_C1,
  D1 = vaccine_scenario_D1,
  extra = vaccine_scenario_extra,
  stop("unknown scenario selected or selection missing")
)

sim_parameters <- selected_scenario |>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize() |>
  within({
    VE = 1-rr_ps
    scenario_nr = seq_along(VE)
  })

message(paste(nrow(sim_parameters), "rows"))

# Constants for simulation -----------------------------------------------

N_sim <- 10000
alpha_ci <- 0.05
alpha_test <- c(0.05, 0.025)

# List of analysis functions ---------------------------------------------

my_analyse <- list(
  # both V and W observed
  iv       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3),
  ps_cov   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE),
  ps_nocov = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE),
  pp       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3),
  # V unobserved
  iv_vunobs       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE),
  ps_cov_vunobs   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, V_unobserved=TRUE),
  ps_nocov_vunobs = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, V_unobserved=TRUE),
  pp_vunobs       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE),
  # W unobserved
  iv_wunobs       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3, W_unobserved=TRUE),
  ps_cov_wunobs   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, W_unobserved=TRUE),
  ps_nocov_wunobs = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, W_unobserved=TRUE),
  pp_wunobs       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, W_unobserved=TRUE),
  # both V and W unobserved
  iv_vwunobs       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE, W_unobserved=TRUE),
  ps_cov_vwunobs   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, V_unobserved=TRUE, W_unobserved=TRUE),
  ps_nocov_vwunobs = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, V_unobserved=TRUE, W_unobserved=TRUE),
  pp_vwunobs       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE, W_unobserved=TRUE)
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

message(paste(length(my_analyse), "analysis functions"))

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  iv               = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_cov           = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_nocov         = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  pp               = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  iv_vunobs        = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_cov_vunobs    = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_nocov_vunobs  = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  pp_vunobs        = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  iv_wunobs        = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_cov_wunobs    = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_nocov_wunobs  = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  pp_wunobs        = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  iv_vwunobs       = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_cov_vwunobs   = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  ps_nocov_vwunobs = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  pp_vwunobs       = summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  iv               = summarise_test(alpha_test, name="test"),
  ps_cov           = summarise_test(alpha_test, name="test"),
  ps_nocov         = summarise_test(alpha_test, name="test"),
  pp               = summarise_test(alpha_test, name="test"),
  iv_vunobs        = summarise_test(alpha_test, name="test"),
  ps_cov_vunobs    = summarise_test(alpha_test, name="test"),
  ps_nocov_vunobs  = summarise_test(alpha_test, name="test"),
  pp_vunobs        = summarise_test(alpha_test, name="test"),
  iv_wunobs        = summarise_test(alpha_test, name="test"),
  ps_cov_wunobs    = summarise_test(alpha_test, name="test"),
  ps_nocov_wunobs  = summarise_test(alpha_test, name="test"),
  pp_wunobs        = summarise_test(alpha_test, name="test"),
  iv_vwunobs       = summarise_test(alpha_test, name="test"),
  ps_cov_vwunobs   = summarise_test(alpha_test, name="test"),
  ps_nocov_vwunobs = summarise_test(alpha_test, name="test"),
  pp_vwunobs       = summarise_test(alpha_test, name="test"),
  iv               = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_cov           = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_nocov         = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  pp               = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  iv_vunobs        = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_cov_vunobs    = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_nocov_vunobs  = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  pp_vunobs        = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  iv_wunobs        = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_cov_wunobs    = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_nocov_wunobs  = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  pp_wunobs        = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  iv_vwunobs       = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_cov_vwunobs   = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  ps_nocov_vwunobs = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich"),
  pp_vwunobs       = summarise_estimator(VE, VE_sandwich, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich")
)

message(paste(length(environment(my_summarise)$summarise_functions), "summarise functions"))

# Run the simulations ----------------------------------------------------

message(paste("setting up cluster with", detectCores(logical=FALSE)-1 , "cores"))

cl <- makeCluster(detectCores(logical=FALSE)-1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

clusterExport(cl = cl, varlist = c("alpha_ci", "alpha_test"))

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

message(paste("running simulations,", nrow(sim_parameters), "scenarios,", N_sim, "replications"))

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

message("stopping cluster")
stopCluster(cl)

# Save results -----------------------------------------------------------

message("saving results")
save(results, main_sessioninfo, nodes_sessioninfo, file=format(Sys.time(), paste0("results_vaccine_scenario_", scenario, "_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")))
