# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

# Define parameter values and derived quantities -------------------------

sim_parameters <- mace_scenario() |>
  mace_scenario_set_sample_size() |>
  true_trt_mace()

# Constants for simulation -----------------------------------------------

N_sim <- 10000
alpha <- 0.05

# List of analysis functions ---------------------------------------------

my_analyse <- list(
  cox_nocov=analyse_mace_cox_nocov,
  cox_cov=analyse_mace_cox_cov,
  ipw_nocov=analyse_mace_ipw_nocov,
  ipw_cov=analyse_mace_ipw_cov
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

my_summarise <- create_summarise_function(
  cox_nocov=summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  cox_cov  =summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  ipw_nocov=summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  ipw_cov  =summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  cox_nocov=summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  cox_cov  =summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  ipw_nocov=summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  ipw_cov  =summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  cox_nocov=summarise_test(alpha=0.05, name="test"),
  cox_cov  =summarise_test(alpha=0.05, name="test"),
  ipw_nocov=summarise_test(alpha=0.05, name="test"),
  ipw_cov  =summarise_test(alpha=0.05, name="test")
)

# Run the simulations ----------------------------------------------------

cl <- makeCluster(41)
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
  generate = generate_mace,
  analyse = my_analyse,
  summarise = my_summarise,
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Save results -----------------------------------------------------------

save(results, main_sessioninfo, nodes_sessioninfo, file=format(Sys.time(), paste0("results_mace_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")))
