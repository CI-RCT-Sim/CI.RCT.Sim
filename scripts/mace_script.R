devtools::load_all()

# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- mace_scenario() |> true_trt_mace()

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 1000

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

my_analyse <- list(
  cox_nocov=analyse_mace_cox_nocov,
  cox_cov=analyse_mace_cox_cov,
  ipw_nocov=analyse_mace_ipw_nocov,
  ipw_cov=analyse_mace_ipw_cov 
)


my_summarise <- create_summarise_function(
  cox_nocov=summarise_func_mace,
  cox_cov=summarise_func_mace,
  ipw_nocov=summarise_func_mace,
  ipw_cov=summarise_func_mace 
)
# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_mace,
  analyse = my_analyse,
  summarise = my_summarise
)

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------
results[,grep("power",colnames(results))]
