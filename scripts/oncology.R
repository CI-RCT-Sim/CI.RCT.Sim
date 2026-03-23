devtools::load_all()

# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- oncology_scenario() |>
  oncology_scenario_set_truevalues()

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 1000
alpha <- 0.05

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

my_analyse <- list(
  rpsftm = analyse_oncology_rpsftm(),
  tse = analyse_oncology_TSE(),
  # gformula = analyse_oncology_gformula(),
  # ipw = analyse_oncology_ipw(),
  itt = analyse_oncology_itt(),
  cens = analyse_oncology_cens()
)

# -------------------------------------------------------------------
# List of summarisation functions
# -------------------------------------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  rpsftm = summarise_estimator(
    est = coef,
    real = eff_true
  ),
  tse = summarise_estimator(
    est = coef,
    real = eff_true
  ),
  # gformula = summarise_estimator(
  #   est = coef,
  #   real = eff_true
  # ),
  # ipw = summarise_estimator(
  #   est = coef,
  #   real = eff_true
  # ),
  itt = summarise_estimator(
    est = coef,
    real = eff_true
  ),
  cens = summarise_estimator(
    est = coef,
    real = eff_true
  )
)

# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_oncology,
  analyse = my_analyse,
  summarise = my_summarise
)

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

save(results, file = format(Sys.time(), "results_test_%Y-%m-%d_%H%M.Rdata"))
