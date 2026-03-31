# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

# Define parameter values and derived quantities -------------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues()

# Constants for simulation -----------------------------------------------

N_sim <- 1000
alpha <- 0.05

# List of analysis functions ---------------------------------------------

my_analyse <- list(
  ## Treatment policy estimands
  ipwtp = analyse_diabetes_ipw(strategy = "treatment_policy"),
  mmrmtp = analyse_diabetes_mmrm(strategy = "treatment_policy"),
  mitp = analyse_diabetes_mi(strategy = "treatment_policy"),
  ## Hypothetical estimands
  ipwhyp = analyse_diabetes_ipw(strategy = "hypothetical"),
  dm = analyse_diabetes_demediation(),
  gcom = analyse_diabetes_gcomputation(),
  mmrmhyp = analyse_diabetes_mmrm(strategy = "hypothetical"),
  mihyp = analyse_diabetes_mi(strategy = "hypothetical")
)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  ## Treatment policy estimands
  ipwtp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  mmrmtp = summarise_estimator(
    est   = coef,
    real  = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0,
    name="est"
  ),
  mitp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  ## Hypothetical estimands
  ipwhyp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  dm = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  gcom = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  mmrmhyp = summarise_estimator(
    est   = coef,
    real  = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0,
    name="est"
  ),
  mihyp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  # rejection rates
  ## Treatment policy estimands
  ipwtp = summarise_estimator(
    alpha,
    name="test"
  ),
  mmrmtp = summarise_test(
    alpha,
    name="test"
  ),
  mitp = summarise_test(
    alpha,
    name="test"
  ),
  ## Hypothetical testimands
  ipwhyp = summarise_test(
    alpha,
    name="test"
  ),
  dm = summarise_test(
    alpha,
    name="test"
  ),
  gcom = summarise_test(
    alpha,
    name="test"
  ),
  mmrmhyp = summarise_test(
    alpha,
    name="test"
  ),
  mihyp = summarise_test(
    alpha,
    name="test"
  ),
)

# Run the simulations ----------------------------------------------------

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

clusterExport(cl = cl, varlist = c("alpha"))

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_diabetes,
  analyse = my_analyse,
  summarise = my_summarise,
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Save results -----------------------------------------------------------

save(results, file = format(Sys.time(), "results_diabetes_%Y-%m-%d_%H%M.Rdata"))
