devtools::load_all()

# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- assumptions_diabetes_rescue() |>
  true_summary_statistics_diabetes_rescue()

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 1000
alpha <- 0.1

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

my_analyse <- list(
  ## Treatment policy estimands
  ipwtp = analyse_diabetes_ipw(estimand = "tp"),
  mmrmtp = analyse_diabetes_rescue_mmrm(strategy = "treatment_policy"),
  mitp = analyse_diabetes_rescue_mi(estimand = "treatment_policy"),
  ## Hypothetical estimands
  ipwhyp = analyse_diabetes_ipw(estimand = "hyp"),
  dm = analyse_diabetes_demediation(),
  gcom = analyse_diabetes_gcomputation(),
  mmrmhyp = analyse_diabetes_rescue_mmrm(strategy = "hypothetical"),
  mihyp = analyse_diabetes_rescue_mi(estimand = "hypothetical")
)

# -------------------------------------------------------------------
# List of summarisation functions
# -------------------------------------------------------------------
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
    null  = 0
  ),
  mmrmtp = summarise_estimator(
    est   = coef,
    real  = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  ),
  mitp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  ),
  ## Hypothetical estimands
  ipwhyp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  ),
  dm = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci[1],
    upper = ci[2],
    null = 0
  ),
  gcom = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  ),
  mmrmhyp = summarise_estimator(
    est   = coef,
    real  = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  ),
  mihyp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  )
)

# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_diabetes_rescue,
  analyse = my_analyse,
  summarise = my_summarise,
  parallel = TRUE
)

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

# results |>
#   subset(select = c(
#     # names(sim_parameters),
#     "ipwtp.bias",
#     "ipwtp.sd_est",
#     "ipwhyp.bias",
#     "ipwhyp.sd_est",
#     "dm.bias",
#     "dm.sd_est",
#     "mmrm.bias",
#     "mmrm.sd_est",
#     "mmrm.coverage",
#     "mmrm.1.mean_ci_width",
#     "gcom.bias",
#     "gcom.sd_est"
#   ))

saveRDS(results, "results.rds")
save.image("final_workspace.RData")
