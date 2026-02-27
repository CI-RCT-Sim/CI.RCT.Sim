devtools::load_all()

# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- assumptions_diabetes_rescue()[1, ] |>
  true_summary_statistics_diabetes_rescue()

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 10 # 00
alpha <- 0.1

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

my_analyse <- list(
  ipwtp = analyse_ipw(estimand = "tp"),
  ipwhyp = analyse_ipw(estimand = "hyp"),
  dm = analyse_diabetes_demediation(),
  mmrm = analyse_diabetes_rescue_mmrm(ci_level = 1 - alpha),
  gcom = analyse_diabetes_gcomputation()
)

# -------------------------------------------------------------------
# List of summarisation functions
# -------------------------------------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  ipwtp = summarise_estimator(est = coef, real = eff_true),
  ipwhyp = summarise_estimator(est = coef, real = eff_true),
  dm = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0
  ),
  # additional custom summary: mean CI width
  dm = function(condition, results, fixed_objects = NULL) {
    data.frame(
      mean_ci_width = mean(results$ci[2] - results$ci[1], na.rm = TRUE)
    )
  },
  mmrm = summarise_estimator(
    est   = coef,
    real  = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0
  ),
  # additional custom summary: mean CI width
  mmrm = function(condition, results, fixed_objects = NULL) {
    data.frame(
      mean_ci_width = mean(results$ci_upper - results$ci_lower, na.rm = TRUE)
    )
  },
  gcom = summarise_estimator(est = coef, real = eff_true)
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

results |>
  subset(select = c(
    # names(sim_parameters),
    "ipwtp.bias",
    "ipwtp.sd_est",
    "ipwhyp.bias",
    "ipwhyp.sd_est",
    "dm.bias",
    "dm.sd_est",
    "mmrm.bias",
    "mmrm.sd_est",
    "mmrm.coverage",
    "mmrm.1.mean_ci_width",
    "gcom.bias",
    "gcom.sd_est"
  ))

saveRDS(results, "results.rds")
save.image("final_workspace.RData")
