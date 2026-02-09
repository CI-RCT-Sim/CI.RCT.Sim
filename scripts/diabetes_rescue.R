library("CI.RCT.Sim")
library("SimDesign")

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
  mmrm = analyse_diabetes_rescue_mmrm(ci_level = 1 - alpha)
)

# -------------------------------------------------------------------
# List of summarisation functions
# -------------------------------------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(

  # bias, SD, coverage etc. for the treatment effect at final visit
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
  }
)

# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

results <- runSimulation(
  design      = sim_parameters,
  replications = N_sim,
  generate    = generate_diabetes_rescue,
  analyse     = my_analyse,
  summarise   = my_summarise
)

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

results |>
  subset(
    select = c(
      names(sim_parameters),
      "mmrm.bias",
      "mmrm.sd_est",
      "mmrm.coverage",
      "mmrm.1.mean_ci_width"
    )
  )
