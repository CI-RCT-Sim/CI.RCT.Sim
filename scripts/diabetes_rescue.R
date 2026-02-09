devtools::load_all()

# Define the parameter valued, calculate derived quantities from the parameters
sim_parameters <- assumptions_diabetes_rescue()[1, ] |>
  true_summary_statistics_diabetes_rescue()

# constants for simulation and aggregation
N_sim <- 10 # 00
alpha <- 0.1

# list of analysis functions
my_analyse <- list(
  # ipw = analyse_ipw(estimand = "tp") # ,
  dm = analyse_diabetes_demediation()
)

# list of summarisation functions
# those are applied to the results of the analysis functions by the same name,
# repeating one name multiple times is allowed
# summarise_estimator and summarise_test are generic summarisation function
# from SimNPH, but own functions can be defined.
my_summarise <- create_summarise_function(
  # ipw = summarise_estimator(est = coef, real = eff_true) # ,
  dm = summarise_estimator(est = coef, real = eff_true)
  # ttest = summarise_test(alpha),
  # lm = function(condition, results, fixed_objects = NULL) {
  #   data.frame(mean_ci_width = mean(results$ci_upper - results$ci_lower))
  # }
)

# run the simulations using SimDesign::runSimulation
results <- runSimulation(
  sim_parameters,
  replications = N_sim,
  generate = generate_diabetes_rescue,
  analyse = my_analyse,
  summarise = my_summarise
)

# inspect the results
results |>
  subset(select = c(names(sim_parameters), "dm.bias", "dm.sd_est"))
