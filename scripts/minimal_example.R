devtools::load_all()

# Define the parameter valued, calculate derived quantities from the parameters
sim_parameters <- assumptions_minimal_example() |>
  true_summary_statistics_minimal_example()

# constants for simulation and aggregation
N_sim <- 1000
alpha <- 0.1

# list of analysis functions
my_analyse <- list(
  lm    = analyse_minimal_example_lm(ci_level=0.9),
  ttest = analyse_minimal_example_t()
)

# list of summarisation functions
# those are applied to the results of the analysis functions by the same name,
# repeating one name multiple times is allowed
# summarise_estimator and summarise_test are generic summarisation function
# from SimNPH, but own functions can be defined.
my_summarise <- create_summarise_function(
  lm    = summarise_estimator(est=coef, real=eff_size, lower=ci_lower, upper=ci_upper, null=0),
  ttest = summarise_test(alpha),
  lm    = function(condition, results, fixed_objects = NULL){
    data.frame(mean_ci_width=mean(results$ci_upper - results$ci_lower))
  }
)

# run the simulations using SimDesign::runSimulation
results <- runSimulation(
  sim_parameters,
  replications = N_sim,
  generate = generate_minimal_example,
  analyse = my_analyse,
  summarise = my_summarise
)

# inspect the results
results |>
  subset(select=c(names(sim_parameters), "lm.bias", "lm.sd_est", "ttest.rejection_0.1", "lm.1.mean_ci_width"))
