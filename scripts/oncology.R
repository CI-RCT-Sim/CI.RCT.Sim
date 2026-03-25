devtools::load_all()

# -------------------------------------------------------------------
# Derive true treatment effect
# -------------------------------------------------------------------

sim_parameters <- oncology_scenario(print = FALSE)[1, ]

N_sim <- 10
alpha <- 0.05

my_analyse <- list(
  truth = function(condition, dat, fixed_objects = NULL) {
    mod <- coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = dat)
    HR <- exp(coef(mod)[1])
    list(HR = HR)
  }
)

my_summarise <- create_summarise_function(
  truth = summarise_estimator(
    est = HR,
    real = exp(beta_death[[1]][6]),
    null = 1
  )
)

preliminary_results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_oncology,
  analyse = my_analyse,
  summarise = my_summarise,
  fixed_objects = list(allow_switch = FALSE, logHR_assumed = NULL, ev_soll = 10000),
  seed = 1234#,
  #parallel=TRUE
)


# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- oncology_scenario(print = FALSE) |>
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
  rpsftm_rc = analyse_oncology_rpsftm(recensor = TRUE),
  rpsftm = analyse_oncology_rpsftm(recensor = FALSE),
  tse_rc = analyse_oncology_TSE(recensor = TRUE),
  tse = analyse_oncology_TSE(recensor = FALSE),
  # gformula = analyse_oncology_gformula(),
  ipw = analyse_oncology_ipw(),
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
  rpsftm_rc = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
  ),
    rpsftm = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
  ),
  tse_rc = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
  ),
    tse = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
  ),
  # gformula = summarise_estimator(
  #   est = HR,
  #   real = preliminary_results$truth.mean_est,
  #   null = 1
  # ),
  ipw = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
  ),
  itt = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
  ),
  cens = summarise_estimator(
    est = HR,
    real = preliminary_results$truth.mean_est,
    lower = low,
    upper = up,
    null = 1
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
  summarise = my_summarise,
  fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, ev_soll = NULL)
)

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

# save(results, file = format(Sys.time(), "results_test_%Y-%m-%d_%H%M.Rdata"))
