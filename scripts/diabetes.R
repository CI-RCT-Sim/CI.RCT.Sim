devtools::load_all()

# -------------------------------------------------------------------
# Derive true treatment effect under treatment policy
# -------------------------------------------------------------------

sim_parameters <- diabetes_scenario() |>
  dplyr::mutate(nfix = 1e6, miss = rep(list(c(-1e5, 0, 0, 0)), 16))

N_sim <- 20

my_analyse <- list(
  tp_mean = function(condition, dat, fixed_objects = NULL) {
    new_dat <- dat |>
      dplyr::group_by(trt) |>
      dplyr::summarise(
        chg = base::mean(y12 - y0),
        .groups = "drop"
      )
    if (nrow(new_dat) != 2) stop("Expected exactly 2 treatment groups")
    list(est = new_dat$chg[2] - new_dat$chg[1])
  }
)

my_summarise <- create_summarise_function(
  tp_mean = summarise_estimator(
    est = est,
    real = 0,
    null = 0
  )
)

preliminary_results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_diabetes,
  analyse = my_analyse,
  summarise = my_summarise,
  parallel = TRUE
)


# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues()

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 1000
alpha <- 0.05

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

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
    null = 0
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
    null = 0
  ),
  ## Hypothetical estimands
  ipwhyp = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0
  ),
  dm = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0
  ),
  gcom = summarise_estimator(
    est = coef,
    real = eff_true,
    lower = ci_lower,
    upper = ci_upper,
    null = 0
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
    null = 0
  )
)

# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

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

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

save(results, file = format(Sys.time(), "results_test_%Y-%m-%d_%H%M.Rdata"))
