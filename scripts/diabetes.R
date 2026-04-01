# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)


# Derive true treatment effect  under treatment policy -------------------

pre_sim_parameters <- diabetes_scenario() |>
  dplyr::mutate(nfix = 1e6, miss = rep(list(c(-1e5, 0, 0, 0)), 16))

pre_N_sim <- 20

pre_my_analyse <- list(
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

pre_my_summarise <- create_summarise_function(
  tp_mean = summarise_estimator(
    est = est,
    real = 0,
    null = 0
  )
)

cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

pre_results <- runSimulation(
  design = pre_sim_parameters,
  replications = pre_N_sim,
  generate = generate_diabetes,
  analyse = pre_my_analyse,
  summarise = pre_my_summarise,
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Define parameter values and derived quantities -------------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues() |>
  dplyr::mutate(tp_eff = pre_results$tp_mean.mean_est)

# Constants for simulation -----------------------------------------------

N_sim <- 1000

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
    real = tp_eff,
    lower = ci_lower,
    upper = ci_upper,
    null = 0,
    name="est"
  ),
  mmrmtp = summarise_estimator(
    est   = coef,
    real  = tp_eff,
    lower = ci_lower,
    upper = ci_upper,
    null  = 0,
    name="est"
  ),
  mitp = summarise_estimator(
    est = coef,
    real = tp_eff,
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

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

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

save(results, main_sessioninfo, nodes_sessioninfo, file=format(Sys.time(), paste0("results_diabetes_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")))
