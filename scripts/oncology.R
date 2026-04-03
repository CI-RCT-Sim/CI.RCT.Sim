# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

# Derive true treatment effect -------------------------------------------

pre_sim_parameters <- oncology_scenario(print = FALSE)

pre_N_sim <- 10

pre_my_analyse <- list(
  truth = function(condition, dat, fixed_objects = NULL) {
    mod <- survival::coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = dat)
    HR <- exp(coef(mod)[1])
    list(HR = HR)
  }
)

pre_my_summarise <- create_summarise_function(
  truth = summarise_estimator(
    est = HR,
    real = exp(beta_death[[1]][6]),
    null = 1
  )
)

cl <- makeCluster(41)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

pre_results <- runSimulation(
  design = pre_sim_parameters,
  replications = pre_N_sim,
  generate = generate_oncology,
  analyse = pre_my_analyse,
  summarise = pre_my_summarise,
  fixed_objects = list(allow_switch = FALSE, logHR_assumed = NULL, ev_soll = 100, allow_random_cens = TRUE),
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Define parameter values and derived quantities -------------------------

sim_parameters <- oncology_scenario(print = FALSE) |>
  oncology_scenario_set_truevalues() |>
  dplyr::mutate(true_eff = pre_results$truth.mean_est)

# Constants for simulation -----------------------------------------------

N_sim <- 10000
alpha <- 0.05

# List of analysis functions ---------------------------------------------

my_analyse <- list(
  rpsftm_rc = analyse_oncology_rpsftm(recensor = TRUE),
  rpsftm = analyse_oncology_rpsftm(recensor = FALSE),
  tse_rc = analyse_oncology_TSE(recensor = TRUE),
  tse = analyse_oncology_TSE(recensor = FALSE),
  gformula = analyse_oncology_gformula(B=200),
  ipw = analyse_oncology_ipw(),
  itt = analyse_oncology_itt(),
  cens = analyse_oncology_cens(),
  describe = function(condition, dat, fixed_objects = NULL) {
    tabulate_helper <- function(dat, var) {
      tmp <- list(
        sum(dat[, var]),
        sum(dat[dat$trt == 0, var]),
        sum(dat[dat$trt == 1, var])
      )

      names(tmp) <- c(var, paste0(var, "_ctrl"), paste0(var, "_trt"))
      tmp
    }
    result <- list(
      n_pat = nrow(dat),
      n_ctrl = sum(dat$trt == 0),
      n_trt = sum(dat$trt == 1),
      n_switch = sum(dat$switch),
      max_followup = max(dat$event_time),
      sufficient_events = sum(dat$ev) >= condition$ev_soll,
      n_random_cens = sum(dat$random_cens)
    )
    result <- c(result, tabulate_helper(dat, "ev"))
    if (!is.null(attr(dat, "followup"))) {
      result$study_time <- attr(dat, "followup")
    } else {
      result$study_time <- NA_real_
    }
    if (hasName(dat, "ice")) {
      result <- c(result, tabulate_helper(dat, "ice"))
    }
    if (hasName(dat, "subgroup")) {
      result <- c(result, tabulate_helper(dat, "subgroup"))
    }
    result
  }
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  rpsftm_rc = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  rpsftm = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  tse_rc = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  tse = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  gformula = summarise_estimator(
    est = HR,
    real = true_eff,
    null = 1,
    name = "est"
  ),
  ipw = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  itt = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  cens = summarise_estimator(
    est = HR,
    real = true_eff,
    lower = low,
    upper = up,
    null = 1,
    name = "est"
  ),
  # rejection rates
  rpsftm_rc = summarise_test(
    alpha,
    name = "test"
  ),
  rpsftm = summarise_test(
    alpha,
    name = "test"
  ),
  tse_rc = summarise_test(
    alpha,
    name = "test"
  ),
  tse = summarise_test(
    alpha,
    name = "test"
  ),
  gformula = summarise_test(
    alpha,
    name = "test"
  ),
  ipw = summarise_test(
    alpha,
    name = "test"
  ),
  itt = summarise_test(
    alpha,
    name = "test"
  ),
  cens = summarise_test(
    alpha,
    name = "test"
  ),
  describe = summarise_describe()
)

# Run the simulations ----------------------------------------------------

cl <- makeCluster(41)
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
  generate = generate_oncology,
  analyse = my_analyse,
  summarise = my_summarise,
  fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, ev_soll = NULL, allow_random_cens = TRUE),
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Save results -----------------------------------------------------------

save(results, main_sessioninfo, nodes_sessioninfo, file=format(Sys.time(), paste0("results_onco_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")))
