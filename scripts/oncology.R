# devtools::install()
# renv::restore()
devtools::load_all()
rm(list=ls())
library(CI.RCT.Sim)
library(parallel)
library(survival)

# Derive true treatment effect -------------------------------------------

sim_parameters <- oncology_scenario() |>
  oncology_scenario_set_truevalues()

#pre_sim_parameters <- oncology_scenario()

pre_N_sim <- 2

pre_my_analyse <- list(
  truth = function(condition, dat, fixed_objects = NULL) {
    mod <- survival::coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = dat)
    HR <- exp(coef(mod)[1])
    list(HR = HR)
  }
)

pre_my_summarise <- create_summarise_function(
  truth = summarise_estimator( #test if this can be an object defined before
    est = HR,
    real = exp(beta_death[[1]][6]), #actually not required, but we need a placeholder here
    null = 1
  )
)

summy<-summarise_estimator( #test if this can be an object defined before
  est = HR,
  real = exp(beta_death[[1]][6]), #actually not required, but we need a placeholder here
  null = 1
)
pre_my_summarise <- create_summarise_function(
  truth = summy
)

cl <- makeCluster(detectCores(logical = FALSE) - 1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
  library("survival")
})

#SimClean()

pre_results <- runSimulation(
  design = sim_parameters,
  replications = pre_N_sim,
  generate = generate_oncology,
  analyse = pre_my_analyse,
  summarise = pre_my_summarise,
  fixed_objects = list(allow_switch = FALSE, logHR_assumed = NULL, ev_soll = 100, allow_random_cens = TRUE),
  parallel = TRUE,
  cl = cl
)
as.data.frame(pre_results)
stopCluster(cl)

# Under H0 the true effect is HR = 1
pre_results[which(sapply(pre_results$beta_death, `[[`, 6) == 0),]$truth.mean_est <- 1

# Define parameter values and derived quantities -------------------------

sim_parameters <- oncology_scenario() |>
  oncology_scenario_set_truevalues() |>
  dplyr::mutate(true_eff = pre_results$truth.mean_est) #redundant, pool in first step

sim_parameters <- sim_parameters |>
  dplyr::mutate(true_eff = pre_results$truth.mean_est) #redundant, pool in first step

# Constants for simulation -----------------------------------------------

N_sim <- 10
alpha <- 0.05

# List of analysis functions ---------------------------------------------

my_analyse <- list(
  rpsftm_rc = analyse_oncology_rpsftm(recensor = TRUE),
  rpsftm = analyse_oncology_rpsftm(recensor = FALSE),
  tse_rc = analyse_oncology_TSE(recensor = TRUE),
  tse = analyse_oncology_TSE(recensor = FALSE),
  gformula = analyse_oncology_gformula(B = 20),
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


summy<-summarise_estimator(
  est = HR,
  real = true_eff,
  lower = low,
  upper = up,
  null = 1,
  name = "est"
)

sumtest<-summarise_test(
  alpha/2,
  name = "test"
)
my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  rpsftm_rc = summy,
  rpsftm = summy,
  tse_rc = summy,
  tse = summy,
  gformula = summy,
  ipw = summy,
  itt = summy,
  cens = summy,
  # rejection rates
  rpsftm_rc = sumtest,
  rpsftm = sumtest,
  tse_rc = sumtest,
  tse = sumtest,
  gformula = sumtest,
  ipw = sumtest,
  itt = sumtest,
  cens = sumtest,
  describe = summarise_describe()
)

# Run the simulations ----------------------------------------------------

cl <- makeCluster(detectCores(logical = FALSE) - 1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})
clusterExport(cl = cl, varlist = c("alpha"))

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

results <- runSimulation(
  design = sim_parameters[1:3,],
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
path="data/"
save(results, main_sessioninfo, nodes_sessioninfo, file = paste(path,format(Sys.time(), paste0("results_onco_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")),sep=""))
