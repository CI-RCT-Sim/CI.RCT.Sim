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

cl <- makeCluster(detectCores(logical=FALSE)-1)
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

# Set true value = 0 under the null hypothesis
pre_results[which(pre_results$eff_true == 0),]$tp_mean.mean_est <- 0

# Define parameter values and derived quantities -------------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues() |>
  dplyr::mutate(tp_eff = pre_results$tp_mean.mean_est)

# Constants for simulation -----------------------------------------------

N_sim <- 10000
alpha <- 0.025

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
  mihyp = analyse_diabetes_mi(strategy = "hypothetical"),
  describe = function(condition, dat, fixed_objects = NULL) {
    result <- list(
      n_pat  = nrow(dat),
      n_ctrl = sum(dat$trt == 0),
      n_trt  = sum(dat$trt == 1),
      n_resc = sum(dat$rescue_start < 14, na.rm = TRUE),
      n_miss = sum(dat$m_start < 13, na.rm = TRUE),
      naive  = mean(dat[which(dat$trt == 1),]$y12-dat[which(dat$trt == 1),]$y0, na.rm = TRUE) -
        mean(dat[which(dat$trt == 0),]$y12-dat[which(dat$trt == 0),]$y0, na.rm = TRUE)
    )
    result
  }

)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

sum_tp <- summarise_estimator(est = coef, real = tp_eff, lower = ci_lower,
  upper = ci_upper, null = 0, name = "est")

sum_hyp <- summarise_estimator(est = coef, real = eff_true, lower = ci_lower,
  upper = ci_upper, null = 0, name = "est")

sum_test <- summarise_test(alpha, name = "test")


my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  ## Treatment policy estimands
  ipwtp = sum_tp,
  mmrmtp = sum_tp,
  mitp = sum_tp,
  ## Hypothetical estimands
  ipwhyp = sum_hyp,
  dmhyp = sum_hyp,
  gcomhyp = sum_hyp,
  mmrmhyp = sum_hyp,
  mihyp = sum_hyp,
  # rejection rates
  ## Treatment policy estimands
  ipwtp = sum_test,
  mmrmtp = sum_test,
  mitp = sum_test,
  ## Hypothetical testimands
  ipwhyp = sum_test,
  dm = sum_test,
  gcom = sum_test,
  mmrmhyp = sum_test,
  mihyp = sum_test,
  mmrmtp = function(condition, results, fixed_objects = NULL){
    data.frame(n_conv = sum(results$converged),
               fallbacks = sum(results$fallback),
               unstruc = sum(results$covariance == "us"))},
  mmrmhyp = function(condition, results, fixed_objects = NULL){
    data.frame(n_conv = sum(results$converged),
               fallbacks = sum(results$fallback),
               unstruc = sum(results$covariance == "us"))},
  describe = summarise_describe()
)

# Run the simulations ----------------------------------------------------

cl <- makeCluster(detectCores(logical=FALSE) - 1)
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
