# script to calculate the empirical true treatment effect

# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

# Derive true treatment effect  under treatment policy -------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues() |>
  dplyr::mutate(nfix = 1e6, miss = rep(list(c(-1e5, 0, 0, 0)), 14))

N_sim <- 40

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

cl <- makeCluster(detectCores(logical=FALSE)-1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
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

# Set true value = 0 under the null hypothesis
results[which(results$eff_true == 0),]$tp_mean.mean_est <- 0

# Save results -----------------------------------------------------------

saveRDS(results, file=format(Sys.time(), paste0("tp_eff_%m-%d.rds")))
