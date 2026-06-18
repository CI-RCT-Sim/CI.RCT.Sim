# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

source("scripts/vaccine_scenario_classes.R")
source("scripts/vaccine_analyse_summarise.R")

# Define parameter values and derived quantities -------------------------

scenario <- Sys.getenv("scenario")

message(paste("Scenario", scenario, "selected"))

selected_scenario <- switch(
  scenario,
  A1 = vaccine_scenario_A1,
  A2 = vaccine_scenario_A2,
  B1 = vaccine_scenario_B1,
  C1 = vaccine_scenario_C1,
  D1 = vaccine_scenario_D1,
  extra = vaccine_scenario_extra,
  stop("unknown scenario selected or selection missing")
)


# get row number to run ---------------------------------------------------

row <- Sys.getenv("row") |>
  strtoi()

# run simulations ---------------------------------------------------------

message(paste0("Running row ", row, " of ", nrow(selected_scenario), " of scenario ", scenario))

sim_parameters <- selected_scenario |>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize() |>
  within({
    VE = 1-rr_ps
    scenario_nr = seq_along(VE)
  }) |>
  _[row, ]

# Constants for simulation -----------------------------------------------

# N_sim <- 10
N_sim <- 10000
alpha_ci <- 0.05
alpha_test <- c(0.05, 0.025)
N_cores <- 124

# get nodename ------------------------------------------------------------

if(Sys.getenv("SLURM_JOB_ID") == ""){
  my_nodename <- Sys.info()["nodename"]
} else {
  my_nodename <- paste0("SlurmJob_", Sys.getenv("SLURM_JOB_ID"))
}

# Run the simulations ----------------------------------------------------

message(paste(length(my_analyse), "analysis functions"))
message(paste(length(environment(my_summarise)$summarise_functions), "summarise functions"))
message(paste("setting up cluster with", N_cores , "cores"))

cl <- makeCluster(N_cores)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
})

clusterExport(cl = cl, varlist = c("alpha_ci", "alpha_test"))

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

message(Sys.time())
message(paste("running simulations,", nrow(sim_parameters), "scenarios,", N_sim, "replications"))

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_vaccine,
  analyse = my_analyse,
  summarise = my_summarise,
  fixed_objects = list(include_unobserved=FALSE),
  parallel = TRUE,
  cl = cl,
  save_details = list(
    compname = my_nodename
  )
)

message(Sys.time())
message("stopping cluster")
stopCluster(cl)

# Save results -----------------------------------------------------------

out_filename <- format(Sys.time(), paste0("results_vaccine_scenario_", scenario, "_row_", formatC(row, width=3, flag="0"), "_", my_nodename, "_%Y-%m-%d_%H%M.Rdata"))
message(paste("saving results to", out_filename))
save(results, main_sessioninfo, nodes_sessioninfo, file=out_filename)


