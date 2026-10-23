library(CI.RCT.Sim)
library(parallel)

# Load the analysis and summarise functions
source("scripts/diabetes_analyse_summarise.R")

# Load the pre-calculated treatment policy effect estimates, originates from the script true_effect.R
tp_eff <- readRDS("tp_eff_10-23.rds")$tp_mean.mean_est

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues() |>
  dplyr::mutate(tp_eff = tp_eff)

# Constants for simulation -----------------------------------------------
N_sim <- 10000
alpha <- 0.025
N_cores <- detectCores(logical=FALSE)-1

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

clusterExport(cl = cl, varlist = c("alpha"))

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

message(Sys.time())
message(paste("running simulations,", nrow(sim_parameters), "scenarios,", N_sim, "replications"))

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_diabetes,
  analyse = my_analyse,
  summarise = my_summarise,
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

out_filename <- format(Sys.time(), paste0("results_diabetes_row_", formatC(row, width=3, flag="0"), "_", my_nodename, "_%Y-%m-%d_%H%M.Rdata"))
message(paste("saving results to", out_filename))
save(results, main_sessioninfo, nodes_sessioninfo, file=out_filename)


